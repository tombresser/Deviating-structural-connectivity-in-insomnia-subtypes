rm(list=ls())

### set enviroment ####
# specify paths and files
path.project <- '/localpath'
path.data <- paste0(path.project, "/data")
path.nkeviewer <- "/localpath/nke-viewer"
groundhog.day <-'2022-02-12'     #yyyy-mm-dd

file.blanken <- 'Blanken_2019_charasteristics_regions_manual.csv'

# libraries
pkgs <- c("tidyverse", "openxlsx", "hrbrthemes")
groundhog.library(pkgs, groundhog.day)

# working dir
setwd(path.project)


### load data source ####
# Blanken et al data
df.blanken <- read.csv(file = paste0(path.data, "/", file.blanken))
rm(list = ls(pattern = "file.*"))

# cleaning
df.blanken <- df.blanken %>%
  mutate(across(.cols = everything(), tolower))

# nke-viewer data
# find all words.csv files
files.csv <- list.files(path = paste0(path.nkeviewer, "/data"),  pattern = ".csv", recursive = TRUE)

# load first file
df.nke <- read.csv(file = paste0(path.nkeviewer, "/data/", files.csv[1]))

# add cluster grouping columns
df.nke <- df.nke %>%
  mutate(cluster_grouping = dirname(files.csv[1]))

# now append all other files row wise
for (n in seq(2, length(files.csv))) {
  print(files.csv[n])
  # load data
  TEMP <- read.csv(file = paste0(path.nkeviewer, "/data/", files.csv[n]))
  
  # add cluster grouping column
  TEMP <- TEMP %>%
    mutate(cluster_grouping = dirname(files.csv[n]))
  
  # add rowwise to existing df
  df.nke <- bind_rows(df.nke, TEMP)
  rm(TEMP)
}

# make sure everything is lower case
df.nke <- df.nke %>%
  mutate(across(.cols = everything(), tolower))


### match insomnia terms to NKE ####
# check for traits with spaces
grep(" ", df.blanken$charasteristic, value = T)

# replace spaces in insomnia terms
terms.blanken <- df.blanken %>%
  mutate(charasteristic = gsub(" ", "_", charasteristic)) %>%
  select(charasteristic) %>%
  pull()

# expand insomnia terms based on subterms, synonyms and expert knowledge (prof.dr. EJW van Someren)
terms.expanded<- c(terms.blanken, "insomnia", "life_events", "stress", "pleasure", "arousal", "mood", "happiness", "salience", 
                   "emotional_memory", "reward", "cognitive", "cognition", "interoception", "loss", "response_inhibition",
                   "negative_emotion", "emotion", "trauma", "gonogo_task", "sadness", "affect", "reflection", "reflect")

# store in list
insomnia.terms <- list(blanken = terms.blanken,
                       expanded = terms.expanded)

# save terms and clean
write_csv(as_tibble(insomnia.terms$expanded), file = "data/Blanken_search_terms.csv")
rm(terms.blanken, terms.expanded)

# check NKE terms
df.nke %>%
  distinct(., TOKEN) 

# create grep pattern from search terms
search.traits  <- paste(insomnia.terms$expanded, collapse = "|")

# label matches as id_trait
df.nke.filtered <- df.nke %>%
  mutate(id_trait = ifelse(grepl(search.traits, TOKEN) == T, 1, 0))

# check number of matches per cluster grouping
df.nke.filtered %>%
  group_by(cluster_grouping) %>%
  summarise(n=sum(id_trait)) %>%
  pivot_wider(names_from = cluster_grouping, values_from = n)


### correct false matches ####
# manual check of labelled id_traits
df.nke.filtered %>%
  filter(id_trait == 1) %>%
  distinct(., TOKEN, .keep_all = T) %>%
  select(TOKEN)

# define wrong matches
wrong <- c("sexual_arousal", "orgasm", "visual_recognition", "word_recognition", "visual_word_recognition",
           "face_recognition", "object_recognition", "recognition")

# document excluded NKE matches
write_csv(as_tibble(wrong), file = "data/NKE_excluded_terms.csv")

# remove them from results
df.nke.filtered <- df.nke.filtered %>% 
  mutate(id_trait = ifelse(TOKEN %in% wrong, 0, id_trait))
rm(wrong)

# document all NKE tokens (terms) that are a "correct" match with Blanken
df.nke.filtered %>%
  filter(id_trait == 1) %>%
  distinct(., TOKEN, .keep_all = T) %>%
  select(TOKEN) %>%
  write_csv(., file = "data/NKE_matched_terms.csv")


### label insomnia domains (>= 40% id_traits) ####
overview <- df.nke.filtered %>%
  group_by(cluster_grouping, CLUSTER, DOMAIN) %>%
  # determine number of insomnia traist within domain
  summarise(id_traits = sum(id_trait),
            # total number of traits within a domain
            traits = n(),
            # perentage traits linked to id
            id_overlap = id_traits / traits
            ) %>%
  # only keep domains with at least 40% insomnia traits
  filter(id_overlap >= 0.4) %>%
  arrange(cluster_grouping) %>%
  # create k-specific domain ID
  mutate(grouping_domain = paste0(cluster_grouping, "_", DOMAIN))

# label insomnia domains in df.nke.filtered
df.nke.filtered <- df.nke.filtered %>%
  # create similar k-specific domain ID as in "overview"
  mutate(grouping_domain = paste0(cluster_grouping, "_", DOMAIN),
    # label "insomnia domains" using k-specific domain IDs (to prevent labeling domains with the same name in a different k)
    id_domain = ifelse(grouping_domain %in% overview$grouping_domain, 1, 0))

# plot insomnia domains across k's
plot.domains <- df.nke.filtered %>%
  select(cluster_grouping, DOMAIN, id_domain) %>%
  mutate(id_domain = as.factor(id_domain)) %>%
  ggplot(aes(y = reorder(DOMAIN, id_domain, function(x){ sum(as.numeric(x)) }), 
             x = cluster_grouping, fill= id_domain)) +
  geom_tile() +
  scale_fill_manual(values=c("#FFFFFF", "#000000")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(title = "NKE domains labelled as insomnia related per k", y = "NKE domain", 
       x = "k", caption = "(Number of domains defined in the brain)")
jpeg(filename = "figures/NKE_insomnia_domains.jpeg")
plot.domains
dev.off()


### add domains with <40% insomnia traits but important traits
# some traits are important but unlikely to reach 40% within a domain
# response inhibition, introception, gonogo_task
extra_domains <- paste(c("response_inhibition", "interoception", "gonogo_task"), collapse = "|")
df.nke.filtered <- df.nke.filtered %>%
  mutate(id_domain = ifelse(id_trait == 1 & grepl(extra_domains, TOKEN) == T, 1, id_domain)) 


# make file with cluster_grouping + insomnia CLUSTER (number associated with insomnia domain)
df.nke.filtered %>%
  filter(id_domain == 1) %>%
  #focus on higher cluster groupings, smaller number could be less specific
  filter(cluster_grouping >= "k25") %>%
  select(cluster_grouping, CLUSTER, DOMAIN) %>%
  group_by(cluster_grouping) %>%
  summarise(clusters = toString(unique(CLUSTER)),
            domains = toString(unique(DOMAIN))) %>%
  write_csv(., "data/insomnia_clusters_k25-50.csv")
