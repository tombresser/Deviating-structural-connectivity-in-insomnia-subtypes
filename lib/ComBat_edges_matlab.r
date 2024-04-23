rm(list=ls())

# libraries
library("tidyverse")
library("neuroCombat")

# data path
log_path <- '/some/path/log'
path.cache <- '/some/path/data'


# -----------------------------------------------------------------------
### log functions because this script is executed from within matlab 
# Create log file name
get_log_file_name <- function(file_name="matlab_combat")
{
  file_name <- paste(file_name,format(Sys.time(), "%Y%m%d"),sep="_")
  file_name <- paste(log_path, paste(file_name,"log", sep="."), sep="/")
  return(file_name)
}

# Assigning log file into a logger variable
logger <- file(get_log_file_name(), open = "at")

# This for redirecting any warning or error outputs
sink(logger, type="message")

# logging function
logit <- function(msg, ...) {
  cat(format(Sys.time(), "%Y-%m-%d %X"), ":", paste(msg, ...), "\n", append = T,
      file = logger)
  }

logit("Function started")


# -----------------------------------------------------------------------
# load data
df.dat <- read.csv(paste0(path.cache, '/tmp_combat_data.csv'),
                  header = FALSE)
df.mod <- read.csv(paste0(path.cache, '/tmp_combat_model.csv'),
                   header = TRUE)

# double check loaded data
mod_check <- colnames(df.mod) %in% c("subtype", "age", "male", "BrainSegVolNotVent", "study")

if (all(mod_check)) {
  # scanner/study vector
  batch <- df.mod %>% 
    pull(study)
  
  # covarites matrix
  mod <- df.mod %>%
    rename("sex" = male) %>%
    mutate(
      sex = recode(sex, `1` = "male", `0` = "female"),
      subtype = as_factor(subtype)) %>%
    model.matrix(~ subtype + age + sex + BrainSegVolNotVent, .)
  
  # in edge data set zeros to NA
  dat <- df.dat %>%
    mutate(across(where(is.numeric), 
                  ~na_if(., 0)))
  logit("Convert zeros to NA in the data")
  
  # exclude edges with 0 or 1 value within a batch
  # ComBat is unable to handle this
  batch     <- as.factor(batch)
  
  batches   <- lapply(levels(batch), function(x)which(batch==x))
  n.batches <- sapply(batches, length)
  nas.batches <- lapply(batches, function(x) colSums(is.na(dat[x,])))
  logic.batches <- map2(n.batches, nas.batches, ~ (.x - .y) <= 2)
  
  logic.edge <- logic.batches %>%
    bind_rows() %>%    # unlist, row per batch
    colSums() == 0    # all bacthes should have at least 1 value in an edge
  
  # index vector of edges suitable for neuroCombat
  good.edges <- as.vector(logic.edge)
  logit(paste0("Harmonizing data using", sum(good.edges), " edges"))
  
  # select edges with sufficient values within each batch
  new.dat <- dat[, good.edges]
  
  
  # columns should be subjects
  if (length(batch) != ncol(new.dat)) {
    # if not, pivot dat
    ID <- paste0("S", seq(1, nrow(new.dat)))
    new.dat <- new.dat %>%
      add_column(ID) %>%                      # add col with subj ID
      pivot_longer(-ID, names_to = "tmp") %>% # edge name to tmp
      pivot_wider(names_from = ID) %>%        # subj ID as cols
      select(!tmp)                            # drop tmp col
    logit("Pivotted data to columns = subjects")
  }
  
  # harmonize data using ComBat
  df.harmonized <- neuroCombat(dat=new.dat, batch=batch, mod=mod, ref.batch = "ercp")

  # transform  harmonized data back to all edges (NA if not in ComBat)
  out <- matrix(NA, nrow = ncol(dat), ncol = nrow(dat))
  out[good.edges,] <- df.harmonized$dat.combat
  
  # rotate back rows = subjects, columns = edges
  if (length(batch) == ncol(out)) {
    # if not, pivot dat
    edge <- paste0("e", seq(1, nrow(out)))
    colnames(out) <- paste0("S", seq(1, ncol(out)))
    out <- out %>%
      as_tibble() %>%
      add_column(edge) %>%                      # add col with edge ID
      pivot_longer(-edge, names_to = "subj") %>% # edge name to tmp
      pivot_wider(names_from = edge) %>%        # edge ID as cols
      select(!subj)                            # drop tmp col
    logit("Pivotted data back to columns = edge and rows = subjects")
  }
  
  # save output
  write_csv(out, file = paste0(path.cache, "/tmp_combat_harmonized.csv"),
            col_names = FALSE)
  
  logit("Data harmonized and saved as", paste0(path.cache, '/tmp_combat_harmonized.csv'))
  
  # prior images
  FIG <- str_replace(path.cache, "cache", "figures")
  png(paste0(FIG, "/", format(Sys.time(), "%Y%m%d"),"_ComBar_delta_prior.png"))
  drawPriorDelta(df.harmonized$estimates, xlim=c(0,2.5))
  dev.off()
  
  FIG <- str_replace(path.cache, "cache", "figures")
  png(paste0(FIG, "/", format(Sys.time(), "%Y%m%d"),"_ComBar_gamma_prior.png"))
  drawPriorGamma(df.harmonized$estimates, xlim=c(-1.5,1.5))
  dev.off()
  
} else {
  logit("Column names of ", paste0(path.cache, '/tmp_combat_model.csv', "did not match the expected covariates"))
}

