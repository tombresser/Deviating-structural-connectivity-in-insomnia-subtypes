rm(list=ls())

# libraries
library("tidyverse")
library("kableExtra")
require("webshot2")
library("flextable")
library("magrittr")

# specify paths and files
path.project <- '/local_path'
file_demo <- paste0(path.project, "/ITQ_demographics.csv")

# set wroking dir
setwd(path.project)

# load data
data_table <- read.csv(file_demo)

## Build table ----
library("table1")

# first look
table1(~ Sex + Age + ISI_score + PSQI_score_PP + sleep_problem_duration + IDSSR_score +
         HADS_depression + HADS_anxiety + BAI_score | factor(insomnia_subtype), 
       data=data_table,
       overall=F)

# ------------------------------
# specify some field and functions to improve table
# total number of people with insomnia
n_insomnia <- data_table %>% 
  filter(insomnia_subtype != "Control group") %>% 
  nrow()

# custom labels
labels <- list(
  variables=list(Sex="Sex",
                 Age="Age (years)",
                 ISI_score="Insomnia severity (ISI)",
                 PSQI_score_PP="Sleep quality (PSQI)",
                 sleep_problem_duration="Insomnia duration (years) ",
                 IDSSR_score="Depression severity (IDS-SR)",
                 HADS_depression="Depression severity (HADS)",
                 HADS_anxiety="Anxiety (HADS)"),
  groups=list( "", c(paste0("Insomnia disorder (N=", n_insomnia, ")"))))

# render functions to specify how to display continuous and categorical variables
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

# split data by group(insomnia subtype)
strata <- split(data_table, data_table$insomnia_subtype)

# create table using specified fields and render functions
table1(strata, labels, groupspan=c(1,5),
       render.continuous=my.render.cont, render.categorical=my.render.cat)  
  # export or copy + special paste as html
