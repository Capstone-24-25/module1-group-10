library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

setwd("~/Pstat 197/module1-group-10")

# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2),
                      show_col_types = FALSE) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_clean <- read_csv('data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            show_col_types = FALSE,
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ scale(log10(.x))[, 1])) %>%
  # reorder columns
  select(group, ados, everything())


#exploratory analysis of outlying values
outliers <- biomarker_clean %>%
  rowwise() %>%
  mutate(across(.cols = c(-group,ados),
                ~ abs(.x) > 3)) %>%
  mutate(num_outlier = sum(c_across(-c(group, ados)))) %>%
  select(group, num_outlier)

summary <- outliers %>%
  group_by(group) %>%
  summarize(outlier_mean = mean(num_outlier),
            outlier_sd = sd(num_outlier),
            outlier_max = max(num_outlier),
            outlier_total = sum(num_outlier),
            outlier_subjects = sum(num_outlier > 0))

outliers
summary

# export as r binary
save(list = 'biomarker_clean', 
     file = 'data/biomarker-clean.RData')

