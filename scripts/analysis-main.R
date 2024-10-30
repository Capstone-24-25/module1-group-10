##------- SET UP and IMPORT DATA -----------------------------------------------------------
library(ggplot2)
library(tidyverse)
setwd("~/Downloads/module1-group-10/data")
biomarker_raw<-read.csv("biomarker-raw.csv")
load('biomarker-clean.RData')

##------- Question 1 -----------------------------------------------------------

### Sampled 5 proteins from the raw dataset and isolated their protein levels
set.seed(1045)
protein_numbers <- sample(3:ncol(biomarker_raw), size = 5, replace=F)
protein_sample <- biomarker_raw[,protein_numbers]
colnames(protein_sample) <- protein_sample[1,]
protein_sample <- protein_sample[-1,]

### Distributions of the protein levels for these 5 proteins is shown below
protein_sample2 <- protein_sample %>% pivot_longer(cols = everything(), 
                                                   names_to = "protein", 
                                                   values_to = "level")
protein_sample2$level <- as.numeric(protein_sample2$level)
protein_sample2 %>% ggplot() + geom_histogram(aes(x = level), bins = 50) + 
  facet_wrap("protein") + xlim(0,8000)

### QQPlots of the protein levels for these 5 proteins is shown below
protein_sample2 %>% ggplot(aes(sample = level)) + stat_qq() + 
  stat_qq_line(color = "red") + facet_wrap("protein") + ylim(0, 10000)

### Shapiro tests for sample of 5 proteins
protein_sample2 %>% group_by(protein) %>% summarise("pval" = 
  shapiro.test(as.numeric(level))$p.value) %>% 
  mutate("normal" = ifelse(pval < 0.05, FALSE, TRUE))

### Import data from biomarker_clean to represent log transformed data 
### for the sampled proteins
log_protein_sample <- biomarker_clean %>% select("PGRP-S", "PACAP-27", "TRAIL R4", "IGFBP-1", "cIAP-2") %>% pivot_longer(cols = everything(), names_to = "protein", values_to = "level")
log_protein_sample$level <- as.numeric(log_protein_sample$level)

### Histograms of Log Transformed Data
log_protein_sample %>% ggplot() + geom_histogram(aes(x = level), bins = 50) + facet_wrap("protein")

### QQplots of Log Transformed Data
log_protein_sample %>% ggplot(aes(sample = level)) + stat_qq() + stat_qq_line(color = "red") + facet_wrap("protein")


##------- Question 2 -----------------------------------------------------------
# Comparing values to scaled values

# get names
var_names <- read_csv('biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()


biomarker_scaled <- read_csv('biomarker-raw.csv', 
                             skip = 2,
                             col_select = -2L,
                             col_names = c('group', 
                                           'empty',
                                           pull(var_names, abbreviation),
                                           'ados'),
                             na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale without trim
  mutate(across(.cols = -c(group, ados), 
                ~scale(log10(.x))[, 1])) %>%
  # reorder columns
  select(group, ados, everything())

# Biomarker scaled data without group and ados columns

biomarker_outlier <- biomarker_scaled %>%
  select(-c(group, ados))


### Counting the number of outliers for each group

# Number of outliers per person are all data points that equal 3
above_three_outliers <- rowSums(biomarker_outlier > 3)
# or negative 3
below_neg_three_outliers <- rowSums(biomarker_outlier < -3)

# The total number of outliers per person is the sum of the two
outlierspp <- above_three_outliers + below_neg_three_outliers
outlierspp

# adding a column to biomarker_clean with number of outliers per person
biomarker_w_outliers <- biomarker_clean
biomarker_w_outliers$outliers <- outlierspp

# Grouping the dataset in to ASD and TD and comparing the number of outliers
biomarker_w_outliers %>%
  group_by(group) %>%
  summarize(TotalOutliers = sum(outliers),
            OutliersMean = mean(outliers),
            OutilersSD = sd(outliers),
            OutliersMax = max(outliers),
            OutlierNumbers = sum(outliers > 0),
            NonOutlierNumbers = sum(outliers == 0))

# There are 1007 outliers for people with ASD and 1372 for TD children

# export as r binary
save(list = 'biomarker_w_outliers', 
     file = 'data/biomarker_w_outliers.RData')
##------- Question 3 -----------------------------------------------------------

# repeat the analysis but carry out the entire selection procedure on a training partition
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

var_names <- read_csv("biomarker-raw.csv", col_names = F, n_max = 2, col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, abbreviation = V2) %>%
  na.omit()

biomarker_data <- read_csv("biomarker-raw.csv", skip = 2, col_select = -2L,
                           col_names = c('group', 'empty', pull(var_names, abbreviation), 'ados'),
                           na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  mutate(across(.cols = -c(group, ados), ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  select(group, ados, everything())

set.seed(101422)
biomarker_split <- initial_split(biomarker_data, prop = 0.8)
biomarker_train <- training(biomarker_split)
biomarker_test <- testing(biomarker_split)

# choose a larger number (more than ten) of top predictive proteins using each selection method
# used the top 20 proteins based on their significance or importance scores (for both t-test and random forest)
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)
predictors <- biomarker_train %>% select(-c(group, ados))
response <- biomarker_train %>% pull(group) %>% factor()

set.seed(101422)
rf_out <- randomForest(x = predictors, y = response, ntree = 1000, importance = TRUE)

proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)

# use a fuzzy intersection instead of a hard intersection to combine the sets of top predictive proteins across selection methods
protein_union <- union(proteins_s1, proteins_s2)
ranked_proteins <- tibble(
  protein = protein_union,
  rank_s1 = match(protein_union, proteins_s1, nomatch = 21),
  rank_s2 = match(protein_union, proteins_s2, nomatch = 21)
) %>%
  mutate(total_rank = rowMeans(across(starts_with("rank")), na.rm = TRUE)) %>%
  arrange(total_rank) %>%
  slice_min(total_rank, n = 20) %>%
  pull(protein)

##------- Question 4 -----------------------------------------------------------
# get names
var_names <- read_csv('biomarker-raw.csv', 
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
biomarker_clean <- read_csv('biomarker-raw.csv', 
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