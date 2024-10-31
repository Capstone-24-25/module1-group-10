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

#Question 4
set.seed(101422)

biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.8)
biomarker_train <- training(biomarker_split)
biomarker_test <- testing(biomarker_split)

test_fn <- function(.df){
  t_test(.df, formula = level ~ group, order = c('ASD', 'TD'), alternative = 'two-sided', var.equal = F)
}

ttests_out <- biomarker_train %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = 'protein', values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(), hm = log(m) + 1/(2*m) - digamma(1), rank = row_number(), p.adj = m*hm*p_value/rank)

proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)

predictors <- biomarker_train %>%
  select(-c(group, ados))
response <- biomarker_train %>% pull(group) %>% factor()

rf_out <- randomForest(x = predictors, y = response, ntree = 1000, importance = T)

proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)
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

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(ranked_proteins)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# Define the full model with all predictors
initial_model <- glm(class ~ 1, data = train_data, family = binomial)
full_model <- glm(class ~ ., data = train_data, family = binomial)

forward_model <- stats::step(initial_model, direction='forward', scope = list(lower = initial_model, upper = full_model))

backward_model <- glm(class ~ ., data = train_data, family = binomial) %>% stats::step(direction = "backward")

#Stepwise regression
both_model <- glm(class ~ ., data = train_data, family = binomial)%>% stats::step(direction = "both", scope = list(lower = initial_model, upper = full_model))

summary(forward_model)
summary(backward_model)
summary(both_model)

par(mfrow = c(2, 2))
plot(forward_model$residuals, main = "Forward Residuals", ylab = "Residuals")
plot(backward_model$residuals, main = "Backward Residuals", ylab = "Residuals")
plot(both_model$residuals, main = "Both-Direction Residuals", ylab = "Residuals")

# forward: AIC=116.67
# class ~ DERM + 14-3-3 protein zeta/delta + IgD + PTN + M2-PK

# backward: AIC=116.24
# class ~ DERM + IgD + PTN + RELT + CK-MB + 14-3-3 protein zeta/delta

# both: AIC=116.24
# class ~ DERM + IgD + PTN + RELT + CK-MB + 14-3-3 protein zeta/delta

# To find the best panel of predictive proteins, one method we used was step wise regression. 
# Using the chosen 20 predictive proteins, forwards, backwards, and both direction step wise 
# regression was applied. The goal was to find the model with the lowest possible AIC score 
# as it would mean it would be best in describing variance in the data. 
# The forward stepwise regression achieved a AIC of 116.67 while the backwards and 
# both direction methods generated a AIC score of 116.24 and had the same protein predictors. 
# Thus we will go with the proteins that appeared in 2 or more of the methods and had 
# achieved the lowest AIC. 6 of the 20 proteins were selected, DERM, IgD, PTN, RELT, 
# CK-MB, 14-3-3 protein zelta/delta, which we will now identify as our core proteins. 
# Looking at the metrics, the null deviance was 170.44 and residual deviance was 102.24 
# meaning that our proteins were really helpful in predicting whether the participant was 
# in the ASD or TD group.

# Looking at the residual graphs, the residuals are randomly scattered around zero with no 
# clear pattern which tells us there is constant variance and confirms the strength of the 
# model. 