# Comparing values to scaled values

# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()


biomarker_scaled <- read_csv('data/biomarker-raw.csv', 
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