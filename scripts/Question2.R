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

### Counting the number of outliers for each group

# Number of outliers per person are all data points that equal 3
above_three_outliers <- rowSums(biomarker_clean == 3)
# or negative 3
below_neg_three_outliers <- rowSums(biomarker_clean == -3)

# The total number of outliers per person is the sum of the two
outlierspp <- above_three_outliers + below_neg_three_outliers
outlierspp

# turning the NAs to zeros
outlierspp[is.na(outlierspp)] <- 0
outlierspp

# adding a column to biomarker_clean with number of outliers per person
biomarker_w_outliers <- biomarker_clean
biomarker_w_outliers$outliers <- outlierspp

# Grouping the dataset in to ASD and TD and comparing the number of outliers
biomarker_w_outliers %>%
  group_by(group) %>%
  summarize(TotalOutliers = sum(outliers))

# There are 1007 outliers for people with ASD and zero for TD children