## code summary:
# calculate F1 of predicted donor-recipient pairings vs true pairings for DNA distance 0.001-3
# select optimal DNA distance threshold for strain matching
# plot AUCROC curve for two models used for FOCUS trial
# plot AUCROC curve for FOCUS trial and Gut Bugs Trial models
# plot precision recall curve for FOCUS trial and Gut Bugs Trial models

## load libraries

library(readxl)
library(tidyverse)
library(caret)
library(openxlsx)

## load data

load("uc_metadata_subset.RData") # uc_metadata

load("strainphlan.Rdata") # strainphlan output

uc_pairings <- read_excel(path = "FOCUS Patient Donor matching metadata.xlsx", skip = 1, col_names = TRUE) # true pairings

## load functions

# function to run strain engraftment analysis with different thresholds
strain_engraftment_analysis <- function(threshold, strain_dna_dist){
  
  # produce strain matches with threshold
  strain_matches <- strain_dna_dist %>% 
    filter(Sample1_details != Sample2_details, # remove same sample comparisons
           Dist_norm <= threshold) # apply genetic threshold for strain match
  
  # generate all strain matches, using only individual donors (not donor batches)
  strain_matches_individual_donors <- strain_matches %>%
    filter(!grepl("Donor_batch", Sample1_details),    # Sample 1 = post-FMT recipient sample
           !grepl("Donor_batch", Sample2_details))  # Sample 2 = donor sample
  
  # find matching strains present in recipient at baseline and week 8
  FMT_BL_wk8_donor_sample_matches <- strain_matches_individual_donors %>% 
    filter(grepl("FMT_wk8|FMT_BL", Sample1_details),    # Sample 1 = recipient samples
           grepl("Donor", Sample2_details)) %>%  # Sample 2 = donor sample
    mutate(Donor_ID = substr(Sample2_details,1,4), # extract donor ID
           Recipient_ID = substr(Sample1_details,1,4), # extract recipient ID
           Timepoint = str_split(Sample1_details, "_", simplify = TRUE)[,3]) %>% # extract timepoint
    select(c(Donor_ID, Sample2, Recipient_ID, Species, Timepoint)) %>% # select relevant columns
    rename(Donor_sample = Sample2) %>% # rename donor sample column
    spread(key = Timepoint, value = Timepoint) %>% # spread the timepoint column
    replace_na(list(BL = "absent", wk8 = "absent")) %>% # replace NAs as strain absent
    mutate(BL = str_replace(BL, "BL", "donor_match"),
           wk8 = str_replace(wk8, "wk8", "donor_match")) %>%
    arrange(Donor_ID, Donor_sample, Recipient_ID, Species) %>%
    group_by(Donor_ID, Recipient_ID, Species) %>%
    filter(any(BL == "donor_match") & any(wk8 == "donor_match")) %>%
    ungroup() %>%
    arrange(Donor_ID, Recipient_ID, Donor_sample, Species)
  
  # baseline-matching strains to subtract from data
  FMT_BL_subtraction <- FMT_BL_wk8_donor_sample_matches %>%
    select(Donor_ID, Recipient_ID, Species) %>%
    distinct()
  
  # subtract matching strains also at baseline
  FMT_engraftment_individual_donors_subBL <- strain_matches_individual_donors %>% 
    filter(grepl("FMT_wk8", Sample1_details),    # Sample 1 = post-FMT recipient sample
           grepl("Donor", Sample2_details)) %>%  # Sample 2 = donor sample
    mutate(Donor_ID = substr(Sample2_details,1,4), # extract donor ID
           Recipient_ID = substr(Sample1_details,1,4)) %>%  # extract recipient ID
    anti_join(FMT_BL_subtraction, by = c("Donor_ID", "Recipient_ID", "Species")) # FMT week 8 donor matches not also in recipient at BL
  
  FMT_engraftment_individual_donors_subBL # return result
}

# function to run strain engraftment analysis with different thresholds, adjusting for donor-matching strains in placebos
strain_engraftment_analysis_placebo_adjusted <- function(threshold, strain_dna_dist, sea_output){
  # sea_output: output from baseline-adjusted strain engraftment analysis
  # this function performs strain engraftment analysis (baseline-adjusted) to find donor matches in placebos at week 8
  # those placebo matches are then subtracted from the original sea_output results
  
  # produce strain matches with threshold
  strain_matches <- strain_dna_dist %>% 
    filter(Sample1_details != Sample2_details, # remove same sample comparisons
           Dist_norm <= threshold) # apply genetic threshold for strain match
  
  # generate all strain matches, using only individual donors (not donor batches)
  strain_matches_individual_donors <- strain_matches %>%
    filter(!grepl("Donor_batch", Sample1_details),    # Sample 1 = post-FMT recipient sample
           !grepl("Donor_batch", Sample2_details))  # Sample 2 = donor sample
  
  # find matching strains present in recipient at baseline and week 8
  placebo_BL_wk8_donor_sample_matches <- strain_matches_individual_donors %>% 
    filter(grepl("Placebo_wk8|Placebo_BL", Sample1_details),    # Sample 1 = recipient samples
           grepl("Donor", Sample2_details)) %>%  # Sample 2 = donor sample
    mutate(Donor_ID = substr(Sample2_details,1,4), # extract donor ID
           Recipient_ID = substr(Sample1_details,1,4), # extract recipient ID
           Timepoint = str_split(Sample1_details, "_", simplify = TRUE)[,3]) %>% # extract timepoint
    select(c(Donor_ID, Sample2, Recipient_ID, Species, Timepoint)) %>% # select relevant columns
    rename(Donor_sample = Sample2) %>% # rename donor sample column
    spread(key = Timepoint, value = Timepoint) %>% # spread the timepoint column
    replace_na(list(BL = "absent", wk8 = "absent")) %>% # replace NAs as strain absent
    mutate(BL = str_replace(BL, "BL", "donor_match"),
           wk8 = str_replace(wk8, "wk8", "donor_match")) %>%
    arrange(Donor_ID, Donor_sample, Recipient_ID, Species) %>%
    group_by(Donor_ID, Recipient_ID, Species) %>%
    filter(any(BL == "donor_match") & any(wk8 == "donor_match")) %>%
    ungroup() %>%
    arrange(Donor_ID, Recipient_ID, Donor_sample, Species)
  
  # baseline-matching strains to subtract from data
  placebo_BL_subtraction <- placebo_BL_wk8_donor_sample_matches %>%
    select(Donor_ID, Recipient_ID, Species) %>%
    distinct()
  
  # subtract matching strains also at baseline
  placebo_engraftment_individual_donors_subBL <- strain_matches_individual_donors %>% 
    filter(grepl("Placebo_wk8", Sample1_details),    # Sample 1 = post-placebo recipient sample
           grepl("Donor", Sample2_details)) %>%  # Sample 2 = donor sample
    mutate(Donor_ID = substr(Sample2_details,1,4), # extract donor ID
           Recipient_ID = substr(Sample1_details,1,4)) %>%  # extract recipient ID
    anti_join(placebo_BL_subtraction, by = c("Donor_ID", "Recipient_ID", "Species")) # placebo week 8 donor matches not also in recipient at BL
  
  # placebo matching strains to omit from results
  placebo_engrafted_spp <- placebo_engraftment_individual_donors_subBL %>%
    select(c(Donor_ID, Species)) %>%
    distinct()
  
  # remove placebo matching strain species from FMT engraftment data
  fmt_engrafted_sub_placebo_spp <- 
    sea_output %>%
    anti_join(placebo_engrafted_spp, by = c("Species", "Donor_ID"))
  
  fmt_engrafted_sub_placebo_spp # return result
}

# function to convert predicted data into a presence absence matrix that includes all donors and recipients
presence_absence_matrix <- function(sea_output, donor_ids, fmtwk8_ids){
  # takes output from strain engraftment analysis and formats as presence absence matrix
  # formatted_sea_output (tibble): presence (1) absence (0) matrix including all individual donors and FMT week 8 recipients
  formatted_sea_output <- sea_output %>%
    select(Donor_ID, Recipient_ID) %>% # select predicted donor-recipient pairings from strain engraftment analysis
    distinct() %>% # get distinct pairings - double ups existed due to multiple donor strains detected in a recipient
    # i.e. was it a pairing or was it not
    mutate(value = 1) %>%  # create 'presence' value column for each pairing
    
    right_join(donor_ids, by = c("Donor_ID" = "Participant_ID")) %>% # join with all donor IDs to make sure none missing
    mutate(value = ifelse(is.na(value), 0, value), # mark any missing donor IDs as absent
           Recipient_ID = ifelse(is.na(Recipient_ID), "1001", Recipient_ID)) %>% # replace missing paired recipient ID with arbitrary ID
    
    right_join(fmtwk8_ids, by = c("Recipient_ID" = "Participant_ID")) %>% # join with all FMT wk8 IDs to make sure none missing
    mutate(value = ifelse(is.na(value), 0, value), # mark any missing donor IDs as absent
           Donor_ID = ifelse(is.na(Donor_ID), "D006", Donor_ID)) %>% # replace missing paired donor ID with arbitrary ID
    
    arrange(Recipient_ID) %>% # arrange recipient IDs (will be the new column headers)
    pivot_wider(names_from = Recipient_ID, values_from = value, values_fill = 0) %>% # convert data to wide format
    arrange(Donor_ID) %>% # now arrange donor IDs (rows)
    select(-Donor_ID) # remove donor ID column before making into vector for confusion matrix
  
  formatted_sea_output
}

# function to loop through above strain engraftment analysis functions using defined thresholds
sea_f1_multiple_thresholds <- function(strain_dna_dist, threshold_min, threshold_max, threshold_interval, donor_ids, fmtwk8ids, true_vector){
  # strain_dna_dist: DNA distances from pairwise strian comparisons (from strainphlan output)
  # threshold_min: minimum DNA dist threshold to test (must be > 0)
  # threshold_max: maximum DNA dist threshold to test
  # threshold_interval: DNA dist threshold intervals to test
  # donor_ids: all individual donor IDs
  # fmtwk8ids: all FMT recipient week 8 IDs
  # true_vector: vector form of true pairings presence/absence matrix
  
  thresholds <- seq(threshold_min, threshold_max, by = threshold_interval) # testing these thresholds
  comparison_df <- tibble(Threshold = numeric(), # initialise empty results dataframe to compare thresholds
                          F1 = numeric(), Precision = numeric(), Recall = numeric(), `Balanced Accuracy` = numeric()) 
  for (threshold in thresholds) {
    
    # run strain engraftment analysis
    sea_output <- strain_engraftment_analysis(threshold, strain_dna_dist)
    
    # format strain engraftment analysis results as presence absence matrix including all donor and recipient IDs
    formatted_sea_output <- presence_absence_matrix(sea_output, donor_ids, fmtwk8_ids)
    
    # convert predicted results to a vector for confusion matrix
    predicted_vector <- as.vector(as.matrix(formatted_sea_output))
    
    # create confusion matrix comparing predicted results with the true results
    confusion_matrix <- confusionMatrix(factor(predicted_vector), factor(true_vector), positive="1")
    
    # extract F1 from confusion matrix output
    threshold_F1 <- confusion_matrix$byClass[[7]]
    threshold_precision <- confusion_matrix$byClass[[5]]
    threshold_recall <- confusion_matrix$byClass[[6]]
    threshold_balanced_accuracy <- confusion_matrix$byClass[[11]]
    threshold_sensitivity <- confusion_matrix$byClass[[1]]
    threshold_specificity <- confusion_matrix$byClass[[2]]
    
    # create final results dataframe
    results <- tibble(
      Threshold = threshold,
      F1 = threshold_F1, Precision = threshold_precision, Recall = threshold_recall, `Balanced Accuracy` = threshold_balanced_accuracy,
      Sensitivity = threshold_sensitivity, Specificity = threshold_specificity)
    
    comparison_df <- comparison_df %>%
      bind_rows(results)
  }
  comparison_df %>% mutate(Group = "baseline_subtraction") %>% arrange(desc(F1)) # arrange results with best F1 at top
}

# function to loop through above strain engraftment analysis functions using defined thresholds, with placebo adjustment
sea_f1_multiple_thresholds_placebo_subtraction <- function(strain_dna_dist, threshold_min, threshold_max, threshold_interval, donor_ids, fmtwk8ids, true_vector){
  # strain_dna_dist: DNA distances from pairwise strian comparisons (from strainphlan output)
  # threshold_min: minimum DNA dist threshold to test (must be > 0)
  # threshold_max: maximum DNA dist threshold to test
  # threshold_interval: DNA dist threshold intervals to test
  # donor_ids: all individual donor IDs
  # fmtwk8ids: all FMT recipient week 8 IDs
  # true_vector: vector form of true pairings presence/absence matrix
  
  thresholds <- seq(threshold_min, threshold_max, by = threshold_interval) # testing these thresholds
  comparison_df <- tibble(Threshold = numeric(), # initialise empty results dataframe to compare thresholds
                          F1 = numeric(), Precision = numeric(), Recall = numeric(), `Balanced Accuracy` = numeric()) 
  for (threshold in thresholds) {
    
    # run strain engraftment analysis
    sea_output <- strain_engraftment_analysis(threshold, strain_dna_dist)
    
    # adjust for donor-matching strains in placebos
    sea_output_placebo_adjusted <- strain_engraftment_analysis_placebo_adjusted(threshold, strain_dna_dist, sea_output)
    
    # format strain engraftment analysis results as presence absence matrix including all donor and recipient IDs
    formatted_sea_output_placebo_adjusted <- presence_absence_matrix(sea_output_placebo_adjusted, donor_ids, fmtwk8_ids)
    
    # convert predicted results to a vector for confusion matrix
    predicted_vector_placebo_adjusted <- as.vector(as.matrix(formatted_sea_output_placebo_adjusted))
    
    # create confusion matrix comparing predicted results with the true results
    confusion_matrix_placebo_adjusted <- confusionMatrix(factor(predicted_vector_placebo_adjusted), factor(true_vector), positive="1")
    
    # extract metrics from confusion matrix output
    threshold_F1 <- confusion_matrix_placebo_adjusted$byClass[[7]]
    threshold_precision <- confusion_matrix_placebo_adjusted$byClass[[5]]
    threshold_recall <- confusion_matrix_placebo_adjusted$byClass[[6]]
    threshold_balanced_accuracy <- confusion_matrix_placebo_adjusted$byClass[[11]]
    threshold_sensitivity <- confusion_matrix_placebo_adjusted$byClass[[1]]
    threshold_specificity <- confusion_matrix_placebo_adjusted$byClass[[2]]
    
    # create final results dataframe
    results <- tibble(
      Threshold = threshold,
      F1 = threshold_F1, Precision = threshold_precision, Recall = threshold_recall, `Balanced Accuracy` = threshold_balanced_accuracy,
      Sensitivity = threshold_sensitivity, Specificity = threshold_specificity)
    
    comparison_df <- comparison_df %>%
      bind_rows(results)
  }
  comparison_df %>% mutate(Group = "baseline_placebo_subtraction") %>% arrange(desc(F1)) # arrange results with best F1 at top
}

## running

# join true pairings with metadata
true_pairings_metadata <- uc_metadata %>%
  left_join(uc_pairings, by = c("Participant_ID" = "...1")) %>%
  filter(Group == "FMT", Timepoint == "wk8") %>% # all FMT recipients
  mutate(Sample1_details = paste0(Participant_ID, "_", Group, "_wk8")) %>%
  select(-c(Sample_ID, Group, Timepoint, Participant_ID))

# format true pairings data
true_pairings_formatted <- true_pairings_metadata %>%
  pivot_longer(-Sample1_details, names_to = "Donor_ID", values_to = "n") %>%
  filter(n!= 0) %>%
  group_by(Sample1_details) %>%
  mutate(all_donors = sum(n)) %>%
  mutate(donor_probability = n/all_donors) %>%
  mutate(Donor_ID = ifelse(Donor_ID == "D6", paste0(str_sub(Donor_ID, 1, 1), "00", str_sub(Donor_ID, 2)),
                           paste0(str_sub(Donor_ID, 1, 1), "0", str_sub(Donor_ID, 2)))) %>%
  ungroup() %>%
  mutate(Recipient_ID = str_extract(Sample1_details, "^[^_]+"))
#save(true_pairings_formatted, file = "true_pairings_formatted.RData")

# filter UC metadata for individual donor IDs (n = 14) and FMT week 8 IDs (n = 32)
# used later to make sure predicted data includes all individual donors and FMT recipients at week 8
donor_ids <- uc_metadata %>%
  filter(Group == "Donor_individual") %>% # individual donors
  select(Participant_ID) %>% # select participant ID
  distinct() # get distinct rows (multiple donor samples)

fmtwk8_ids <- uc_metadata %>%
  filter(Group == "FMT" & Timepoint == "wk8") %>% # FMT recipients
  select(Participant_ID) %>% # select participant ID
  distinct() # get distinct rows (multiple donor samples)

# format true data into presence absence matrix with same ordering of donor IDs and recipient IDs
formatted_true_pairings <- true_pairings_formatted %>%
  select(Donor_ID, Recipient_ID, n) %>% # donor-recipient matches are already distinct
  
  # same code as in predicted data formatting function
  right_join(donor_ids, by = c("Donor_ID" = "Participant_ID")) %>% # join with all donor IDs to make sure none missing
  mutate(n = ifelse(is.na(n), 0, n), # mark any missing donor IDs as absent
         Recipient_ID = ifelse(is.na(Recipient_ID), "1001", Recipient_ID)) %>% # replace missing paired recipient ID with arbitrary ID
  
  right_join(fmtwk8_ids, by = c("Recipient_ID" = "Participant_ID")) %>% # join with all FMT wk8 IDs to make sure none missing
  mutate(n = ifelse(is.na(n), 0, n), # mark any missing donor IDs as absent
         Donor_ID = ifelse(is.na(Donor_ID), "D006", Donor_ID)) %>% # replace missing paired donor ID with arbitrary ID
  
  arrange(Recipient_ID) %>% # arrange recipient IDs (will be the new column headers)
  pivot_wider(names_from = Recipient_ID, values_from = n, values_fill = 0) %>% # convert data to wide format
  arrange(Donor_ID) %>% # now arrange donor IDs (rows)
  select(-Donor_ID) # remove donor ID column before making into vector for confusion matrix

# convert presence absence matrix (true pairings) into vector
true_vector <- as.vector(as.matrix(formatted_true_pairings))

# run strain engraftment analysis using normalised DNA distance strain matching thresholds 0.001-3, at 0.001 intervals
# these were the most common DNA distance values from the pairwise strain comparisons (see Figure 1)
# use caret confusion matrix package to compare predicted data with true data
# reference: true results (true donor recipient pairings)
# data: predicted results (strain engraftment results for each threshold)

# results for baseline-adjusted strain engraftment analysis
sea_f1_multiple_thresholds_res <- sea_f1_multiple_thresholds(strain_dna_dist, threshold_min = 0.001, threshold_max = 3, threshold_interval = 0.001,
                                                             donor_ids, fmtwk8ids, true_vector)

# results for baseline- and subsequent placebo-adjusted strain engraftment analysis
sea_f1_multiple_thresholds_res_placebo_adjusted <- sea_f1_multiple_thresholds_placebo_subtraction(strain_dna_dist, 
                                                                                                  threshold_min = 0.001, threshold_max = 3, threshold_interval = 0.001,
                                                                                                  donor_ids, fmtwk8ids, true_vector)

# join data for plotting
sea_f1_multiple_thresholds_both_methods <- sea_f1_multiple_thresholds_res %>%
  bind_rows(sea_f1_multiple_thresholds_res_placebo_adjusted) %>%
  arrange(desc(F1)) # order by F1

# F1 value for each threshold tested
threshold_testing_plot <- sea_f1_multiple_thresholds_both_methods %>%
  filter(Threshold >= 0.002) %>% # plot thresholds >= 0.002 to remove long tail 
  ggplot(aes(x = Threshold, y = F1, colour = Group)) +
  geom_line() +
  scale_colour_manual(values = c("royalblue", "orange"),
                      labels = c("Placebo\n+baseline subtraction", "Baseline subtraction"), name = "Adjustment method") +
  scale_x_continuous(breaks = seq(0, 3, by = 1),   # major grid lines
                     minor_breaks = seq(0, 3, by = 0.2)) + # minor grid lines
  xlab("DNA distance threshold (median normalised)") +
  ylab("F1 (comparison of\npredicted vs true pairings)") +
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.position = "top",
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(color="black", size=0.4))

# select optimal threshold
sea_f1_multiple_thresholds_both_methods[1:10,] # 0.267 (with no placebo adjustment) where F1 = 0.557

# run strain engraftment analysis with optimal threshold to predict FMT donor-recipient pairings
sea_output_best <- strain_engraftment_analysis(0.267, strain_dna_dist)

# format strain engraftment analysis results as presence absence matrix including all donor and recipient IDs
formatted_sea_output_best <- presence_absence_matrix(sea_output_best, donor_ids, fmtwk8_ids)

# convert predicted results to a vector for confusion matrix
predicted_vector_best <- as.vector(as.matrix(formatted_sea_output_best))

# create confusion matrix comparing predicted results with the true results
confusion_matrix_best <- confusionMatrix(factor(predicted_vector_best), factor(true_vector), positive="1") # tell caret that 1 = presence, otherwise 0 is default

# compare with confusion matrix for previously used threshold (0.2), which used placebo adjustment in strain engraftment analysis
sea_output_past <- strain_engraftment_analysis(0.2, strain_dna_dist)
sea_output_placebo_adjusted_past <- strain_engraftment_analysis_placebo_adjusted(0.2, strain_dna_dist, sea_output_past)
formatted_sea_output_placebo_adjusted_past <- presence_absence_matrix(sea_output_placebo_adjusted_past, donor_ids, fmtwk8_ids)
predicted_vector_placebo_adjusted_past <- as.vector(as.matrix(formatted_sea_output_placebo_adjusted_past))
confusion_matrix_placebo_adjusted_past <- confusionMatrix(factor(predicted_vector_placebo_adjusted_past), factor(true_vector), positive="1")

# evaluate efficacy of blinded predictions
sea_f1_multiple_thresholds_res_placebo_adjusted %>% filter(Threshold == 0.200) # F1 score = 0.538; precision = 0.805; recall = 0.405

# save data
#write.xlsx(sea_f1_multiple_thresholds_both_methods, "threshold_testing_df.xlsx") # joined F1/threshold table for supplementary information
#save(sea_output_best, formatted_sea_output_best, 
#     confusion_matrix_best, file = "sea_output_optimised.RData") # optimised strain engraftment analysis results


# plot AUCROC for both methods
# need to test all DNA dist thresholds for this plot
max(strain_dna_dist$Dist_norm) # 317.08
# for DNA dist thresholds from 4-318, test these at intervals of 1

# extra results for baseline-adjusted strain engraftment analysis
sea_threshold_testing_extra_res <- sea_f1_multiple_thresholds(strain_dna_dist, threshold_min = 4, threshold_max = 318, threshold_interval = 1,
                                                             donor_ids, fmtwk8ids, true_vector)

# extra results for baseline- and subsequent placebo-adjusted strain engraftment analysis
sea_threshold_testing_extra_res_placebo_adjusted <- sea_f1_multiple_thresholds_placebo_subtraction(strain_dna_dist, 
                                                                                                  threshold_min = 4, threshold_max = 318, threshold_interval = 1,
                                                                                                  donor_ids, fmtwk8ids, true_vector)

# join data for plotting
all_thresholds_both_methods <- sea_threshold_testing_extra_res %>% 
  bind_rows(sea_threshold_testing_extra_res_placebo_adjusted) %>% # extra results for 4-max
  bind_rows(sea_f1_multiple_thresholds_both_methods) %>% # original results for 0-3
  arrange(desc(F1)) # order by F1

sea_f1_multiple_thresholds_both_methods <- read_xlsx("threshold_testing_df.xlsx") # load data
# plot AUCROC 
aucroc_plot <- sea_f1_multiple_thresholds_both_methods %>%
  ggplot(aes(x = Specificity, y = Sensitivity, colour = Group)) +
  geom_line() +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed") +
  xlim(1, 0) + ylim(0, 1) +
  
  scale_colour_manual(values = c("royalblue", "orange"),
                      labels = c("Placebo\n+baseline subtraction", "Baseline subtraction"), name = "Adjustment method") +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(color="black", size=0.4))

# plot FOCUS vs Gut Bugs models
uc_threshold_testing <- sea_f1_multiple_thresholds_both_methods %>%
  mutate(Trial = "FOCUS")

# load gutbugs data
gb_threshold_testing <- read_xlsx("threshold_testing_df.xlsx") %>%
  mutate(Trial = "Gut Bugs")

# join data
joined_threshold_testing <- uc_threshold_testing %>%
  bind_rows(gb_threshold_testing) %>%
  filter(Group == "baseline_subtraction") # both models ended up using only baseline subtraction for optimal threshold

# joined AUCROC
joined_aucroc_plot <- joined_threshold_testing %>%
  ggplot(aes(x = Specificity, y = Sensitivity, colour = Trial)) +
  geom_line() +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed") +
  xlim(1, 0) + ylim(0, 1) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) +
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.position = "top",
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(color="black", size=0.4))

# joined precision recall plot
joined_precision_recall_plot <- joined_threshold_testing %>%
  ggplot(aes(x = Recall, y = Precision, colour = Trial)) +
  geom_line() +
  xlim(0, 1) + ylim(0, 1) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) +
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.position = "top",
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(color="black", size=0.4))