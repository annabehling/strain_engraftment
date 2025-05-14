## code summary:
# perform strain engraftment analysis using optimal threshold
# make FMT donor-recipient pairing predictions
# calculate donor recall
# calculate donor engraftment efficiency
# identify donors with high and low engraftment efficiency

## load libraries

library(tidyverse)
library(ggbeeswarm)
library(dunn.test)

## load data

load("gb_metadata_subset.RData") # uc_metadata

load("sea_output_optimised.RData") # optimised sea results

load("gb_pairings.RData") # true pairings

load("strainphlan.Rdata") # strainphlan output

## load palettes

large_colour_palette_6 <- c("#009E73", "#F0E442", "#56B4E9", "#CC79A7", "#D55E00", "#000000", "#E69F00", "#0072B2", "#999999")

## running

# add participant sex to metadata
gb_metadata_plus <- gb_metadata %>%
  mutate(Sex = ifelse(substr(Participant_ID,2,2) == "F", "Female", "Male"))

# quantify recipients and donors in predicted pairings
n_pred_donors <- sea_output_best %>%
  select(Donor_ID) %>% # individual donors
  distinct() # 9/9

n_pred_recip <- sea_output_best %>%
  select(Recipient_ID) %>% # FMT recipients
  distinct() # 39/39

# count donor range for recipients
donor_range <- sea_output_best %>%
  select(Recipient_ID, Donor_ID) %>%
  distinct() %>%
  group_by(Recipient_ID) %>%
  count()

# format data for plot
best_donor_batch_predictions_prob <- sea_output_best %>%
  select(Species, Donor_ID, Recipient_ID) %>% # get distinct engrafted strain species for each donor-recipient pairing
  distinct() %>%
  group_by(Recipient_ID, Donor_ID) %>% 
  count() %>% # count n distinct engrafted species for each pairing
  ungroup() %>%
  group_by(Recipient_ID) %>% # for each recipient,
  mutate(all_strains = sum(n)) %>% # count number of distinct strain species matches across all their predicted donors
  mutate(donor_probability = n/all_strains) %>% # donors with more strain matches in a recipient have higher confidence
  ungroup()

# plot optimal predicted donor-recipient pairings as proportional stacked bar (not all are true positives)
unblind_pairing_predictions_plot <- best_donor_batch_predictions_prob %>%
  left_join(gb_metadata_plus %>% select(Participant_ID, Sex) %>% distinct(), by = c("Recipient_ID" = "Participant_ID")) %>%
  ggplot(aes(x = Recipient_ID, y = donor_probability, fill = Donor_ID)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust=0.5), colour="white") +
  facet_grid(~Sex, scales = "free", space = "free") +
  scale_fill_manual(values = large_colour_palette_6, name = "Donor ID") +
  xlab("Recipient ID") + ylab("Distribution of predicted donor matching strains") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), # rotate x axis labels 90 degrees
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(color = "black"))

# join true pairings with metadata
true_pairings_metadata <- gb_metadata %>%
  left_join(gb_pairings, by = c("Participant_ID" = "...1")) %>%
  filter(Group == "FMT", Timepoint == "wk6") %>% # all FMT recipients
  mutate(Sample1_details = paste0(Participant_ID, "_", Group, "_wk6")) %>%
  select(-c(Sample_ID, Group, Timepoint, Participant_ID))

# format true pairings data
true_pairings_formatted <- true_pairings_metadata %>%
  pivot_longer(-Sample1_details, names_to = "Donor_ID", values_to = "n") %>%
  filter(n!= 0) %>%
  group_by(Sample1_details) %>%
  mutate(all_donors = sum(n)) %>%
  mutate(donor_probability = n/all_donors) %>%
  ungroup() %>%
  mutate(Recipient_ID = str_extract(Sample1_details, "^[^_]+"))
#save(true_pairings_formatted, file = "true_pairings_formatted.RData") # true pairings data

# count number of times each donor was used
true_pairings_count <- true_pairings_formatted %>%
  select(Donor_ID, Recipient_ID) %>%
  group_by(Donor_ID) %>% # count true recipients for each donor
  count()
mean(true_pairings_count$n) # 17
max(true_pairings_count$n) # 24
min(true_pairings_count$n) # 5

# calculate donor recall (true positive rate: TP/(TP+FN))
donor_recall <- true_pairings_formatted %>% # true pairings
  mutate(Recipient_ID = substr(Sample1_details,1,4)) %>%
  select(Recipient_ID, Donor_ID, n) %>%
  left_join(best_donor_batch_predictions_prob %>% select(Recipient_ID, Donor_ID, n), by = c("Recipient_ID", "Donor_ID")) %>% # predicted pairings
  rename(true_n = n.x, # rename count columns
         predicted_n = n.y) %>%
  mutate(predicted_n = ifelse(!is.na(predicted_n), 1, 0)) %>% # replace engrafted strain count with presence / absence, na = no strain match
  select(-Recipient_ID) %>% # deselect recipient data to focus on donors
  group_by(Donor_ID) %>% # group by donor ID
  mutate(donor_donations = sum(true_n),
         donor_detections = sum(predicted_n), # only true positives captured since prediction data was left joined with the true data by donor and recipient ID
         donor_true_detection_rate = donor_detections/donor_donations*100) %>%
  ungroup() %>% # ungroup donor data
  select(-c(true_n, predicted_n)) %>% # deselect temporary columns
  distinct() # get distinct rows

# plot violin plot to show average true prediction rate across all donors
donor_recall_plot <- donor_recall %>%
  ggplot(aes(x = "", y = donor_true_detection_rate)) +
  geom_violin() +
  geom_quasirandom(alpha = 0.7, shape = 21, size = 5, aes(fill = Donor_ID)) +
  scale_fill_manual(values = large_colour_palette_6, name = "Donor ID") +
  ylim(0, 100) +
  xlab(NULL) + ylab("Donor recall (%)") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(color = "black"))

# donor recall stats
mean(donor_recall$donor_true_detection_rate) # mean = 79.7
sd(donor_recall$donor_true_detection_rate) # sd = 22.1

# calculate donor engraftment efficiency
# first count the number of unique strain species across all donated samples for each donor
all_donor_strains <- gb_metadata %>%
  inner_join(sample_strains, by = "Sample_ID") %>%
  filter(Group == "Donor") %>%
  select(Participant_ID, Species) %>%
  distinct() %>% 
  group_by(Participant_ID) %>% 
  count() %>% 
  rename(Donor_ID = Participant_ID, Total = n)

# then calculate the number of donor engrafted strains in each recipient
# do only for subset of predicted donor-recipient pairings that were true positive predictions (true pairings)
true_pairings <- true_pairings_formatted %>%
  select(Donor_ID, Recipient_ID) %>%
  mutate(Pairing = paste0(Donor_ID, "_", Recipient_ID)) %>%
  select(Pairing)

engraftment_counts <- sea_output_best %>% # using optimal predicted pairings data
  select(Species, Donor_ID, Recipient_ID) %>% # get distinct engrafted strain species for each donor-recipient pairing
  distinct() %>%
  group_by(Recipient_ID, Donor_ID) %>% 
  count() %>% # count n distinct engrafted species for each pairing
  ungroup() %>%
  spread(Donor_ID, n, fill = 0) %>% # need to respread then gather to add 0 counts
  gather(Donor_ID, Engrafted, -Recipient_ID) %>%
  mutate(Pairing = paste0(Donor_ID, "_", Recipient_ID)) %>%
  right_join(true_pairings) %>% # join with true pairings data to also capture recipients without any donor engraftment detected
  mutate(Engrafted = ifelse(is.na(Engrafted), 0, Engrafted), # fill in missing data
         Recipient_ID = ifelse(is.na(Recipient_ID), str_split(Pairing, "_", simplify = TRUE)[, 2], Recipient_ID),
         Donor_ID = ifelse(is.na(Donor_ID), str_split(Pairing, "_", simplify = TRUE)[, 1], Donor_ID)) 

engraftment_efficiency <- left_join(engraftment_counts, all_donor_strains) %>% 
  mutate(Efficiency = (Engrafted/Total)*100) %>%
  group_by(Donor_ID) %>%
  mutate(mean_engraftment_efficiency = mean(Efficiency, na.rm = TRUE)) %>%
  ungroup()

engraftment_efficiency_plot <- engraftment_efficiency %>%
  left_join(gb_metadata_plus %>% select(Participant_ID, Sex) %>% distinct(), by = c("Donor_ID" = "Participant_ID")) %>%
  ggplot(aes(x = Donor_ID, y = Efficiency, fill = Donor_ID)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.7, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  facet_grid(~Sex, scales = "free", space = "free") +
  scale_fill_manual(values = large_colour_palette_6, name = "Donor ID") +
  scale_colour_manual(values = large_colour_palette_6, name = "Donor ID") +
  ylim(0, 29) + # leave space for significance bars
  xlab("Donor ID") + ylab("Strain engraftment efficiency (%)") + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# descriptive stats
engraftment_efficiency_stats <- engraftment_efficiency %>% 
  group_by(Donor_ID) %>% 
  summarise(Mean = mean(Efficiency),
            Median = median(Efficiency),
            Min = min(Efficiency),
            Max = max(Efficiency))
max(engraftment_efficiency_stats$Mean) # 14.4
min(engraftment_efficiency_stats$Mean) # 1.4

# count distinct donor engrafted strains
engrafted_strains <- sea_output_best %>% # using optimal predicted pairings data
  select(Species, Donor_ID, Recipient_ID) %>%
  distinct() %>%
  mutate(Pairing = paste0(Donor_ID, "_", Recipient_ID)) %>%
  inner_join(true_pairings) %>% # inner joining because recipients without any engraftment aren't relevant for this count
  select(Species) %>%
  distinct() # 57

# test if engraftment efficiency is different between donors
engraftment_efficiency %>% pull(Efficiency) %>% shapiro.test()
ggplot(engraftment_efficiency, aes(x = Efficiency)) + geom_histogram() # data not normally distributed
kruskal.test(Efficiency ~ Donor_ID, data = engraftment_efficiency) # use non-parametric test
# p < 0.001 ***

# run post-hoc testing to see which donors differ
dunn_test <- dunn.test(engraftment_efficiency$Efficiency, engraftment_efficiency$Donor_ID, method = "bonferroni")
donor_comparisons <- dunn_test$comparisons
donor_pvals <- dunn_test$P.adj

sig_donor_comparisons <- donor_comparisons[donor_pvals < 0.05]
sig_pvals <- donor_pvals[donor_pvals < 0.05]
data.frame(Comparison = sig_donor_comparisons, P.adj = sig_pvals)

# save data
#save(engraftment_efficiency, file = "engraftment_efficiency.RData") # engraftment efficiency data

# test correlation between donor recall and engraftment efficiency
donor_recall_formatted <- donor_recall %>%
  select(Donor_ID, donor_true_detection_rate) %>%
  left_join(gb_metadata, by = c("Donor_ID" = "Participant_ID")) %>%
  filter(Group == "Donor") %>%
  select(Donor_ID, Sample_ID, donor_true_detection_rate) %>%
  rename(donor_recall = donor_true_detection_rate)

# join donor recall and engraftment efficiency
correlation_test_data <- donor_recall_formatted %>%
  select(Donor_ID, donor_recall) %>%
  distinct() %>%
  inner_join(engraftment_efficiency) %>%
  select(Donor_ID, donor_recall, mean_engraftment_efficiency) %>%
  distinct()

# plot donor recall and mean engraftment efficiency correlation
correlation_test_plot <- correlation_test_data %>%
  ggplot(aes(x = mean_engraftment_efficiency, y = donor_recall), colour = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "gray", alpha = 0.2) +
  geom_hline(linetype = "dashed", yintercept = 100) +
  geom_point(aes(fill = Donor_ID), alpha = 0.7, size = 4, shape = 21) +
  scale_fill_manual(values = large_colour_palette_6, name = "Donor ID") +
  xlim(0, NA) +
  xlab("Mean strain engraftment efficiency (%)") + ylab("Donor recall (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        #panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# test pearson's correlation
ee_recall_cor <- cor.test(correlation_test_data$mean_engraftment_efficiency, correlation_test_data$donor_recall) # r = 0.724, p = 0.028 *
ee_recall_cor$estimate^2 # R^2 = 0.524

# plot true FMT donor-recipient pairings for reference
true_pairings_plot <- true_pairings_formatted %>%
  left_join(gb_metadata_plus %>% select(Participant_ID, Sex) %>% distinct(), by = c("Recipient_ID" = "Participant_ID")) %>%
  ggplot(aes(x = Recipient_ID, y = donor_probability, fill = Donor_ID)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~Sex, scales = "free", space = "free") +
  scale_fill_manual(values = large_colour_palette_6, name = "Donor ID") +
  xlab("Recipient ID") + ylab("Donor batch distributions") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), # rotate x axis labels 90 degrees
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(color = "black"))
