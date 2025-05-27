## code summary:
# perform strain engraftment analysis using optimal threshold
# make FMT donor-recipient pairing predictions
# calculate donor recall
# calculate donor engraftment efficiency
# identify donors with high and low engraftment efficiency

## load libraries

library(readxl)
library(tidyverse)
library(ggbeeswarm)
library(dunn.test)

## load data

load("uc_metadata_subset.RData") # uc_metadata

load("sea_output_optimised.RData") # optimised sea results

uc_pairings <- read_excel(path = "FOCUS Patient Donor matching metadata.xlsx", skip = 1, col_names = TRUE) # true pairings

load("strainphlan.Rdata") # strainphlan output

## load palettes

large_colour_palette_3 <- c("hotpink", "maroon", "orange", "pink2" ,"#1B9E77", "turquoise", "khaki2", "cornflowerblue", "#3B9AB2", "#EBCC2A",
                            "lightblue", "#BEAED4", "lightgreen", "slategray3")

## running

# quantify recipients and donors in predicted pairings
n_pred_donors <- sea_output_best %>%
  select(Donor_ID) %>% # individual donors
  distinct() # 14/14

n_pred_recip <- sea_output_best %>%
  select(Recipient_ID) %>% # FMT recipients
  distinct() # 29/32

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
  ggplot(aes(x = Recipient_ID, y = donor_probability, fill = Donor_ID)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust=0.5), colour="white") +
  scale_fill_manual(values = large_colour_palette_3, name = "Donor ID") +
  xlab("Recipient ID") + ylab("Distribution of predicted donor matching strains") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), # rotate x axis labels 90 degrees
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(color = "black"))

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

# count number of times each donor was used
true_pairings_count <- true_pairings_formatted %>%
  select(Donor_ID, Recipient_ID) %>%
  group_by(Donor_ID) %>% # count true recipients for each donor
  count()
mean(true_pairings_count$n) # 12
max(true_pairings_count$n) # 29
min(true_pairings_count$n) # 1

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
  scale_fill_manual(values = large_colour_palette_3, name = "Donor ID") +
  ylim(0, 100) +
  xlab(NULL) + ylab("Donor recall (%)") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(color = "black"))

# donor recall stats
mean(donor_recall$donor_true_detection_rate) # mean = 41.7
sd(donor_recall$donor_true_detection_rate) # sd = 25.3

# investigate correlation between donor recall and usage data
recall_use_data <- donor_recall %>% # join recall and usage data
  full_join(true_pairings_count, by = "Donor_ID") %>%
  rename(donor_recall = donor_true_detection_rate)

# plot donor recall and usage data correlation
correlation_test_plot <- recall_use_data %>%
  ggplot(aes(x = n, y = donor_recall), colour = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "gray", alpha = 0.2) +
  geom_point(aes(fill = Donor_ID), alpha = 0.7, size = 4, shape = 21) +
  scale_fill_manual(values = large_colour_palette_3, name = "Donor ID") +
  ylim(0, 100) +
  xlab("Number of FMT recipients") + ylab("Donor recall (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        #panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# test pearson's correlation
recall_use_cor <- cor.test(recall_use_data$n, recall_use_data$donor_recall) # r = 0.260, p = 0.370
recall_use_cor$estimate^2 # R^2 = 0.0674


# calculate donor engraftment efficiency
# first count the number of unique strain species across all donated samples for each donor
all_donor_strains <- uc_metadata %>%
  inner_join(sample_strains, by = "Sample_ID") %>%
  filter(Group == "Donor_individual") %>%
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

engraftment_efficiency_plot <- ggplot(engraftment_efficiency, aes(x = Donor_ID, y = Efficiency, fill = Donor_ID)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.7, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = large_colour_palette_3, name = "Donor ID") +
  scale_colour_manual(values = large_colour_palette_3, name = "Donor ID") +
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
max(engraftment_efficiency_stats$Mean) # 14.1
min(engraftment_efficiency_stats$Mean) # 0

# count distinct donor engrafted strains
engrafted_strains <- sea_output_best %>% # using optimal predicted pairings data
  select(Species, Donor_ID, Recipient_ID) %>%
  distinct() %>%
  mutate(Pairing = paste0(Donor_ID, "_", Recipient_ID)) %>%
  inner_join(true_pairings) %>% # inner joining because recipients without any engraftment aren't relevant for this count
  select(Species) %>%
  distinct() # 19

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
# consistently low = D033, D040
# consistently high = D043, D053

# save data
#save(engraftment_efficiency, file = "engraftment_efficiency.RData") # engraftment efficiency data

# test correlation between donor recall and engraftment efficiency
donor_recall_formatted <- donor_recall %>%
  select(Donor_ID, donor_true_detection_rate) %>%
  left_join(uc_metadata, by = c("Donor_ID" = "Participant_ID")) %>%
  filter(Group == "Donor_individual") %>%
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
  scale_fill_manual(values = large_colour_palette_3, name = "Donor ID") +
  xlim(NA, 15) +
  xlab("Mean strain engraftment efficiency (%)") + ylab("Donor recall (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        #panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# test pearson's correlation
ee_recall_cor <- cor.test(correlation_test_data$mean_engraftment_efficiency, correlation_test_data$donor_recall) # r = 0.904, p < 0.001 ***
ee_recall_cor$estimate^2 # R^2 = 0.817

# plot true FMT donor-recipient pairings for reference
true_pairings_plot <- true_pairings_formatted %>%
  ggplot(aes(x = Recipient_ID, y = donor_probability, fill = Donor_ID)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = large_colour_palette_3, name = "Donor ID") +
  xlab("Recipient ID") + ylab("Donor batch distributions") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), # rotate x axis labels 90 degrees
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(color = "black"))
