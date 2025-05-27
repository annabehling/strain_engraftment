## code summary:
# choose DNA distance threshold for strain matching
# perform strain engraftment analysis (baseline- and placebo-adjusted results)
# make FMT donor-recipient pairing predictions
# investigate predicted D043 (high engraftment) engrafted strains

## load libraries

library(tidyverse)
library(ape)
library(phangorn)
library(pheatmap)

## load data

load("uc_metadata_subset.RData") # uc_metadata

## load palettes

large_colour_palette_3 <- c("hotpink", "maroon", "orange", "pink2" ,"#1B9E77", "turquoise", "khaki2", "cornflowerblue", "#3B9AB2", "#EBCC2A",
                            "lightblue", "#BEAED4", "lightgreen", "slategray3")


palette <- c("#44AA99", "#cc7fbf", "#DDCC77", "#88CCEE", "#CC6677", "#6699CC", "#888888")

## running

## strainphlan3

# read in alignment sequences
fasta_files <- list.files(path = "strainphlan3/alignments/", pattern = "\\.aln")

# create empty tibbles to store results
strain_counts <- tibble()
sample_strains <- tibble()
strain_dna_dist <- tibble()

# for each species with strainphlan data, calculate number of strains/samples and pairwise DNA distances
for (f in fasta_files) {
  
  seq <- read.dna(paste("strainphlan3/alignments/", f, sep = ""), format="fasta")
  
  strain_counts <- 
    tibble(Species = f, Total = nrow(seq)) %>%
    mutate(Species = str_remove(Species, ".StrainPhlAn3_concatenated.aln")) %>% 
    bind_rows(strain_counts)
  
  sample_strains <- 
    tibble(Sample_ID = rownames(seq), Species = f) %>% 
    mutate(Species = str_remove(Species, ".StrainPhlAn3_concatenated.aln")) %>% 
    bind_rows(sample_strains)
  
  phyDat_seq <- phyDat(seq, type = "DNA", levels = NULL)
  
  dna_dist <- dist.ml(phyDat_seq, model="JC69")
  
  dna_dist_long <- 
    as.matrix(dna_dist) %>%
    data.frame() %>%
    rownames_to_column("Sample1") %>%
    gather(Sample2, Dist, -Sample1) %>%
    mutate(Dist_pseudo = Dist+1e-6, # add small pseudocount to solve issue re. median normalisation when median is 0
           Dist_norm = Dist_pseudo / median(Dist_pseudo)) %>%
    mutate(Species = str_remove(f, ".StrainPhlAn3_concatenated.aln"))
  
  strain_dna_dist <- rbind(dna_dist_long, strain_dna_dist) # rbind to create a super long table for all strains
  
}

# add metadata for sample 1 and 2
strain_dna_dist <- strain_dna_dist %>% 
  mutate(Sample_ID = Sample1) %>% 
  left_join(uc_metadata) %>% 
  mutate(Sample1_details = paste0(Participant_ID, "_", Group, "_", Timepoint)) %>% 
  select(-Sample_ID, -Participant_ID, -Group, -Timepoint) %>% 
  mutate(Sample_ID = Sample2) %>% 
  left_join(uc_metadata) %>% 
  mutate(Sample2_details = paste0(Participant_ID, "_", Group, "_", Timepoint),
         Sample1_details = str_replace(Sample1_details, "batch_Donor", "batch"),
         Sample1_details = str_replace(Sample1_details, "individual_Donor", "individual"),
         Sample2_details = str_replace(Sample2_details, "batch_Donor", "batch"),
         Sample2_details = str_replace(Sample2_details, "individual_Donor", "individual")) %>% 
  select(-Sample_ID, -Participant_ID, -Group, -Timepoint)

# summary descriptions
sum(strain_counts$Total)           # 1,181 individual strains profiles 
n_distinct(sample_strains$Species) # belonging to 36 species

strain_counts_sample <- sample_strains %>% 
  group_by(Sample_ID) %>% 
  count() 

mean(strain_counts_sample$n) # average of 8 strains identified/sample
range(strain_counts_sample$n) # minimum of 1 - maximum of 17

# save data
#save(strain_dna_dist, strain_counts, sample_strains, strain_counts_sample, file = "strainphlan.Rdata")

## donor strain engraftment analysis

# density plot to decide the threshold to use as a "strain-match"
# remove strain matches where sample 1 or sample 2 is donor batch (misleading in terms of intra- or inter subject matches)
density_plot_stringent <- strain_dna_dist %>% 
  filter(Sample1 != Sample2) %>% # removes matches to the same sample
  filter(!grepl("Donor_batch", Sample1_details), # removes matches to donor batches
         !grepl("Donor_batch", Sample2_details)) %>%
  mutate(Subject_comparison = ifelse(substr(Sample1_details, 1,4) == substr(Sample2_details, 1,4), "intra", "inter")) %>% 
  filter(Dist_norm < 5) %>% # remove large tail 
  ggplot(aes(x = Dist_norm, fill = Subject_comparison)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0.2, linetype="dashed") + # draw line at strain match threshold
  theme_bw() +
  scale_fill_manual(values = c("slategray3","indianred3"), name = "Subject comparison") +
  xlab("DNA distance\n(median normalised)") + ylab("Density") +
  theme_bw() +
  theme(axis.text.x = element_text(),
        legend.direction = "horizontal",
        legend.position = "top",
        #legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4)) #formats background

## identify donor matching strains in recipients

# the majority of intra-subject strain comparisons have a Dist_norm < 0.2 while the majority of inter-subject 
# strain comparisons have a Dist_norm > 0.2 - therefore will use 0.2 as our strain identity threshold
strain_matches <- strain_dna_dist %>% 
  filter(Sample1_details != Sample2_details, # remove same sample comparisons
         Dist_norm <=0.2) # apply genetic threshold for strain match # less stringent would be 0.3

# generate all strain matches, using only individual donors (not donor batches)
strain_matches_individual_donors <- strain_matches %>%
  filter(!grepl("Donor_batch", Sample1_details),    # Sample 1 = post-FMT recipient sample
         !grepl("Donor_batch", Sample2_details))  # Sample 2 = donor sample

# first identify donor matching strains at BL (these won't count as donor matching if present at week 8)
# FMT recipients
FMT_BL_wk8_donor_sample_matches <- strain_matches_individual_donors %>% 
  filter(grepl("FMT_wk8|FMT_BL", Sample1_details),    # Sample 1 = recipient samples
         grepl("Donor", Sample2_details)) %>%  # Sample 2 = donor sample
  mutate(Donor_ID = substr(Sample2_details,1,4),
         Recipient_ID = substr(Sample1_details,1,4),
         Timepoint = str_split(Sample1_details, "_", simplify = TRUE)[,3]) %>% 
  select(c(Donor_ID, Sample2, Recipient_ID, Species, Timepoint)) %>%
  rename(Donor_sample = Sample2) %>%
  spread(key = Timepoint, value = Timepoint) %>% # spread the timepoint column
  replace_na(list(BL = "absent", wk8 = "absent")) %>%
  mutate(BL = str_replace(BL, "BL", "donor_match"),
         wk8 = str_replace(wk8, "wk8", "donor_match")) %>%
  arrange(Donor_ID, Donor_sample, Recipient_ID, Species) %>%
  group_by(Donor_ID, Recipient_ID, Species) %>%
  filter(any(BL == "donor_match") & any(wk8 == "donor_match")) %>%
  ungroup() %>%
  arrange(Donor_ID, Recipient_ID, Donor_sample, Species)

# placebo recipients
placebo_BL_wk8_donor_sample_matches <- strain_matches_individual_donors %>% 
  filter(grepl("Placebo_wk8|Placebo_BL", Sample1_details),    # Sample 1 = recipient samples
         grepl("Donor", Sample2_details)) %>%  # Sample 2 = donor sample
  mutate(Donor_ID = substr(Sample2_details,1,4),
         Recipient_ID = substr(Sample1_details,1,4),
         Timepoint = str_split(Sample1_details, "_", simplify = TRUE)[,3]) %>% 
  select(c(Donor_ID, Sample2, Recipient_ID, Species, Timepoint)) %>%
  rename(Donor_sample = Sample2) %>%
  spread(key = Timepoint, value = Timepoint) %>% # spread the timepoint column
  replace_na(list(BL = "absent", wk8 = "absent")) %>%
  mutate(BL = str_replace(BL, "BL", "donor_match"),
         wk8 = str_replace(wk8, "wk8", "donor_match")) %>%
  arrange(Donor_ID, Donor_sample, Recipient_ID, Species) %>%
  group_by(Donor_ID, Recipient_ID, Species) %>%
  filter(any(BL == "donor_match") & any(wk8 == "donor_match")) %>%
  ungroup() %>%
  arrange(Donor_ID, Recipient_ID, Donor_sample, Species)

# FMT matching strains to omit from results
FMT_BL_subtraction <- FMT_BL_wk8_donor_sample_matches %>%
  select(Donor_ID, Recipient_ID, Species) %>%
  distinct()

# placebo matching strains to omit from results
placebo_BL_subtraction <- placebo_BL_wk8_donor_sample_matches %>%
  select(Donor_ID, Recipient_ID, Species) %>%
  distinct()

# find baseline-adjusted donor strain matches in FMT and placebo recipients
# get strain engraftment results at week 8 post-intervention
# FMT recipients
FMT_engraftment_individual_donors_subBL <- strain_matches_individual_donors %>% 
  filter(grepl("FMT_wk8", Sample1_details),    # Sample 1 = post-FMT recipient sample
         grepl("Donor", Sample2_details)) %>%  # Sample 2 = donor sample
  mutate(Donor_ID = substr(Sample2_details,1,4),
         Recipient_ID = substr(Sample1_details,1,4)) %>% 
  anti_join(FMT_BL_subtraction, by = c("Donor_ID", "Recipient_ID", "Species")) # subtract matching strains also at baseline

# placebo recipients
placebo_engraftment_individual_donors_subBL <- strain_matches_individual_donors %>% 
  filter(grepl("Placebo_wk8", Sample1_details),    # Sample 1 = post-placebo recipient sample
         grepl("Donor", Sample2_details)) %>%  # Sample 2 = donor sample
  mutate(Donor_ID = substr(Sample2_details,1,4),
         Recipient_ID = substr(Sample1_details,1,4)) %>% 
  anti_join(placebo_BL_subtraction, by = c("Donor_ID", "Recipient_ID", "Species")) # subtract matching strains also at baseline

# plot baseline-adjusted donor strain matching results as heatmap
# define custom scale for both plots (FMT and placebo)
colour_breaks <- c(seq(0, 0.998, length.out = 2), seq(0.999, 4, length.out = 50))
colour_values <- c("#dddddd", colorRampPalette(c("#3B9AB2", "slategray3", "#EBCC2A", "indianred3"))(51))

# combine FMT and placebo data for plotting
all_recipients_engraftment_individual_donors_sub_BL <- FMT_engraftment_individual_donors_subBL %>%
  bind_rows(placebo_engraftment_individual_donors_subBL) %>%
  mutate(Group = case_when(grepl("FMT", Sample1_details) ~ "FMT",
                           grepl("Placebo", Sample1_details) ~ "Placebo"))

all_recipients_engraftment_individual_donors_sub_BL_matrix_data <- all_recipients_engraftment_individual_donors_sub_BL %>%
  select(Recipient_ID, Donor_ID, Species, Group) %>% # get distinct donor matching strain species for each pairing in each group
  distinct() %>% # count the number of recipient strains that match each donor
  group_by(Recipient_ID, Donor_ID, Group) %>% 
  count() %>% # count n distinct donor matching strain species for each pairing in each group
  ungroup() %>%
  spread(Donor_ID, n) %>%
  arrange(Group, Recipient_ID) 

all_recipients_engraftment_individual_donors_sub_BL_matrix <- all_recipients_engraftment_individual_donors_sub_BL_matrix_data %>%
  select(-Group) %>%
  column_to_rownames("Recipient_ID") %>%
  replace(is.na(.), 0)

sub_BL_heatmap_annotation <- all_recipients_engraftment_individual_donors_sub_BL_matrix_data %>% 
  select(Group, Recipient_ID) %>%
  column_to_rownames("Recipient_ID")

sub_BL_heatmap_colours <- list(Group = c("FMT" = "#b2182b", "Placebo" = "#2166ac"))

# plot baseline-adjusted heatmap
all_recipients_BLcorrected_heatmap <- pheatmap(all_recipients_engraftment_individual_donors_sub_BL_matrix,
                                               color = colour_values, 
                                               breaks = colour_breaks,
                                               border_color = "white", 
                                               gaps_row = 28, # add a gap between FMT and placebo recipients
                                               annotation_row = sub_BL_heatmap_annotation,
                                               annotation_colors = sub_BL_heatmap_colours,
                                               annotation_names_row = F,
                                               show_rownames = T, 
                                               show_colnames = T,
                                               cluster_cols = F, 
                                               cluster_rows = F, 
                                               scale="none",
                                               cellheight = 9, 
                                               cellwidth = 9,
                                               legend = T)

# placebo-adjusted FMT strain engraftment

# deduct placebo matching strain species from FMT engraftment data
placebo_engrafted_spp <- placebo_engraftment_individual_donors_subBL %>%
  select(c(Donor_ID, Species)) %>%
  distinct()

fmt_engrafted_sub_placebo_spp <- FMT_engraftment_individual_donors_subBL %>% 
  anti_join(placebo_engrafted_spp, by = c("Species", "Donor_ID"))

# plot baseline-adjusted and placebo-adjusted donor strain matching results for FMT recipient as heatmap
FMT_engraftment_individual_donors_sub_BL_placebo_matrix <- fmt_engrafted_sub_placebo_spp %>%
  select(Recipient_ID, Donor_ID, Species) %>% # get distinct donor matching strain species for each pairing
  distinct() %>% # count the number of recipient strains that match each donor
  group_by(Recipient_ID, Donor_ID) %>% 
  count() %>% # count n distinct donor matching strain species for each pairing in each group
  ungroup() %>%
  spread(Donor_ID, n) %>% 
  arrange(Recipient_ID) %>% 
  column_to_rownames("Recipient_ID") %>%
  replace(is.na(.), 0)

FMT_recipients_BLplacebocorrected_heatmap <- pheatmap(FMT_engraftment_individual_donors_sub_BL_placebo_matrix,
                                                      color = colour_values, 
                                                      breaks = colour_breaks,
                                                      border_color = "white", 
                                                      show_rownames = T, 
                                                      show_colnames = T,
                                                      cluster_cols = F, 
                                                      cluster_rows = F, 
                                                      scale="none",
                                                      cellheight = 9, 
                                                      cellwidth = 9,
                                                      legend = T)

# predict FMT donor-recipient pairings from baseline- and placebo-adjusted strain matching

# format data for plot
blind_pairing_predictions <- fmt_engrafted_sub_placebo_spp %>%
  select(Species, Donor_ID, Recipient_ID) %>% # get distinct engrafted strain species for each donor-recipient pairing
  distinct() %>%
  group_by(Recipient_ID, Donor_ID) %>% 
  count() %>% # summarise number of predicted engraftment events between each donor-recipient pairing 
  ungroup() %>%
  group_by(Recipient_ID) %>% # for each recipient,
  mutate(all_strains = sum(n)) %>% # count number of strain matches across all their predicted donors
  mutate(donor_probability = n/all_strains) # donors with more strain matches in a recipient have higher confidence

# plot as proportional stacked bar
blind_pairing_predictions_plot <- blind_pairing_predictions %>%
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

# count number of recipients with donor matching strains
n_recipients <- blind_pairing_predictions %>%
  select(Recipient_ID) %>%
  distinct()

# count donor range for recipients
donor_range <- blind_pairing_predictions %>%
  select(Recipient_ID, Donor_ID) %>%
  group_by(Recipient_ID) %>%
  count()

# investigate predicted D043 engrafted strains in FMT recipients (week 8)
D043_FMT_engraftment <- strain_matches_individual_donors %>% 
  filter(grepl("FMT_wk8", Sample1_details),    # Sample 1 = post-FMT recipient sample
         grepl("D043_Donor_individual", Sample2_details)) # Sample 2 = donor sample

# check how many recipients had Prevotella copri from D043 
D043_Pcopri_recipients <- D043_FMT_engraftment %>%
  mutate(Recipient_ID = substr(Sample1_details,1,4)) %>% 
  select(Species, Recipient_ID) %>%
  filter(Species == "s__Prevotella_copri") %>%
  select(Recipient_ID) %>%
  distinct()
nrow(D043_Pcopri_recipients) # 15

# plot predicted D043 engrafted strains
D043_FMT_engraftment_plot <- D043_FMT_engraftment %>%
  mutate(Recipient_ID = substr(Sample1_details,1,4)) %>% 
  select(Species, Recipient_ID) %>%
  mutate(Species = str_remove(Species, "s__")) %>%
  mutate(Species = str_replace_all(Species, "_", " ")) %>%
  group_by(Recipient_ID, Species) %>%
  count() %>% # number of each species matching in each recipient
  mutate(n = 1) %>%
  ggplot(aes(x = Recipient_ID, y = n, fill = Species)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = palette, name = "Species") +
  xlab("FMT recipient ID") + ylab("Engrafted D043\ndonor strains") +
  guides(fill = guide_legend(ncol = 2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), # rotate x axis labels 45 degrees
        legend.direction = "horizontal",
        legend.position = "top",
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))


# for supplementary, show prediction results for strain match <= 0.3
density_plot_relaxed <- strain_dna_dist %>% 
  filter(Sample1 != Sample2) %>% # removes matches to the same sample
  filter(!grepl("Donor_batch", Sample1_details), # removes matches to donor batches
         !grepl("Donor_batch", Sample2_details)) %>%
  mutate(Subject_comparison = ifelse(substr(Sample1_details, 1,4) == substr(Sample2_details, 1,4), "intra", "inter")) %>% 
  filter(Dist_norm < 5) %>% # remove large tail 
  ggplot(aes(x = Dist_norm, fill = Subject_comparison)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0.3, linetype="dashed") + # draw line at strain match threshold
  theme_bw() +
  scale_fill_manual(values = c("slategray3","indianred3"), name = "Subject comparison") +
  xlab("DNA distance\n(median normalised)") + ylab("Density") +
  theme_bw() +
  theme(axis.text.x = element_text(),
        legend.direction = "horizontal",
        legend.position = "top",
        #legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4)) #formats background

# rerun above code with strain matching threshold of 0.3 to show effect of relaxed threshold on donor recipient pairings
