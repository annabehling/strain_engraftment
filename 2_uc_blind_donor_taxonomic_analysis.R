## code summary:
# process metaphlan3 data
# calculate donor alpha diversity
# calculate beta diversity and run permanovas
# calculate donor P/B ratio

## load libraries

library(tidyverse)
library(vegan)
library(ggbeeswarm)
library(ggrepel)
library(dunn.test)

## load data

load("uc_metadata_subset.RData") # uc_metadata

## load palettes

large_colour_palette_3 <- c("hotpink", "maroon", "orange", "pink2" ,"#1B9E77", "turquoise", "khaki2", "cornflowerblue", "#3B9AB2", "#EBCC2A",
                            "lightblue", "#BEAED4", "lightgreen", "slategray3")

## running

# process metaphlan3 data
metaphlan <- read_tsv("metaphlan3/uc_metaphlan3_profiles.tsv", skip = 1) 

metaphlan_long <- metaphlan %>% 
  select(-NCBI_tax_id) %>% 
  rename(Taxa = clade_name) %>% 
  gather(Sample_ID, RA, -Taxa) %>% # covert to long format so it's easier to edit the Sample IDs
  mutate(Sample_ID = str_remove(Sample_ID, ".metaphlan"), # remove suffix from Sample ID
         Taxa = str_remove(Taxa, ".*\\|")) # simplify taxa name (removes everything before the last | delimiter)

# isolate taxonomic ranks (and the proportion of reads of unknown taxonomy)
species <- filter(metaphlan_long, grepl("s__", Taxa)) %>% spread(Sample_ID, RA)
genus <- filter(metaphlan_long, grepl("g__", Taxa)) %>% spread(Sample_ID, RA)
family <- filter(metaphlan_long, grepl("f__", Taxa)) %>% spread(Sample_ID, RA)
order <- filter(metaphlan_long, grepl("o__", Taxa)) %>% spread(Sample_ID, RA)
phylum <- filter(metaphlan_long, grepl("p__", Taxa)) %>% spread(Sample_ID, RA)

unknown <- filter(metaphlan_long, grepl("UNKNOWN", Taxa))

# note that the sum of columns for any given sample don't sum to 100 (as we have not included the unknown counts)
colSums(species[,-1])
colSums(genus[,-1])
colSums(family[,-1])
colSums(order[,-1])
colSums(phylum[,-1])

# therefore we need to renormalise the relative abundances to sum to 1 (better for downstream analyses)
species[,-1] <- lapply(species[,-1], function(x){ x/sum(x, na.rm=TRUE)})
genus[,-1] <- lapply(genus[,-1], function(x){ x/sum(x, na.rm=TRUE)})
family[,-1] <- lapply(family[,-1], function(x){ x/sum(x, na.rm=TRUE)})
order[,-1] <- lapply(order[,-1], function(x){ x/sum(x, na.rm=TRUE)})
phylum[,-1] <- lapply(phylum[,-1], function(x){ x/sum(x, na.rm=TRUE)})

# check samples now sum to 1
colSums(species[,-1])
colSums(genus[,-1])
colSums(family[,-1])
colSums(order[,-1])
colSums(phylum[,-1])

# transpose data (format required for diversity metrics)
species <- column_to_rownames(species, "Taxa") %>% t() %>% as.data.frame()
genus <- column_to_rownames(genus, "Taxa") %>% t() %>% as.data.frame()
family <- column_to_rownames(family, "Taxa") %>% t() %>% as.data.frame()
order <- column_to_rownames(order, "Taxa") %>% t() %>% as.data.frame()
phylum <- column_to_rownames(phylum, "Taxa") %>% t() %>% as.data.frame()

# save taxonomy tables
#save(metaphlan, metaphlan_long, species, genus, family, order, phylum, unknown, file = "metaphlan.Rdata")

# calculate alpha diversity
alpha_div <- tibble(Sample_ID = rownames(species),
                    Phyla = specnumber(phylum),                      # no. phyla detected
                    Family = specnumber(family),                     # no. family detected
                    Genera = specnumber(genus),                      # no. genera detected
                    Species = specnumber(species),                   # no. species detected
                    Shannon = vegan::diversity(species, index = "shannon")) # Shannon diversity index
#save(alpha_div, file = "alpha_div.Rdata")

# plot donor alpha diversity
donors_alpha_div_plot <- alpha_div %>%
  left_join(uc_metadata, by = "Sample_ID") %>%
  filter(Group == "Donor_individual") %>%
  ggplot(aes(x = Participant_ID, y = Shannon, fill = Participant_ID)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(fatten = 0.5, alpha = 0.6, outlier.colour = NA, position = position_dodge2(width = 0.75, preserve = "single")) +
  ylim(NA, 3.7) + # make room for significance bars
  scale_fill_manual(values = large_colour_palette_3) +
  scale_colour_manual(values = large_colour_palette_3) +
  xlab("Donor ID") + ylab("Shannon's diversity index") + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), # rotate x axis labels 45 degrees
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# test if alpha diversity is different between donors
donor_alpha_div <- alpha_div %>% left_join(uc_metadata, by = "Sample_ID") %>% filter(Group == "Donor_individual")
donor_alpha_div %>% pull(Shannon) %>% shapiro.test() # p = 0.157
ggplot(donor_alpha_div, aes(x = Shannon)) + geom_histogram() # data normally distributed
anova_ad <- aov(Shannon ~ Participant_ID, data = donor_alpha_div) # use parametric test
summary(anova_ad) 
# p = 0.00226

# run post-hoc testing to see which donors differ
tukey_test_ad <- TukeyHSD(anova_ad)
donor_comparisons_ad <- rownames(tukey_test_ad$Participant_ID)
donor_pvals_ad <- tukey_test_ad$Participant_ID[, 4]

sig_donor_comparisons_ad <- donor_comparisons_ad[donor_pvals_ad < 0.05]
sig_pvals_ad <- donor_pvals_ad[donor_pvals_ad < 0.05]
data.frame(Comparison = sig_donor_comparisons_ad, P.adj = sig_pvals_ad)
# D033-D006 0.01610623
# D041-D006 0.02820954
# D033-D023 0.02614135
# D041-D023 0.04129165
# D053-D033 0.02240203
# D059-D033 0.02550641
# D053-D041 0.03922878
# D059-D041 0.04463722


# calculate beta diversity
beta_div <- vegdist(species, method = "bray") # Bray Curtis dissimilarity index

beta_div_long <- beta_div %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample1") %>%
  gather(Sample2, BrayCurtis, -Sample1) %>%
  filter(!(Sample1 == Sample2))  # remove rows where samples are compared to themselves
#save(beta_div, beta_div_long, file = "beta_div.Rdata")

# perform PCoA
pca_res <- cmdscale(beta_div) # perform PCoA
pca_df <- as.data.frame(pca_res) # convert PCoA results to df format
colnames(pca_df) <- c("PC1", "PC2")

# get mean values for donors with multiple samples
pca_df_donor_mean <- pca_df %>% # take raw principal coordinate data
  rownames_to_column("Sample_ID") %>% # convert row names to sample ID column
  left_join(uc_metadata, by = "Sample_ID") %>%
  group_by(Participant_ID, Group) %>%
  mutate(Sample_ID_new = if_else(str_detect(Group, "D"), paste0(first(Sample_ID), "_mean"), Sample_ID), # replaced with representative sample ID for grouping
         PC1_new = if_else(str_detect(Group, "D"), mean(PC1), PC1),
         PC2_new = if_else(str_detect(Group, "D"), mean(PC2), PC2)) %>% # get mean PC1 and PC2 coordinates for donors only
  ungroup() %>% # ungroup data
  select(Sample_ID_new, PC1_new, PC2_new) %>% # select new columns only
  distinct() %>% # get distinct rows so one set of mean values per donor
  mutate(Sample_ID_new = str_remove(Sample_ID_new, "_mean")) %>% # return donor IDs to normal
  rename(PC1 = PC1_new, PC2 = PC2_new) %>% # return PC column names to normal
  column_to_rownames("Sample_ID_new")

# data for donor vs recipient baseline plot
pca_donormean_BL <- pca_df_donor_mean %>%
  rownames_to_column("Sample_ID") %>%
  left_join(uc_metadata) %>%
  filter(Timepoint == "BL" | Timepoint == "Donor") %>%
  mutate(Facet = "Baseline")

# data for donor vs recipient week 8 plot
pca_donormean_wk8 <- pca_df_donor_mean %>%
  rownames_to_column("Sample_ID") %>%
  left_join(uc_metadata) %>%
  filter(Timepoint == "wk8" | Timepoint == "Donor") %>%
  mutate(Facet = "Week 8")

# join data
pca_plot_data_donormean <- bind_rows(pca_donormean_BL, pca_donormean_wk8)

# plot beta diversity 
beta_div_plot <- pca_plot_data_donormean %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Group), colour = "black", pch = 21, size = 2, alpha = 0.7) +
  geom_label_repel(data = pca_plot_data_donormean %>% filter(Participant_ID == "D043"), # label D043 outlier
                   aes(label = Participant_ID), min.segment.length = 0,
                   nudge_x = -0.2,
                   nudge_y = 0.2,
                   size = 2.5) +
  facet_wrap(~Facet) +
  scale_fill_manual(values = c("darkslateblue", "#fdc407", "#b2182b", "#2166ac"), labels = c("Donor batch", "Donor individual", "FMT", "Placebo")) +
  xlab("PC1") + ylab("PC2") +
  guides(fill = guide_legend(ncol = 2, override.aes = list(label = ""))) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(),
        legend.direction = "horizontal",
        legend.position = "top",
        legend.key.size = unit(0.5, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
#ggsave(file = "betadiv.svg", 
#       plot = beta_div_plot, units = "cm", width=10, height=10) 

# run permanovas
# working with species relative abundance data (metaphlan output) that has been renormalised to account for 'unknown' values not included
# one sample (need donor means) per individual
# don't include donor data except for comparison with all recipients at baseline

# need to get average values for each donor for permanova (shouldn't be multiple samples for one individual)
mean_donor_spp_ra <- species %>% # working with metaphlan species relative abundance data
  rownames_to_column("Sample_ID") %>% # convert row names to sample ID column
  left_join(uc_metadata, by = "Sample_ID") %>% # join with metadata
  filter(Timepoint == "Donor") %>% # get donor data
  group_by(Participant_ID, Group) %>% # group by donor ID and group (batch / individual - some IDs feature in both but are separate samples)
  summarise(across(2:477, mean, na.rm = TRUE)) %>% # summarise mean species RAs for each donor
  ungroup() %>%
  left_join(uc_metadata, by = c("Participant_ID", "Group")) %>%
  group_by(Participant_ID, Group) %>%
  mutate(Sample_ID_new = if_else(str_detect(Group, "D"), first(Sample_ID), Sample_ID)) %>% # replace IDs for mean values with representative sample ID from group
  ungroup() %>%
  select(c(481, 3:478)) %>% # select new sample ID, spp columns
  distinct() %>%
  rename(Sample_ID = Sample_ID_new)

# join with donor mean data with recipient data
all_participant_spp_ra <- species %>% # working with metaphlan species relative abundance data
  rownames_to_column("Sample_ID") %>% # convert row names to sample ID column
  left_join(uc_metadata, by = "Sample_ID") %>% # join with metadata
  filter(Timepoint != "Donor") %>% # get recipient data
  select(1:477) %>% # select sample ID and species columns
  bind_rows(mean_donor_spp_ra)

# return sample ID to row names
all_participant_spp_ra_formatted <- all_participant_spp_ra %>% # working with donor mean and recipient metaphlan species relative abundance data
  column_to_rownames("Sample_ID") 

# format metadata
all_participant_metadata <- uc_metadata %>%
  filter(Sample_ID %in% all_participant_spp_ra$Sample_ID) %>%
  mutate(Donor_recipient = if_else(str_detect(Timepoint, "Donor"), "Donor", "Recipient")) # add extra column for differentiating donors from all recipients

## all recipients at BL only vs individual donor means
# make metadata 
BLrecipients_donors_metadata <- all_participant_metadata %>%
  filter(Group == "Donor_individual" | Timepoint == "BL") %>%
  select(Sample_ID, Donor_recipient)

# subset RA data
BLrecipients_donors_permanova_data <- all_participant_spp_ra_formatted %>%
  rownames_to_column("Sample_ID") %>%
  filter(Sample_ID %in% BLrecipients_donors_metadata$Sample_ID) %>%
  column_to_rownames("Sample_ID")

# run permanova, features: Donor_recipient
# all recipients at BL only vs individual donor means
set.seed(1234)
adonis2(BLrecipients_donors_permanova_data ~ Donor_recipient, 
        data = BLrecipients_donors_metadata, 
        by = "margin", # don't need to worry about order of features specified after '~' above
        permutations = 999, 
        method = "bray") # p = 0.001, R2 = 0.04041

## FMT vs placebo at baseline
# make metadata 
BL_FMT_placebo_metadata <- all_participant_metadata %>%
  filter(Group == "FMT" | Group == "Placebo") %>%
  filter(Timepoint == "BL") %>%
  select(Sample_ID, Group)

# subset RA data
BL_FMT_placebo_permanova_data <- all_participant_spp_ra_formatted %>%
  rownames_to_column("Sample_ID") %>%
  filter(Sample_ID %in% BL_FMT_placebo_metadata$Sample_ID) %>%
  column_to_rownames("Sample_ID")

# run permanova, features: Group
set.seed(1234)
adonis2(BL_FMT_placebo_permanova_data ~ Group, 
        data = BL_FMT_placebo_metadata, 
        by = "margin", # don't need to worry about order of features specified after '~' above
        permutations = 999, 
        method = "bray") # p = 0.847, R2 = 0.01415

## FMT vs placebo at week 8
# make metadata 
wk8_FMT_placebo_metadata <- all_participant_metadata %>%
  filter(Group == "FMT" | Group == "Placebo") %>%
  filter(Timepoint == "wk8") %>%
  select(Sample_ID, Group)

# subset RA data
wk8_FMT_placebo_permanova_data <- all_participant_spp_ra_formatted %>%
  rownames_to_column("Sample_ID") %>%
  filter(Sample_ID %in% wk8_FMT_placebo_metadata$Sample_ID) %>%
  column_to_rownames("Sample_ID")

# run permanova, features: Group
set.seed(1234)
adonis2(wk8_FMT_placebo_permanova_data ~ Group, 
        data = wk8_FMT_placebo_metadata, 
        by = "margin", # don't need to worry about order of features specified after '~' above
        permutations = 999, 
        method = "bray") # p = 0.001, R2 = 0.11988

# calculate Prevotella/Bacteroides ratio for donors
donors_pb_ratio <- genus %>%
  rownames_to_column("Sample_ID") %>% # convert row names to sample ID column
  select(Sample_ID, g__Prevotella, g__Bacteroides) %>% # select relevant columns
  mutate(Prevotella_pseudo = g__Prevotella + 1e-6, # create pseudo counts with +1e-6 to avoid zeros in calculation
         Bacteroides_pseudo = g__Bacteroides + 1e-6) %>%
  mutate(PB_ratio = Prevotella_pseudo/Bacteroides_pseudo) %>% # calculate P/B ratio
  inner_join(uc_metadata) %>% # join with meta data
  filter(Group == "Donor_individual")

# plot donor P/B ratio
donors_pb_ratio_plot <- donors_pb_ratio %>%
  mutate(colour = if_else(Participant_ID == "D043", "D043", "Other")) %>%
  ggplot(aes(x = Participant_ID, y = PB_ratio, fill = Participant_ID)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(fatten = 0.5, alpha = 0.6, outlier.colour = NA, position = position_dodge2(width = 0.75, preserve = "single")) +
  scale_y_continuous(trans='log10', breaks = c(1, 0.01, 1e-4, 1e-6), labels = c(1, 0.01, 1e-04, 1e-06)) +
  scale_fill_manual(values = large_colour_palette_3) +
  scale_colour_manual(values = large_colour_palette_3) +
  geom_hline(yintercept = 1, linetype = 'dotted') +
  xlab("Donor ID") + ylab("Prevotella/Bacteroides ratio") + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), # rotate x axis labels 45 degrees
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# test if P/B ratio is different between donors
donors_pb_ratio %>% pull(PB_ratio) %>% shapiro.test() # p = 4.731e-09 significantly different from normal distribution
ggplot(donors_pb_ratio, aes(x = PB_ratio)) + geom_histogram() # data not normally distributed
kruskal.test(PB_ratio ~ Participant_ID, data = donors_pb_ratio) # use non-parametric test
# p = 0.0285

# run post-hoc testing to see which donors differ
dunn_test_pb <- dunn.test(donors_pb_ratio$PB_ratio, donors_pb_ratio$Participant_ID, method = "bonferroni")
donor_comparisons_pb <- dunn_test_pb$comparisons
donor_pvals_pb <- dunn_test_pb$P.adj

sig_donor_comparisons_pb <- donor_comparisons_pb[donor_pvals_pb < 0.05]
sig_pvals_pb <- donor_pvals_pb[donor_pvals_pb < 0.05]
data.frame(Comparison = sig_donor_comparisons_pb, P.adj = sig_pvals_pb)
# none significantly different after adjusting p values
