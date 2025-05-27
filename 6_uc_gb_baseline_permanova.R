## code summary:
# calculate beta diversity for all (FMT and placebo) recipient species RA at baseline in FOCUS and Gut Bugs Trials
# perform PERMANOVA of all (FMT and placebo) recipient species RA at baseline in FOCUS and Gut Bugs Trials
# perform PERMANOVA on same dataset but apply different colours to Gut Bugs recipients with and without MetS

## load libraries

library(tidyverse)
library(vegan)

## load data

# load metadata

load("uc_metadata_subset.RData") # uc_metadata

load("gb_metadata_subset.RData") # gb_metadata

# load metaphlan species data

load("metaphlan.RData") # load FOCUS metaphlan3 data
uc_species <- species # rename df

load("metaphlan.RData") # load Gut Bugs metaphlan3 data
gb_species <- species # rename df

# join species data for the two cohorts
joined_species <- bind_rows(uc_species, gb_species) %>% replace(is.na(.), 0) # replace missing data with 0 relative abundance

# subset metadata for baseline samples for the two cohorts
uc_baseline_samples <- uc_metadata %>%
  filter(Timepoint == "BL")

gb_baseline_samples <- gb_metadata %>%
  filter(Timepoint == "BL")

# join metadata for the two cohorts
joined_baseline_samples <- bind_rows(uc_baseline_samples, gb_baseline_samples)

# filter the species relative abundance for just the baseline samples
joined_species_baseline <- joined_species %>%
  rownames_to_column("Sample_ID") %>%
  filter(Sample_ID %in% joined_baseline_samples$Sample_ID) %>%
  column_to_rownames("Sample_ID")

# calculate beta diversity for subset of species relative abundance data
joined_beta_div <- vegdist(joined_species_baseline, method = "bray") # Bray Curtis dissimilarity index

# convert to long format
joined_beta_div_long <- joined_beta_div %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample1") %>%
  gather(Sample2, BrayCurtis, -Sample1) %>%
  filter(!(Sample1 == Sample2))  # remove rows where samples are compared to themselves
#save(beta_div, beta_div_long, file = "beta_div.Rdata")

# perform PCoA
joined_pca_res <- cmdscale(joined_beta_div) # perform PCoA
joined_pca_df <- as.data.frame(joined_pca_res) # convert PCoA results to df format
colnames(joined_pca_df) <- c("PC1", "PC2")

# rejoin pcoa data with metadata
uc_pcoa_metadata <- joined_pca_df %>%
  rownames_to_column("Sample_ID") %>%
  inner_join(uc_metadata) %>%
  filter(Timepoint == "BL") %>%
  mutate(Trial = "FOCUS")

gb_pcoa_metadata <- joined_pca_df %>%
  rownames_to_column("Sample_ID") %>%
  left_join(gb_metadata) %>%
  filter(Timepoint == "BL") %>%
  mutate(Trial = "Gut Bugs")

# join data
joined_pcoa_metadata <- bind_rows(uc_pcoa_metadata, gb_pcoa_metadata)

# plot beta diversity 
joined_beta_div_plot <- joined_pcoa_metadata %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Trial), colour = "black", pch = 21, size = 3, alpha = 0.7) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) +
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

# data for permanova

# subset of FOCUS and Gut Bugs metadata for recipients (FMT and placebo) at baseline only
joined_BL_metadata <- joined_pcoa_metadata %>%
  select(Sample_ID, Trial)

# subset of RA data
joined_species_baseline

# run permanova, features: Donor_recipient
# all recipients at BL only vs individual donor means
set.seed(1234)
adonis2(joined_species_baseline ~ Trial, 
        data = joined_BL_metadata, 
        by = "margin", # don't need to worry about order of features specified after '~' above
        permutations = 999, 
        method = "bray") # p = 0.001***, R2 = 0.0626

# remake figure but with MetS data for Gut Bugs participant to see if that explains their spread

load("sample_metadata.Rdata") # expanded metadata

extended_gb_metadata <- joined_pcoa_metadata %>%
  left_join(metadata %>% select(Sample_ID, Timepoint, MetS), by = c("Sample_ID", "Timepoint")) %>%
  mutate(Trial_info = paste0(Trial, "_", MetS)) %>%
  mutate(Trial_info_plot = case_when(Trial_info == "FOCUS_NA" ~ "FOCUS",
                                     Trial_info == "Gut Bugs_NA" ~ "Gut Bugs (no MetS info)",
                                     Trial_info == "Gut Bugs_Yes" ~ "Gut Bugs (MetS)",
                                     Trial_info == "Gut Bugs_No" ~ "Gut Bugs (no MetS)"))

joined_beta_div_plot_mets <- extended_gb_metadata %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Trial_info_plot), colour = "black", pch = 21, size = 3, alpha = 0.7) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#999999", "lightblue"), name = "Trial") +
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
