## code summary
# plot donor and recipient alpha diversity, relative to the engraftment efficiency (EE) of each pairing
# correlate mean donor engraftment efficiency with mean donor alpha diversity
# correlate engraftment efficiency with mean donor-recipient species beta diversity
# correlate engraftment efficiency with mean donor-recipient functional beta diversity
# find which factor best predicts engraftment efficiency
# find which functional category dissimilarity best predicts engraftment efficiency

## load libraries

library(tidyverse)
library(ggbeeswarm)
library(vegan)
library(ggnewscale)
library(lme4)
library(lmerTest)
library(viridis)

## load data

load("uc_metadata_subset.RData") # uc_metadata

load("sea_output_optimised.RData") # optimised sea results

load("alpha_div.Rdata") # alpha diversity data

load("engraftment_efficiency.RData") # engraftment efficiency data

load("metaphlan.RData") # metaphlan3 data

load("HQ-gene-mapping-filtered.RData") # HQ MAGs, HQ genes, HQ gene clusters

load("true_pairings_formatted.RData") # true pairings

kneaddata_read_counts <- read.table("kneaddata/kneaddata_read_count_table_nounderscore.tsv", # kneaddata sample read count data
                                    sep = "\t", header = TRUE, fill = TRUE)

## load palettes

large_colour_palette_3 <- c("hotpink", "maroon", "orange", "pink2" ,"#1B9E77", "turquoise", "khaki2", "cornflowerblue", "#3B9AB2", "#EBCC2A",
                            "lightblue", "#BEAED4", "lightgreen", "slategray3")

large_colour_palette_4a <- c("darkgreen", "#3B9AB2", "pink2", "#BF5B17", "lightgreen",  "#1B9E77", "turquoise", 
                             "orange", "lightblue", "indianred3", "#BEAED4", "#386CB0", "khaki2", "chartreuse4", "sienna4", "#EBCC2A", "#7570B3","maroon", 
                             "slategray3", "goldenrod", "darksalmon", "coral", "cornflowerblue", "steelblue4", "sandybrown", "firebrick4", 
                             "royalblue",   "thistle", "yellowgreen", "darkblue", "gray86")

## load functions

# format relative abundance data for PCoA plot
vegan_beta_div_pca <- function(species, index){
  beta_div <- vegdist(species, method = index) # "bray", "jaccard", "euclidean"
  pca_res <- cmdscale(beta_div) # perform PCA
  pca_df <- as.data.frame(pca_res) # convert PCA results to df format
  colnames(pca_df) <- c("PC1", "PC2")
  
  pca_df
}

## running

## alpha diversity

# prepare alpha diversity data for individual donors and their FMT recipients at baseline
donor_recipient_alpha_div <- uc_metadata %>%
  inner_join(alpha_div, by = "Sample_ID") %>% # join with metadata by sample ID to get donor ID
  filter(Group == "Donor_individual") %>% # select donor individual
  group_by(Participant_ID) %>% # group by donor ID
  mutate(mean_shannon = mean(Shannon)) %>% # calculate mean Shannon diversity index across samples for each donor
  ungroup() %>%
  inner_join(engraftment_efficiency, by = c("Participant_ID" = "Donor_ID"), relationship = "many-to-many") %>%
  select(Participant_ID, Sample_ID, Efficiency, mean_engraftment_efficiency, mean_shannon, Recipient_ID) %>%
  rename(Donor_ID = Participant_ID,
         Mean_shannon_donor = mean_shannon) %>%
  # join with metadata by recipient ID to get BL sample IDs for each recipient of donors
  left_join(uc_metadata %>% filter(Timepoint == "BL"), by = c("Recipient_ID" = "Participant_ID")) %>%
  rename(Donor_sample = Sample_ID.x, # rename columns
         Recipient_sample = Sample_ID.y) %>%
  select(-c(Group, Timepoint)) %>% # remove unneeded data
  left_join(alpha_div, by = c("Recipient_sample" = "Sample_ID")) %>% # join with alpha_div by recipient sample ID to get recipient alpha diversity
  select(-c(Phyla, Family, Genera, Species)) %>% # remove unneeded columns
  rename(Shannon_recipient = Shannon)

# plot mean donor alpha diversity vs recipient alpha diversity, relative to engraftment efficiency for each pairing
donor_recipient_alpha_div_plot <- donor_recipient_alpha_div %>%
  ggplot(aes(x = Mean_shannon_donor, y = Shannon_recipient)) +
  geom_point(aes(fill = Donor_ID, size = Efficiency), alpha = 0.7, colour = "black", shape = 21) +
  scale_size_continuous(range = c(2, 8), name = "Engraftment\nefficiency (%)") +
  scale_fill_manual(values = large_colour_palette_3, name = "Donor ID") +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  coord_fixed() +
  xlab("Mean donor α-diversity (Shannon index)") + ylab("Recipient baseline α-diversity (Shannon index)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# calculate proportion of donor-recipient pairings with 0 EE for donors with mean alpha diversity < or > 2.5 
high_alpha_donors_ee <- donor_recipient_alpha_div %>%
  select(Donor_ID, Mean_shannon_donor, Recipient_ID, Shannon_recipient, Efficiency) %>%
  distinct() %>%
  filter(Mean_shannon_donor > 2.5) # 136
nrow(high_alpha_donors_ee %>% filter(Efficiency == 0)) # 61/136 pairings with 0 EE = 44.9%

low_alpha_donors_ee <- donor_recipient_alpha_div %>%
  select(Donor_ID, Mean_shannon_donor, Recipient_ID, Shannon_recipient, Efficiency) %>%
  distinct() %>%
  filter(Mean_shannon_donor < 2.5) # 37
nrow(low_alpha_donors_ee %>% filter(Efficiency == 0)) # 32/37 pairings with 0 EE = 86.5%

# test if donor alpha diversity is significantly higher than recipient baseline alpha diversity
donor_alpha_div <- donor_recipient_alpha_div %>%
  select(Donor_ID, Mean_shannon_donor) %>%
  distinct() %>%
  rename(alpha_div = Mean_shannon_donor,
         Participant_ID = Donor_ID) %>%
  mutate(Group = "Donor")
recipient_alpha_div <- donor_recipient_alpha_div %>%
  select(Recipient_ID, Shannon_recipient) %>%
  distinct() %>%
  rename(alpha_div = Shannon_recipient,
         Participant_ID = Recipient_ID) %>%
  mutate(Group = "Recipient")
single_donor_recipient_alpha_div <- bind_rows(donor_alpha_div, recipient_alpha_div)
single_donor_recipient_alpha_div %>% pull(alpha_div) %>% shapiro.test() # p = 0.067 not significantly different from normal distribution
ggplot(single_donor_recipient_alpha_div, aes(x = alpha_div)) + geom_histogram() # data skewed - above result due to small sample sizes?
wilcox.test(alpha_div ~ Group, data = single_donor_recipient_alpha_div) # use non-parametric test
# p = 0.00412

# plot donor and FMT recipient baseline alpha diversity
donor_recipient_alpha_box <- single_donor_recipient_alpha_div %>%
  ggplot(aes(x = Group, y = alpha_div, fill = Group)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2) +
  geom_boxplot(alpha = 0.7, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + # add 15% extra space above points to add significance bars if needed
  scale_fill_manual(values = c("#fdc407", "#b2182b")) +
  xlab(NULL) + ylab("α-diversity (Shannon index)") + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# plot correlation between mean donor alpha diversity and mean donor engraftment efficiency
donor_shannon_ee_correlation_test_plot <- donor_recipient_alpha_div %>%
  select(Donor_ID, Mean_shannon_donor, mean_engraftment_efficiency) %>%
  distinct() %>%
  ggplot(aes(x = Mean_shannon_donor, y = mean_engraftment_efficiency), colour = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "gray", alpha = 0.2) +
  geom_point(aes(fill = Donor_ID), alpha = 0.7, size = 4, shape = 21) +
  scale_fill_manual(values = large_colour_palette_3, name = "Donor ID") +
  xlab("Mean donor α-diversity (Shannon index)") + ylab("Mean strain engraftment\nefficiency (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        #panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# test pearson's correlation
donor_recipient_alpha_div_cor <- donor_recipient_alpha_div %>%
  select(Donor_ID, Mean_shannon_donor, mean_engraftment_efficiency) %>%
  distinct() 

ee_donoralphadiv_cor <- cor.test(donor_recipient_alpha_div_cor$Mean_shannon_donor, 
                                 donor_recipient_alpha_div_cor$mean_engraftment_efficiency) # r = -0.017, p = 0.955
ee_donoralphadiv_cor$estimate^2 # R^2 = 0.000283


## beta diversity (species)

# investigate correlation between donor-recipient species beta diversity and engraftment efficiency
# get sample data for all donors and their recipients at baseline
engraftment_efficiency_samples <- engraftment_efficiency %>%
  select(Donor_ID, Recipient_ID, Pairing, Efficiency) %>%
  pivot_longer(cols = -c(Pairing, Efficiency), 
               names_to = "Type", values_to = "Participant_ID") %>%
  mutate(Donor_ID = substr(Pairing, 1, 4)) %>%
  select(Donor_ID, Participant_ID, Efficiency) %>%
  mutate(Efficiency = ifelse(Donor_ID == Participant_ID, 100, Efficiency)) %>%
  distinct() %>%
  left_join(uc_metadata %>% filter((Group == "FMT" & Timepoint == "BL") | Group == "Donor_individual") %>%
              select(Participant_ID, Sample_ID), by = "Participant_ID", relationship = "many-to-many") %>% # join with metadata by Participant ID to get sample ID
  filter(!(Donor_ID == Participant_ID)) %>% # remove 27 rows where donor ID was duplicated in Participant ID for other plots
  rename(Recipient_sample_ID = Sample_ID, # rename columns
         Recipient_ID = Participant_ID) %>%
  left_join(uc_metadata %>% filter(Group == "Donor_individual") %>% select(Participant_ID, Sample_ID), 
            by = c("Donor_ID" = "Participant_ID"), relationship = "many-to-many") %>% # join with metadata to get individual donor ID samples
  rename(Donor_sample_ID = Sample_ID) # rename column

# calculate Bray-Curtis dissimilarity without converting to PCoA
# all individual donor and FMT baseline samples
indivdonor_recipientsBL_species <- species %>%
  rownames_to_column("Sample_ID") %>%
  left_join(uc_metadata, by = "Sample_ID") %>%
  filter((Group == "FMT" & Timepoint == "BL") | Group == "Donor_individual") %>%
  select(-c(Group, Timepoint, Participant_ID)) %>%
  column_to_rownames("Sample_ID")

beta_div <- vegdist(indivdonor_recipientsBL_species, method = "bray") 
beta_div_table <- as.data.frame(as.table(as.matrix(beta_div))) # convert matrix format to table format
colnames(beta_div_table) <- c("Sample_1", "Sample_2", "beta_diversity")
  
# get beta diversity for all true donor-recipient pairings
donor_recipient_beta_div <- beta_div_table %>%
  rename(Donor_sample_ID = Sample_1, # arbitrarily assign sample 1 as donor,
         Recipient_sample_ID = Sample_2) %>% # sample 2 as recipient
  inner_join(engraftment_efficiency_samples) %>%
  group_by(Donor_ID, Recipient_ID) %>%
  mutate(mean_beta_div = mean(beta_diversity)) %>% # average beta div for each donor-recipient pairing in cases where donor had two samples
  select(Donor_ID, Recipient_ID, Efficiency, mean_beta_div) %>% # select relevant columns
  distinct() %>% # get distinct rows
  mutate(Efficiency = na_if(Efficiency, 0)) # replace 0 effiency with NA

# plot correlation between donor-recipient species beta diversity and engraftment efficiency
donor_recipient_betadiv_ee_correlation_test_plot <- donor_recipient_beta_div %>%
  ggplot(aes(x = mean_beta_div, y = Efficiency), colour = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "gray", alpha = 0.2) +
  geom_point(aes(fill = Donor_ID), alpha = 0.7, size = 4, shape = 21) +
  scale_fill_manual(values = large_colour_palette_3, name = "Donor ID") +
  scale_x_continuous(limits = c(0.5, 1.0)) +
  scale_y_continuous(limits = c(0, 40)) +
  xlab("Mean donor-recipient species β-diversity (Bray-Curtis index)") + ylab("Strain engraftment efficiency (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        #panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# test pearson's correlation
ee_betadiv_spp_cor <- cor.test(donor_recipient_beta_div$mean_beta_div, donor_recipient_beta_div$Efficiency) # r = -0.303, p = 0.00621
ee_betadiv_spp_cor$estimate^2 # R^2 = 0.0921

length(which(donor_recipient_beta_div$mean_beta_div >= 0.6 & donor_recipient_beta_div$mean_beta_div <= 0.9)) # 113/173 pairings have 0.6-0.9 BC dissimilarity


## beta diversity (functions)

# investigate correlation between donor-recipient functions beta diversity and engraftment efficiency
# get functional data for all donors and their FMT recipients at baseline
# filter the HQ MAG data to get just MAG ID and classification
mag_spp <- 
  HQ_mags %>%
  select(MAG_ID, classification)

# join the mag species data with the genes on high quality MAGs, extract sample ID and species data
HQ_genes_spp <- 
  HQ_genes %>%
  inner_join(mag_spp) %>%
  mutate(Species = str_split(classification, "s__", simplify = TRUE)[,2]) %>% # extract species of MAG
  mutate(across(where(is.character), ~ na_if(.,""))) %>% # replace empty cells (no species) with NA
  select(-c(classification, MAG_ID)) %>%
  mutate(Sample_ID = str_split(Gene_ID, "_", simplify = TRUE)[,1]) %>% # extract sample ID from gene ID
  mutate(Species_clean = str_replace_all(Species, "_.+?", "")) %>% # '?' for non-greedy matching
  select(-c(Species, Contig_ID)) %>% # deselect columns
  rename(Species = Species_clean) # removed underscores from genera and species

# get number of HQ genes per samples
HQ_genes_count <- HQ_genes_spp %>%
  select(Sample_ID, Gene_ID) %>%
  group_by(Sample_ID) %>%
  summarise(HQ_gene_count = n()) %>%
  ungroup()

# get sample data for donors and FMT recipients at baseline
indivdonor_recipientsBL_samples <- uc_metadata %>%
  filter((Group == "FMT" & Timepoint == "BL") | Group == "Donor_individual") %>%
  select(Participant_ID, Sample_ID)

# samples with HQ genes
HQ_genes_samples <- indivdonor_recipientsBL_samples %>%
  inner_join(HQ_genes_count)

# find samples with no HQ genes
no_HQ_genes_samples <- indivdonor_recipientsBL_samples %>%
  anti_join(HQ_genes_count) %>%
  inner_join(uc_metadata) %>%
  relocate(Sample_ID, .before = Participant_ID) # 11
#write.xlsx(no_HQ_genes_samples, "no_HQgenes_samples.xlsx") # table format for supplementary information 

# filter HQ gene and MAG data (big dataset) for donors and FMT recipients at baseline
indivdonor_recipientsBL_HQ_genes <-
  HQ_genes_spp %>%
  inner_join(indivdonor_recipientsBL_samples, by = "Sample_ID") 

# filter HQ_cluster data for donors and their recipients at baseline
indivdonor_recipientsBL_HQ_clusters <- indivdonor_recipientsBL_HQ_genes %>%
  left_join(HQ_clusters, by = "Cluster_ID")

# find functional metric (GO, KEGG, COG) with highest annotations in HQ clusters across this subset
colSums(!is.na(indivdonor_recipientsBL_HQ_clusters)) # COG_Functional_cat. = 333194/386865 = 86.1%

# summarise number of genes in each sample with each functional category annotation
indivdonor_recipientsBL_cog_formatted <- indivdonor_recipientsBL_HQ_clusters %>%
  mutate(single_category = strsplit(COG_Functional_cat., "")) %>%
  unnest(single_category) %>% # separate multi-categories
  select(Sample_ID, single_category) %>%
  group_by(Sample_ID, single_category) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(HQ_genes_count, by = "Sample_ID") %>% # join with HQ gene count per sample
  mutate(function_rel_abun = count/HQ_gene_count) %>% # normalise gene function prevalence by total gene count
  select(-c(count, HQ_gene_count)) %>%
  pivot_wider(names_from = single_category, values_from = function_rel_abun) %>% # wide format
  mutate(across(everything(), ~replace_na(.x, 0))) %>% # missing functions are given 0 value
  column_to_rownames("Sample_ID")

# note that the sample rel abundance sums to slightly higher than 1 because genes with multiple function annotations were counted more than once
rowSums(indivdonor_recipientsBL_cog_formatted, 1:ncol(indivdonor_recipientsBL_cog_formatted))

# renormalise function relative abundance values to equal 1 for each sample
indivdonor_recipientsBL_cog_formatted <- indivdonor_recipientsBL_cog_formatted %>% t() %>% as.data.frame() # transpose data so samples are columns
indivdonor_recipientsBL_cog_renorm <- indivdonor_recipientsBL_cog_formatted %>% mutate(across(everything(), ~ . / sum(.))) # renormalise column values by column sums
indivdonor_recipientsBL_cog_renorm <- indivdonor_recipientsBL_cog_renorm %>% t() %>% as.data.frame() # retranspose data so samples are columns

# check row sums again, samples should sum to 1
rowSums(indivdonor_recipientsBL_cog_renorm, 1:ncol(indivdonor_recipientsBL_cog_renorm)) 

# investigate correlation between donor-recipient functions beta diversity and engraftment efficiency
# calculate Bray-Curtis dissimilarity without converting to PCoA
beta_div_functions <- vegdist(indivdonor_recipientsBL_cog_renorm, method = "bray") 
beta_div_functions_table <- as.data.frame(as.table(as.matrix(beta_div_functions))) # convert matrix format to table format
colnames(beta_div_functions_table) <- c("Sample_1", "Sample_2", "beta_diversity")

# get beta diversity for all true donor-recipient pairings
donor_recipient_beta_div_functions <- beta_div_functions_table %>%
  rename(Donor_sample_ID = Sample_1, # arbitrarily assign sample 1 as donor,
         Recipient_sample_ID = Sample_2) %>% # sample 2 as recipient
  inner_join(engraftment_efficiency_samples) %>%
  group_by(Donor_ID, Recipient_ID) %>%
  mutate(mean_beta_div = mean(beta_diversity)) %>% # average beta div for each donor-recipient pairing in cases where donor had two samples
  select(Donor_ID, Recipient_ID, Efficiency, mean_beta_div) %>% # select relevant columns
  distinct() %>% # get distinct rows
  mutate(Efficiency = na_if(Efficiency, 0)) # replace 0 effiency with NA

# note: not all samples had HQ genes, therefore some samples have missing COG function beta diversity data
# those samples will be removed from plot

# plot correlation between donor-recipient functions beta diversity and engraftment efficiency
donor_recipient_betadiv_functions_ee_correlation_test_plot <- donor_recipient_beta_div_functions %>%
  ggplot(aes(x = mean_beta_div, y = Efficiency), colour = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "gray", alpha = 0.2) +
  geom_point(aes(fill = Donor_ID), alpha = 0.7, size = 4, shape = 21) +
  scale_fill_manual(values = large_colour_palette_3[c(1, 3:14)], name = "Donor ID") +
  xlim(NA, 0.175) +
  scale_y_continuous(limits = c(0, 40)) +
  xlab("Mean donor-recipient functions β-diversity (Bray-Curtis index)") + ylab("Strain engraftment efficiency (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        #panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# test pearson's correlation
ee_betadiv_func_cor <- cor.test(donor_recipient_beta_div_functions$mean_beta_div, donor_recipient_beta_div_functions$Efficiency) # r = -0.448, p = 0.000335
ee_betadiv_func_cor$estimate^2 # R^2 = 0.200


## find the metric that is most predictive of engraftment efficiency

# engraftment efficiency data
engraftment_efficiency_samples # engraftment efficiency for true pairings (all individual donor samples, FMT recipients @ BL)

# alpha diversity
alpha_div_data <- alpha_div %>% 
  select(Sample_ID, Shannon) # alpha diversity for all samples

# species beta diversity 
beta_div_table # species beta diversity for true pairings

# functions beta diversity
beta_div_functions_table # functions beta diversity for true pairings

# Prevotella/Bacteroides ratio
pb_ratio <- genus %>%
  rownames_to_column("Sample_ID") %>% # convert row names to sample ID column
  select(Sample_ID, g__Prevotella, g__Bacteroides) %>% # select relevant columns
  mutate(Prevotella_pseudo = g__Prevotella + 1e-6, # create pseudo counts with +1e-6 to avoid zeros in calculation
         Bacteroides_pseudo = g__Bacteroides + 1e-6) %>%
  mutate(PB_ratio = Prevotella_pseudo/Bacteroides_pseudo) %>% # calculate P/B ratio
  inner_join(uc_metadata) %>% # join with meta data
  select(Sample_ID, PB_ratio) # P/B ratio for all samples

# donor batch size
batch_size <- true_pairings_formatted %>%
  select(Recipient_ID, all_donors) %>%
  distinct() %>%
  rename(Batch_size = all_donors)

# join data
ee_factors <- engraftment_efficiency_samples %>%
  # alpha diversity
  left_join(alpha_div_data, by = c("Donor_sample_ID" = "Sample_ID")) %>% # donor alpha diversity (Shannon index)
  group_by(Donor_ID) %>%
  mutate(Mean_donor_alpha_diversity = mean(Shannon)) %>%
  ungroup() %>%
  select(-Shannon) %>%
  left_join(alpha_div_data, by = c("Recipient_sample_ID" = "Sample_ID")) %>% # recipient alpha diversity (Shannon index)
  rename(Recipient_alpha_diversity = Shannon) %>%
  
  # species beta diversity
  left_join(beta_div_table, by = c("Donor_sample_ID" = "Sample_1", "Recipient_sample_ID" = "Sample_2")) %>% # donor-recipient species beta diversity 
  group_by(Donor_ID, Recipient_ID) %>%
  mutate(Mean_species_beta_diversity = mean(beta_diversity)) %>%
  ungroup() %>%
  select(-beta_diversity) %>%
  
  # Prevotella/Bacteroides ratio
  left_join(pb_ratio, by = c("Donor_sample_ID" = "Sample_ID")) %>% # donor Prevotella/Bacteroides ratio
  group_by(Donor_ID) %>%
  mutate(Mean_donor_PB_ratio = mean(PB_ratio)) %>%
  ungroup() %>%
  select(-PB_ratio) %>%
  left_join(pb_ratio, by = c("Recipient_sample_ID" = "Sample_ID")) %>% # recipient Prevotella/Bacteroides ratio
  rename(Recipient_PB_ratio = PB_ratio) %>%
  
  # functions beta diversity
  left_join(beta_div_functions_table, by = c("Donor_sample_ID" = "Sample_1", "Recipient_sample_ID" = "Sample_2")) %>% # donor-recipient functions beta diversity 
  group_by(Donor_ID, Recipient_ID) %>%
  mutate(Mean_functions_beta_diversity = mean(beta_diversity)) %>%
  ungroup() %>%
  select(-beta_diversity) %>%
  
  # calculate mean values for each donor-recipient pairing
  select(-c(Donor_sample_ID, Recipient_sample_ID)) %>%
  group_by(Donor_ID, Recipient_ID) %>%
  distinct() %>%
  ungroup() %>%
  
  # donor batch size
  left_join(batch_size, by = "Recipient_ID")

# where EE = 0, make NA
ee_factors <- ee_factors %>%
  mutate(Efficiency = na_if(Efficiency, 0))

# fit reduced linear mixed model (no interactions between fixed effects)
ee_factors_lm_model <- lmer(Efficiency ~ Mean_donor_alpha_diversity + Recipient_alpha_diversity + 
                   Mean_species_beta_diversity +
                   Mean_donor_PB_ratio + Recipient_PB_ratio + 
                   Mean_functions_beta_diversity +
                   Batch_size + (1| Donor_ID) + (1| Recipient_ID), 
                 data = ee_factors)
summary(ee_factors_lm_model)
# Mean_functions_beta_diversity p = 0.00145, Mean_species_beta_diversity = p = 0.0165

# fit full linear mixed model (with interactions between donor-recipient alpha diversity and PB ratio)
ee_factors_lm_model_full <- lmer(Efficiency ~ Mean_donor_alpha_diversity * Recipient_alpha_diversity + 
                                   Mean_species_beta_diversity +
                                   Mean_donor_PB_ratio * Recipient_PB_ratio + 
                                   Mean_functions_beta_diversity +
                                   Batch_size + (1| Donor_ID) + (1| Recipient_ID), 
                                 data = ee_factors)
summary(ee_factors_lm_model_full)
# Mean_functions_beta_diversity p = 0.00158, Mean_species_beta_diversity = p = 0.01031

# compare fit of models
anova(ee_factors_lm_model_full, ee_factors_lm_model) # p = 0.4641
# full model doesn't fit better - non-significant p value
# keep reduced model

# get confidence interval of fixed effect coefficients
confint(ee_factors_lm_model)

# extract fixed effect size
#ee_factors_fixedef <- fixef(ee_factors_lm_model_full)
ee_factors_fixedef <- fixef(ee_factors_lm_model)
ee_factors_fixedef_df <- data.frame(predictor = names(ee_factors_fixedef), effect = ee_factors_fixedef) %>%
  mutate(plot_names = case_when(predictor == "Mean_donor_alpha_diversity" ~ "Mean donor\nα-diversity", # rename predictors for plotting
                                predictor == "Recipient_alpha_diversity" ~ "Recipient\nα-diversity",
                                predictor == "Mean_species_beta_diversity" ~ "Mean species\nβ-diversity",
                                predictor == "Mean_donor_PB_ratio" ~ "Mean donor\nPrevotella/Bacteroides ratio",
                                predictor == "Recipient_PB_ratio" ~ "Recipient\nPrevotella/Bacteroides ratio",
                                predictor == "Mean_functions_beta_diversity" ~ "Mean functions\nβ-diversity",
                                predictor == "Batch_size" ~ "Batch size",
                                predictor == "(Intercept)" ~ "Intercept",
                                predictor == "Mean_donor_alpha_diversity:Recipient_alpha_diversity" ~ "Mean donor:recipient\nα diversity",
                                predictor == "Mean_donor_PB_ratio:Recipient_PB_ratio" ~ "Mean donor:recipient\nPrevotella/Bacteroides ratio")) %>%
  mutate(abs_effect = abs(effect)) %>% # get absolute effect size for plotting order
  arrange(abs_effect)
ee_factors_fixedef_df$plot_names <- factor(ee_factors_fixedef_df$plot_names, levels = ee_factors_fixedef_df$plot_names) # order bars for plotting

# plot fixed effect size
ee_factors_fixedef_plot <- ee_factors_fixedef_df %>%
  #filter(!(plot_names %in% "Intercept")) %>% # remove intercept
  ggplot(aes(x = plot_names, y = effect)) +
  geom_col() +
  coord_flip() + # plot horizontal bars
  ylab("Effect Size") + xlab(NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))


# run second linear mixed model with donor - recipient renormalised relative abundances of each COG category
# to find the effect of specific function dissimilarity on engraftment efficiency
# renormalised COG relative abundance data
cog_rel_abun <- indivdonor_recipientsBL_cog_renorm %>% # 23 categories + NA
  rownames_to_column("Sample_ID") 

functional_betadiv_factors <- engraftment_efficiency_samples %>% # donor-recipient pairings with engraftment efficiency
  left_join(cog_rel_abun, by = c("Donor_sample_ID" = "Sample_ID")) %>% # get functional data for donor samples
  rename_with(~ paste0(., "_donor"), c(C:W)) %>% # rename functional data columns for donors
  left_join(cog_rel_abun, by = c("Recipient_sample_ID" = "Sample_ID")) %>% # get functional data for recipient samples
  rename_with(~ paste0(., "_recipient"), c(C:W)) %>% # rename functional data columns for recpients
  mutate(across(ends_with("_donor"), ~ . - get(gsub("_donor", "_recipient", cur_column())), 
                .names = "{gsub('_donor', '_diff', col)}")) %>% # calculate donor-recipient difference in relative abundance
  select(Donor_ID, Recipient_ID, Efficiency, contains("diff")) %>%
  group_by(Donor_ID, Recipient_ID, Efficiency) %>% # calculate mean donor-recipient difference in relative abundance, across donor samples
  summarise(across(ends_with("_diff"), mean, .names = "mean_{col}")) %>%
  ungroup()

# fit reduced linear mixed model (no interactions between fixed effects)
ee_functions_lm_model <- lmer(Efficiency ~ 
                                mean_C_diff + mean_D_diff + mean_E_diff + mean_F_diff + mean_G_diff + mean_H_diff + mean_I_diff +
                                mean_J_diff + mean_K_diff + mean_L_diff + mean_M_diff + mean_N_diff + mean_O_diff + mean_P_diff + 
                                mean_Q_diff + mean_S_diff + mean_T_diff + mean_U_diff + mean_V_diff + mean_NA_diff + mean_B_diff + 
                                mean_Z_diff + mean_A_diff + mean_W_diff + 
                                (1| Donor_ID) + (1| Recipient_ID), data = functional_betadiv_factors)
summary(ee_functions_lm_model) 
# nothing significant


# check if sequencing depth correlates with species beta diversity

# format readcount data
n_read_counts <- kneaddata_read_counts %>%
  select(Sample, final.pair1, final.pair2, final.orphan1, final.orphan2) %>%
  mutate(total_readcounts = rowSums(across(where(is.numeric)))) %>%
  select(Sample, total_readcounts) %>%
  rename(Sample_ID = Sample)

# perform PCoA of beta diversity data
spp_beta_div_pcoa <- cmdscale(beta_div) # perform PCoA on output of vegdist()
spp_beta_div_pcoa_df <- as.data.frame(spp_beta_div_pcoa) # convert PCoA results to df format
colnames(spp_beta_div_pcoa_df) <- c("PC1", "PC2")

# rejoin pcoa data with read count per sample data
spp_beta_div_pcoa_readcount <- spp_beta_div_pcoa_df %>%
  rownames_to_column("Sample_ID") %>%
  inner_join(n_read_counts)

# plot beta diversity 
spp_beta_div_plot <- spp_beta_div_pcoa_readcount %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = total_readcounts), colour = "black", pch = 21, size = 2) +
  scale_fill_viridis(name = "Sample readcount") + # fill by sample readcount total
  xlab("PC1") + ylab("PC2") +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# data for permanova

# subset of metadata
spp_permanova_metadata <- spp_beta_div_pcoa_readcount %>%
  select(Sample_ID, total_readcounts)
nrow(spp_permanova_metadata) # 59

# subset of RA data
indivdonor_recipientsBL_species # input to vegdist

# run permanova, features: total_readcounts
# all samples
set.seed(1234)
adonis2(indivdonor_recipientsBL_species ~ total_readcounts, 
        data = spp_permanova_metadata, 
        by = "margin", # don't need to worry about order of features specified after '~' above
        permutations = 999, 
        method = "bray") # p = 0.001 ***, R2 = 0.0645


# check if sequencing depth correlates with functional beta diversity

# perform PCoA of beta diversity data
func_beta_div_pcoa <- cmdscale(beta_div_functions) # perform PCoA on output of vegdist()
func_beta_div_pcoa_df <- as.data.frame(func_beta_div_pcoa) # convert PCoA results to df format
colnames(func_beta_div_pcoa_df) <- c("PC1", "PC2")

# rejoin pcoa data with read count per sample data
func_beta_div_pcoa_readcount <- func_beta_div_pcoa_df %>%
  rownames_to_column("Sample_ID") %>%
  inner_join(n_read_counts)

# plot beta diversity 
func_beta_div_plot <- func_beta_div_pcoa_readcount %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = total_readcounts), colour = "black", pch = 21, size = 2) +
  scale_fill_viridis(name = "Sample readcount") + # fill by sample readcount total
  xlab("PC1") + ylab("PC2") +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# data for permanova

# subset of metadata
func_permanova_metadata <- func_beta_div_pcoa_readcount %>%
  select(Sample_ID, total_readcounts)
nrow(func_permanova_metadata) # 48

# subset of RA data
indivdonor_recipientsBL_cog_renorm # input to vegdist

# run permanova, features: Donor_recipient
# all recipients at BL only vs individual donor means
set.seed(1234)
adonis2(indivdonor_recipientsBL_cog_renorm ~ total_readcounts, 
        data = func_permanova_metadata, 
        by = "margin", # don't need to worry about order of features specified after '~' above
        permutations = 999, 
        method = "bray") # p = 0.019, R2 = 0.0660
