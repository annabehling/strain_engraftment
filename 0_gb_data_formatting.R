## code summary:
# format Gut Bugs Trial data necessary to repeat analysis performed on FOCUS study
# reduced metadata for donors, baseline and week 6 FMT and placebo samples
# FMT donor-recipient pairings and presence/absence matrix

## load libraries

library(readxl)
library(tidyverse)
library(ape)
library(phangorn)

## load data

gb_metadata_batches <- read_excel(path = "metadata_batches.xlsx", col_names = TRUE) # true pairings

## running

# format metadata

gb_metadata <- gb_metadata_batches %>% # make reduced metadata for donors, baseline and week 6 FMT and placebo samples
  select(Sample_ID, Participant_ID, Group, Timepoint) %>%
  filter(Timepoint == "Donor" | Timepoint == "Baseline" | Timepoint == "Week 6") %>%
  mutate(Sample_ID = str_replace(Sample_ID, "6WK", "6wk"), # fix formatting discrepancies in sample IDs
         Sample_ID = str_replace(Sample_ID, "DF16a_D4", "DF16_D4a"), # these were formatted differently in meta data vs strain data
         Sample_ID = str_replace(Sample_ID, "DF16b_D4", "DF16_D4b"),
         Sample_ID = str_replace(Sample_ID, "DF17a_D2", "DF17_D2a"),
         Sample_ID = str_replace(Sample_ID, "DF17b_D2", "DF17_D2b")) %>%
  mutate(Timepoint = str_replace(Timepoint, "Baseline", "BL"), # rename timepoints to match FOCUS formatting
         Timepoint = str_replace(Timepoint, "Week 6", "wk6")) 

# save metadata
#write.xlsx(gb_metadata_subset.xlsx") # table format for supplementary information
#save(gb_metadata_subset.RData")

# count participants
gb_metadata %>% group_by(Group, Timepoint) %>% summarise(n = n())

# count recipients with both timepoints
gb_FMT_both_timepoints <- gb_metadata %>%
  filter(Group == "FMT") %>%
  group_by(Participant_ID) %>%
  summarize(has_t1 = any(Timepoint == "BL"),
            has_t2 = any(Timepoint == "wk6")) %>%
  filter(has_t1 & has_t2) %>%
  pull(Participant_ID) # 39

gb_placebo_both_timepoints <- gb_metadata %>%
  filter(Group == "Placebo") %>%
  group_by(Participant_ID) %>%
  summarize(has_t1 = any(Timepoint == "BL"),
            has_t2 = any(Timepoint == "wk6")) %>%
  filter(has_t1 & has_t2) %>%
  pull(Participant_ID) # 44

# format FMT batch data

donor_batches <- gb_metadata_batches %>%
  filter(Group == "Donor") %>%
  select(Participant_ID, Capsule_batch) %>%
  distinct() %>%
  rename(Donor_ID = Participant_ID)

FMT_batches <- gb_metadata_batches %>%
  filter(Group == "FMT" & Timepoint == "Week 6") %>%
  select(Participant_ID, Capsule_batch) %>%
  distinct() %>%
  rename(Recipient_ID = Participant_ID)

donor_recipient_pairings <- donor_batches %>% # FMT donor-recipient pairings
  left_join(FMT_batches, by = "Capsule_batch", relationship = "many-to-many") %>%
  select(-Capsule_batch) # participant IDs don't require formatting
#save(donor_recipient_pairings, file = "gb_pairings.RData")

# format FMT batch data as presence absence matrix
donor_recipient_pairings_matrix <- donor_recipient_pairings %>%
  mutate(value = 1) %>%
  arrange(Recipient_ID) %>% # arrange recipient IDs (will be the new column headers)
  pivot_wider(names_from = Recipient_ID, values_from = value, values_fill = 0) %>% # convert data to wide format
  arrange(Donor_ID) 

# transpose data so donor IDs are column names
gb_pairings <- donor_recipient_pairings_matrix %>% 
  column_to_rownames("Donor_ID") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("...1")
#save(gb_pairings, file = "gb_pairings.RData")


# format strainphlan data

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
# strain data for more timepoints than we need 
# produces strain DNA dist measures for just the samples in metadata
strain_dna_dist <- strain_dna_dist %>% 
  mutate(Sample_ID = Sample1) %>% 
  left_join(gb_metadata) %>% 
  mutate(Sample1_details = paste0(Participant_ID, "_", Group, "_", Timepoint)) %>% 
  select(-Sample_ID, -Participant_ID, -Group, -Timepoint) %>% 
  mutate(Sample_ID = Sample2) %>% 
  left_join(gb_metadata) %>% 
  mutate(Sample2_details = paste0(Participant_ID, "_", Group, "_", Timepoint)) %>% 
  select(-Sample_ID, -Participant_ID, -Group, -Timepoint) %>%
  filter(!((Sample1_details == "NA_NA_NA") | (Sample2_details == "NA_NA_NA"))) # remove matches with Sample1 or Sample2 NAs

# summary descriptions (note: captures all timepoints out to 4 years)
sum(strain_counts$Total)           # 10,799 individual strains profiles 
n_distinct(sample_strains$Species) # belonging to 110 species

strain_counts_sample <- sample_strains %>% 
  group_by(Sample_ID) %>% 
  count() 

mean(strain_counts_sample$n) # average of 25 strains identified/sample
range(strain_counts_sample$n) # minimum of 7 - maximum of 39

# save data
#save(strain_dna_dist, strain_counts, sample_strains, strain_counts_sample, file = "strainphlan.Rdata")


# format metaphlan3 data from previous study

# already renormalised metaphlan3 relative abundance data
load("metaphlan3/species.Rdata")
species <- species %>%
  rownames_to_column("Sample_ID") %>% # row names to column
  right_join(gb_metadata %>% select(Sample_ID), by = "Sample_ID") %>% # join RA data with metadata to subset
  column_to_rownames("Sample_ID")

load("metaphlan3/genus.Rdata")
genus <- genus %>%
  rownames_to_column("Sample_ID") %>% # row names to column
  right_join(gb_metadata %>% select(Sample_ID), by = "Sample_ID") %>% # join RA data with metadata to subset
  column_to_rownames("Sample_ID")

load("metaphlan3/family.Rdata")
family <- family %>%
  rownames_to_column("Sample_ID") %>% # row names to column
  right_join(gb_metadata %>% select(Sample_ID), by = "Sample_ID") %>% # join RA data with metadata to subset
  column_to_rownames("Sample_ID")

load("metaphlan3/order.Rdata")
order <- order %>%
  rownames_to_column("Sample_ID") %>% # row names to column
  right_join(gb_metadata %>% select(Sample_ID), by = "Sample_ID") %>% # join RA data with metadata to subset
  column_to_rownames("Sample_ID")

load("metaphlan3/phylum.Rdata")
phylum <- phylum %>%
  rownames_to_column("Sample_ID") %>% # row names to column
  right_join(gb_metadata %>% select(Sample_ID), by = "Sample_ID") %>% # join RA data with metadata to subset
  column_to_rownames("Sample_ID")

# save all as one Rdata object
#save(species, genus, family, order, phylum, file = "metaphlan.Rdata")