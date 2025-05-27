## code summary:
# read in sample data and meta data
# format data
# select samples for donors, recipients at baseline and week 8 for assembly

## load libraries

library(readxl)
library(tidyverse)
library(openxlsx)

## load data

uc_metadata <- read_excel(path = "metadata.xlsx", sheet = "All_formatting") # FOCUS study sample metadata
uc_sampledata <- read.table("PRJEB26357.runinfo_ftp.tsv", header = TRUE, sep = "\t") # data download (nf-core) sample metadata

## running

# checking sample count discrepancies
anti_join(uc_metadata, uc_sampledata, by = c("Sample_ID" = "sample_alias")) %>% select(Sample_ID)
# present in meta data but not samples:
# MET102   

anti_join(uc_sampledata, uc_metadata, by = c("sample_alias" = "Sample_ID")) %>% select(sample_alias)
# present in samples but not meta data:
# 0 samples

# subset for available data
uc_metadata_available <- 
  uc_sampledata %>%
  select(sample_alias, id, single_end) %>%
  left_join(uc_metadata, by = c("sample_alias" = "Sample_ID"))

# format meta data
uc_metadata_formatted <- 
  uc_metadata_available %>%
  group_by(Participant_ID) %>%
  mutate(Group = ifelse(
    is.na(Group) & any(grepl("P8", Timepoint, fixed = TRUE)), "Placebo", # recipients with P8 timepoint are placebo recipients
    ifelse(is.na(Group) & !any(grepl("P8", Timepoint, fixed = TRUE)), "FMT", # recipients with no P8 timepoints are FMT recipients
           Group)
  )) %>%
  ungroup()

# check which have single ends (cloud pipeline uses paired ends)
single_ends <-
  uc_metadata_formatted %>%
  filter(single_end == "true")
# 0 samples

## sample selection

# select only Donor, Baseline and (blinded) Week 8 samples at first
uc_samples_firstrun <- 
  uc_metadata_formatted %>%
  filter(Timepoint %in% c("Donor", "Tx0", "P8") | (Timepoint == "Tx8" & Group == "FMT")) %>%
  filter(single_end == "false") %>% # only paired end samples
  select(id) %>%
  mutate(one = paste0(id, "_1.fastq.gz"),
         two = paste0(id, "_2.fastq.gz"))

uc_samples_forward <-
  uc_samples_firstrun %>%
  pull(one)

uc_samples_reverse <- 
  uc_samples_firstrun %>%
  pull(two)

uc_output_data <- data.frame(c(uc_samples_forward, uc_samples_reverse))

# save subsetted sample data with fastq suffix for assembly
#write.table(uc_output_data, "uc_samples_subset.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# get metadata for subsetted samples only
uc_metadata <- 
  uc_metadata_formatted %>%
  filter(Timepoint %in% c("Donor", "Tx0", "P8") | (Timepoint == "Tx8" & Group == "FMT")) %>%
  filter(single_end == "false") %>%
  select(-c(single_end, sample_alias)) %>%
  rename(Sample_ID = id) %>%
  mutate(Sample_ID = str_remove(Sample_ID, "_"),
         Timepoint = case_when(Timepoint == "Tx0" ~ "BL",
                               Timepoint == "Tx8" ~ "wk8",
                               Timepoint == "P8" ~ "wk8",
                               Timepoint == "Donor" ~ "Donor"))

# save metadata
#write.xlsx(uc_metadata, "uc_metadata_subset.xlsx") # table format for supplementary information
#save(uc_metadata, file = "uc_metadata_subset.RData")

# count participants
uc_metadata %>% group_by(Group, Timepoint) %>% summarise(n = n())

# count recipients with both timepoints
uc_FMT_both_timepoints <- uc_metadata %>%
  filter(Group == "FMT") %>%
  group_by(Participant_ID) %>%
  summarize(has_t1 = any(Timepoint == "BL"),
            has_t2 = any(Timepoint == "wk8")) %>%
  filter(has_t1 & has_t2) %>%
  pull(Participant_ID) # 32

uc_placebo_both_timepoints <- uc_metadata %>%
  filter(Group == "Placebo") %>%
  group_by(Participant_ID) %>%
  summarize(has_t1 = any(Timepoint == "BL"),
            has_t2 = any(Timepoint == "wk8")) %>%
  filter(has_t1 & has_t2) %>%
  pull(Participant_ID) # 20

# samples with low QC read counts that did not succeed binning on Terra
failed_binning_samples <- c("ERX2605398_ERR2589122",
                            "ERX2605448_ERR2589172",
                            "ERX2605450_ERR2589174",
                            "ERX2605461_ERR2589185",
                            "ERX2605467_ERR2589191",
                            "ERX2605511_ERR2589235",
                            "ERX2605563_ERR2589287",
                            "ERX2605597_ERR2589321")

failed_binning_samples <- as.data.frame(failed_binning_samples)

failed_binning_samples_metadata <- failed_binning_samples %>%
  mutate(Sample_ID = str_remove(failed_binning_samples, "_")) %>%
  select(-failed_binning_samples) %>%
  inner_join(uc_metadata)
#write.xlsx(failed_binning_samples_metadata, "failed_binning_samples.xlsx") # table format for supplementary information 
