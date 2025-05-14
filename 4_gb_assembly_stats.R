## code summary:
# generate kneaddata summary data
# calculate contig counts
# isolate high quality MAGs, genes and gene clusters
# generate MAGs and gene cluster summary data
# identify phylum and species distributions of high quality MAGs
# calculate COG functional annotation coverage

## load libraries

library(tidyverse)

## running

## read counts

# process kneaddata data
kneaddata <- read_tsv("kneaddata/gutbugs_kneaddata_read_counts.tsv")

kneaddata_summary <- kneaddata %>% 
  mutate(Sample = str_remove(Sample, "_")) %>%
  mutate(Raw_reads = `raw pair1`+`raw pair2`,
         QC_reads = `final pair1`+`final pair2`+`final orphan1`+`final orphan2`,
         Human_reads = (`trimmed pair1`+`trimmed pair2`+`trimmed orphan1`+`trimmed orphan2`)-QC_reads,
         Percent_human = (Human_reads/Raw_reads)*100,
         Lost_reads = Raw_reads-QC_reads,
         Percent_lost = (Lost_reads/Raw_reads)*100) %>% 
  select(Sample_ID = Sample, Raw_reads, QC_reads, Human_reads, Lost_reads, Percent_lost, Percent_human) 

kneaddata_summary %>% 
  gather(Variable, Reads, -Sample_ID) %>% 
  group_by(Variable) %>% 
  summarise(Mean = mean(Reads),
            SD = sd(Reads),
            Min = min(Reads),
            Max = max(Reads))
# mean QC reads = 45717486 = 45.7 million
# sd QC reads = 6085116 = 6.1 million


## contig counts

contig_counts <- read.table("contigs/total_contig_counts.txt") # load contig counts for each sample

mean(contig_counts[, 1]) # mean = 96603.5
min(contig_counts[, 1]) # min = 27729
max(contig_counts[, 1]) # max = 202760


## MAG counts

# process MAGs and gene clusters data
#mags <- read_tsv("gene_catalog_mapping/MAGS.tsv", na = c("", "NA", "NaN")) %>%
#  filter(!(is.na(MAG_ID))) %>%
#  select(-`...1`) %>%
#  distinct() 

#genes <- read_tsv("gene_catalog_mapping/Individual_mapped_genes.tsv", na = c("", "NA", "NaN")) %>%
#  select(-`...1`)

#gene_clusters <- read_tsv("gene_catalog_mapping/Mapped_genes_cluster.tsv", na = c("", "NA", "NaN")) %>%
#  select(-`...1`) %>%
#  arrange(seed_eggNOG_ortholog) %>% 
#  distinct(Cluster_ID, .keep_all = T)
#save(mags, genes, gene_clusters, file = "gene-mapping-filtered.RData")

# isolate high-quality mags
#HQ_mags <- filter(mags, Completeness >90, Contamination <5)
#HQ_genes <- filter(genes, MAG_ID %in% HQ_mags$MAG_ID)
#HQ_clusters <- filter(gene_clusters, Cluster_ID %in% HQ_genes$Cluster_ID)
#save(HQ_genes, HQ_mags, HQ_clusters, file = "HQ-gene-mapping-filtered.RData")

# note: above section commented out as individual gene mapping files are zipped due to size
load("gene_catalog_mapping/gene-mapping-filtered.RData") # load filtered data
load("HQ-gene-mapping-filtered.RData") # load filtered HQ data

# grep ">" gene_catalog/nr.fa | wc -l to get genes in gene catalogue # 2,939,003 genes
n_distinct(HQ_mags$MAG_ID) # 4,189 distinct HQ MAGs
n_distinct(mags$MAG_ID) # 20,941 distinct total MAGs
n_distinct(HQ_mags$MAG_ID)/n_distinct(mags$MAG_ID)*100 # 20.0% of total MAGs are HQ
nrow(HQ_genes) # 9,425,258 genes on HQ MAGs
nrow(HQ_clusters) # 1,090,166 clusters associated with genes on HQ MAGs

# genes/MAG
HQ_genes_mag <- HQ_genes %>% group_by(MAG_ID) %>% count()
mean(HQ_genes_mag$n) # 2,250
sd(HQ_genes_mag$n) # 568
range(HQ_genes_mag$n) # 975 - 6,006

# genes/cluster
HQ_genes_cluster <- HQ_genes %>% select(Cluster_ID, Gene_ID) %>% distinct() %>% group_by(Cluster_ID) %>% count()
mean(HQ_genes_cluster$n) # 8.6
sd(HQ_genes_cluster$n) # 17.8
range(HQ_genes_cluster$n) # 1 - 331

# clusters/MAG
HQ_clusters_mag <- HQ_genes %>% select(Cluster_ID, MAG_ID) %>% distinct() %>% group_by(MAG_ID) %>% count()
mean(HQ_clusters_mag$n) # 2,238
sd(HQ_clusters_mag$n) # 563
range(HQ_clusters_mag$n) # 974 - 5,976

# MAGs/cluster
HQ_mags_cluster <- HQ_genes %>% select(Cluster_ID, MAG_ID) %>% distinct() %>% group_by(Cluster_ID) %>% count()
mean(HQ_mags_cluster$n) # 8.6
sd(HQ_mags_cluster$n) # 17.6
range(HQ_mags_cluster$n) # 1 - 238


## MAG taxa representations

# phyla distribution
mag_phyla <- 
  HQ_mags %>%
  select(MAG_ID, classification) %>%
  mutate(Phylum_tmp = str_split(classification, "p__", simplify = TRUE)[,2]) %>% 
  mutate(Phylum = str_split(Phylum_tmp, ";", simplify = TRUE)[,1]) # extract phylum of MAG

mag_phyla_counts <- data.frame(table(mag_phyla$Phylum))
mag_phyla_pct <- mag_phyla_counts %>%
  rename(Phylum = Var1,
         n_MAGs = Freq) %>%
  mutate(pct_MAGs = n_MAGs/nrow(HQ_mags)*100) %>%
  mutate(Phylum_alt = c(NA, "Actinobacteria", "Bacteroidetes", "Desulfobacterota", "Firmicutes", "Firmicutes", 
                        "Firmicutes", "Fusobacteria", "Proteobacteria", "Verrucomicrobia")) %>% # entering manually based on uniprot synonyms
  group_by(Phylum_alt) %>%
  summarise(total_n = sum(n_MAGs),
            total_pct = sum(pct_MAGs)) %>%
  arrange(-total_pct)

# species distribution
mag_species <- 
  HQ_mags %>%
  select(MAG_ID, classification) %>%
  mutate(Species = str_split(classification, "s__", simplify = TRUE)[,2]) %>% # extract species of MAG
  mutate(Species_clean = str_replace_all(Species, "_.+?", "")) %>% # '?' for non-greedy matching
  select(-Species) # deselect columns

mag_species_counts <- data.frame(table(mag_species$Species_clean))
mag_species_pct <- mag_species_counts %>%
  rename(Species = Var1,
         n_MAGs = Freq) %>%
  mutate(pct_MAGs = n_MAGs/nrow(HQ_mags)*100) %>%
  arrange(-pct_MAGs)


## functional annotations

# COG functional categories
COG_res <- HQ_clusters %>%
  filter(!is.na(COG_Functional_cat.))
nrow(COG_res)/nrow(HQ_clusters)*100 # 78.7% of HQ gene clusters have a COG functional category classification
