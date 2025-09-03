setwd("~/UoL_Bioinformatics/LIFE_703/Analysis/GProfiler")

# Load libraries
library(dplyr)
library(tidyr)

# Read the files
data_3dpf <- read.csv("gProfiler_acalliptera_3dpf.csv", header = TRUE)
sig_data_3dpf_BP <- data_3dpf[data_3dpf$source == "GO:BP" & data_3dpf$adjusted_p_value < 0.05, ]

data_7dpf <- read.csv("gProfiler_acalliptera_7dpf.csv", header = TRUE)
sig_data_7dpf_BP <- data_7dpf[data_7dpf$source == "GO:BP" & data_7dpf$adjusted_p_value < 0.05, ]

data_12dpf <- read.csv("gProfiler_acalliptera_12dpf.csv", header = TRUE)
sig_data_12dpf_BP <- data_12dpf[data_12dpf$source == "GO:BP" & data_12dpf$adjusted_p_value < 0.05, ]

# List of Ensembl IDs for GO:0007626 at 12dpf
loco_genes <- c("ENSACLG00000020178", "ENSACLG00000020255", "ENSACLG00000009594", "ENSACLG00000012229", 
                "ENSACLG00000019877", "ENSACLG00000008961", "ENSACLG00000009697", "ENSACLG00000002345", 
                "ENSACLG00000003960", "ENSACLG00000012029", "ENSACLG00000012432", "ENSACLG00000003503", 
                "ENSACLG00000020717", "ENSACLG00000020892", "ENSACLG00000008238", "ENSACLG00000020960", 
                "ENSACLG00000011036", "ENSACLG00000017825", "ENSACLG00000016160", "ENSACLG00000023399", 
                "ENSACLG00000021317", "ENSACLG00000000589", "ENSACLG00000022202", "ENSACLG00000005391", 
                "ENSACLG00000018100", "ENSACLG00000018790", "ENSACLG00000022562", "ENSACLG00000002267", 
                "ENSACLG00000013295", "ENSACLG00000000385", "ENSACLG00000013440", "ENSACLG00000013880", 
                "ENSACLG00000017438", "ENSACLG00000026321", "ENSACLG00000020738", "ENSACLG00000000624", 
                "ENSACLG00000022210", "ENSACLG00000024963", "ENSACLG00000015505", "ENSACLG00000025753", 
                "ENSACLG00000026171", "ENSACLG00000014846", "ENSACLG00000015371", "ENSACLG00000001872", 
                "ENSACLG00000021519", "ENSACLG00000013961", "ENSACLG00000016270", "ENSACLG00000012153", 
                "ENSACLG00000016249", "ENSACLG00000020155", "ENSACLG00000011861", "ENSACLG00000027041", 
                "ENSACLG00000005750", "ENSACLG00000022006", "ENSACLG00000026638", "ENSACLG00000008156")

# Check 3dpf: Find terms and count genes
matches_3dpf <- sig_data_3dpf_BP %>%
  filter(adjusted_p_value < 0.05) %>%
  mutate(EnsemblID = strsplit(intersections, ",")) %>%
  unnest(EnsemblID) %>%
  filter(EnsemblID %in% loco_genes) %>%
  group_by(term_id, term_name, adjusted_p_value) %>%
  summarise(Intersection_size = n(), Genes = paste(EnsemblID, collapse = ","), .groups = "drop") %>%
  filter(Intersection_size >= 10) %>%
  mutate(stage = "3dpf") %>%
  dplyr::select(term_id, term_name, adjusted_p_value, Intersection_size, Genes, stage)

# Check 7dpf: Find terms and count genes
matches_7dpf <- sig_data_7dpf_BP %>%
  filter(adjusted_p_value < 0.05) %>%
  mutate(EnsemblID = strsplit(intersections, ",")) %>%
  unnest(EnsemblID) %>%
  filter(EnsemblID %in% loco_genes) %>%
  group_by(term_id, term_name, adjusted_p_value) %>%
  summarise(Intersection_size = n(), Genes = paste(EnsemblID, collapse = ","), .groups = "drop") %>%
  filter(Intersection_size >= 10) %>%
  mutate(stage = "7dpf") %>%
  dplyr::select(term_id, term_name, adjusted_p_value, Intersection_size, Genes, stage)

# Combine results
high_overlap_terms <- bind_rows(matches_3dpf, matches_7dpf)

# Save results
write.table(high_overlap_terms, "Tables/High_overlap_above10_locomotory_terms.tsv", sep = "\t", row.names = FALSE)