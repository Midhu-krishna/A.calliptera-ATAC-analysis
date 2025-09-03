setwd("~/UoL_Bioinformatics/LIFE_703/Analysis/GProfiler")

# Load libraries
library(ggplot2)
library(dplyr)
library(ggVennDiagram)
library(tidyr)

# Read the files 
data_3dpf <- read.csv("gProfiler_acalliptera_3dpf.csv", header = TRUE)
data_7dpf <- read.csv("gProfiler_acalliptera_7dpf.csv", header = TRUE)
data_12dpf <- read.csv("gProfiler_acalliptera_12dpf.csv", header = TRUE)

# Combine data
all_data <- bind_rows(
  mutate(data_3dpf, stage = "3dpf"),
  mutate(data_7dpf, stage = "7dpf"),
  mutate(data_12dpf, stage = "12dpf")
)

write.table(all_data, "all_data_combined.tsv", sep = "\t", row.names = FALSE)

# Filter significant GO terms only
all_data_mod <- all_data[
  all_data$source %in% c("GO:MF", "GO:BP", "GO:CC") &
    all_data$adjusted_p_value < 0.05, ]

# Select only term name and stage
term_stage <- all_data_mod[, c("term_name", "stage")]
write.table(term_stage, "term_stage.tsv", sep = "\t", row.names = FALSE)

# Count in how many stages each term appears
term_counts <- term_stage %>%
  group_by(term_name) %>%
  summarise(stage_count = n_distinct(stage), .groups = "drop")

# Split based on stage_count
terms_1_stage <- term_counts %>% filter(stage_count == 1)
terms_2_stages <- term_counts %>% filter(stage_count == 2)

# Join with original to get stage information
unique_1_stage <- inner_join(terms_1_stage, term_stage, by = "term_name")
unique_2_stages <- inner_join(terms_2_stages, term_stage, by = "term_name")

# Print and save both
print("Terms unique to 1 stage:")
print(unique_1_stage)

print("Terms shared by exactly 2 stages:")
print(unique_2_stages)

write.table(unique_1_stage, "unique_terms_1stage.tsv", sep = "\t", row.names = FALSE)
write.table(unique_2_stages, "unique_terms_2stages.tsv", sep = "\t", row.names = FALSE)



# Count total number of term_name entries (including duplicates) in each stage
term_total_per_stage <- term_stage %>%
  group_by(stage) %>%
  summarise(total_terms = n(), .groups = "drop")

# Count number of unique GO term names in each stage
unique_counts_per_stage <- unique_1_stage %>%
  group_by(stage) %>%
  summarise(total_terms = n(), .groups = "drop")


# Print the counts
print(unique_counts_per_stage)
print(term_total_per_stage)