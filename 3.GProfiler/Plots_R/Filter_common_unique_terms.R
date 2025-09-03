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

# Filter MF GO terms only
MF_data <- all_data[
  all_data$source %in% c("GO:BP") &
    all_data$adjusted_p_value < 0.05, ]

# Filter to include stage and q value
term_stage <- MF_data[, c("source", "stage", "term_name", "adjusted_p_value")]

# Seperate the data by stage
dpf3 <- term_stage %>% filter(stage == "3dpf")
dpf7 <- term_stage %>% filter(stage == "7dpf")
dpf12 <- term_stage %>% filter(stage == "12dpf")

# Create list of term names for comparison
terms_3dpf <- dpf3$term_name
terms_7dpf <- dpf7$term_name
terms_12dpf <- dpf12$term_name

############ Comparing the three stages #########

# 3dpf vs 7dpf
common_3vs7 <- intersect(terms_3dpf, terms_7dpf)
unique_3v7 <- setdiff(terms_3dpf, terms_7dpf)
unique_7v3 <- setdiff(terms_7dpf, terms_3dpf)

commonterms_3v7 <- term_stage %>%
  filter(term_name %in% common_3vs7 & stage %in% c("3dpf", "7dpf"))

uniqueterms_3v7 <- term_stage %>%
  filter((term_name %in% unique_3v7 & stage == "3dpf") |
  (term_name %in% unique_7v3 & stage == "7dpf"))


# 3dpf vs 12dpf
common_3vs12 <- intersect(terms_3dpf, terms_12dpf)
unique_3v12 <- setdiff(terms_3dpf, terms_12dpf)
unique_12v3 <- setdiff(terms_12dpf, terms_3dpf)

commonterms_3v12 <- term_stage %>%
  filter(term_name %in% common_3vs12 & stage %in% c("3dpf", "12dpf"))

uniqueterms_3v12 <- term_stage %>%
  filter((term_name %in% unique_3v12 & stage == "3dpf") |
           (term_name %in% unique_12v3 & stage == "12dpf"))

# 7dpf vs 12dpf
common_7vs12 <- intersect(terms_7dpf, terms_12dpf)
unique_7v12 <- setdiff(terms_7dpf, terms_12dpf)
unique_12v7 <- setdiff(terms_12dpf, terms_7dpf)

commonterms_7v12 <- term_stage %>%
  filter(term_name %in% common_7vs12 & stage %in% c("7dpf", "12dpf"))

uniqueterms_7v12 <- term_stage %>%
  filter((term_name %in% unique_7v12 & stage == "7dpf") |
           (term_name %in% unique_12v7 & stage == "12dpf"))

# Terms common in all three stages
common_all <- intersect(intersect(terms_3dpf, terms_7dpf), terms_12dpf)

common_all_stages <- term_stage %>%
  filter(term_name %in% common_all)

# Terms unique to one stage only
unique_all <- term_stage %>%
  group_by(term_name) %>%
  filter(n_distinct(stage) == 1) %>%
  ungroup()

# Save results as seperate tables
write.table(unique_all, "Tables/Unique_all_stages.tsv", sep = "\t", row.names = FALSE)
write.table(common_all_stages, "Tables/Common_all_stages.tsv", sep = "\t", row.names = FALSE)

write.table(commonterms_3v7, "Tables/Commonterms_3v7.tsv", sep = "\t", row.names = FALSE)
write.table(commonterms_3v12, "Tables/Commonterms_3v12.tsv", sep = "\t", row.names = FALSE)
write.table(commonterms_7v12, "Tables/Commonterms_7v12.tsv", sep = "\t", row.names = FALSE)

write.table(uniqueterms_3v7, "Tables/Uniqueterms_3v7.tsv", sep = "\t", row.names = FALSE)
write.table(uniqueterms_3v12, "Tables/Uniqueterms_3v12.tsv", sep = "\t", row.names = FALSE)
write.table(uniqueterms_7v12, "Tables/Uniqueterms_7v12.tsv", sep = "\t", row.names = FALSE)

