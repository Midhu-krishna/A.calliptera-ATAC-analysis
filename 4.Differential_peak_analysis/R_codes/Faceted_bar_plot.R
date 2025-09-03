# =====================================================
# Faceted Bar Plot of Selected GO Terms Across Stages
# =====================================================

# Setting the working directory
setwd("~/UoL_Bioinformatics/Research Project/Analysis/Differential_peak_analysis")


# load libraries
library(ggplot2)
library(tidyverse)

# 1. Define GO terms to plot
selected_terms <- c(
  "head development",
  "eye morphogenesis",
  "retina development in camera-type eye",
  "retina morphogenesis in camera-type eye"
)

# 2. Simple function to read and count genes
read_count <- function(file, stage_name) {
  df <- read_csv(file, show_col_types = FALSE)
  df <- df %>% filter(term_name %in% selected_terms)
  
  if (nrow(df) == 0) {
    # If no genes for this stage, create zero counts
    df <- data.frame(term_name = selected_terms,
                     overlap_size = 0,
                     stage = stage_name)
  } else {
    df <- df %>% 
      group_by(term_name) %>% 
      summarise(overlap_size = n()) %>% 
      mutate(stage = stage_name)
  }
  return(df)
}

# 3. Apply to all files
df_3vs7_down  <- read_count("genes_3vs7_down_BP.csv", "3vs7_down")
df_3vs12_down <- read_count("genes_3vs12_down_BP.csv", "3vs12_down")

# 4. Combine all into one table
go_data <- rbind(df_3vs7_down, df_3vs12_down)

# 6. Faceted bar plot
ggplot(go_data, aes(x = term_name, y = overlap_size, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge") +  # side-by-side bars
  facet_wrap(~ stage, ncol = 2) +                   # one facet per stage
  labs(
    title = "Selected GO Terms Across Differential Peaks",
    x = "GO Term",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # tilt x-axis labels
  )