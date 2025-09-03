# Setting the working directory
setwd("~/UoL_Bioinformatics/Research Project/Analysis/Peak_annotation")

# Loading necessary libraries
library(tidyr)
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)

########################### Pie chart  ############################

# Read the file
peak_data <- read.table("1b-i.peakannot.counts2", header = FALSE, sep = "", stringsAsFactors = FALSE)
colnames(peak_data) <- c("species", "feature", "count")
head(peak_data)

# Calculate the total number of peaks to use for percentages
total_peaks <- sum(peak_data$count)

# Group by feature and calculate the sum of counts for each feature type
feature_counts <- peak_data %>%
  group_by(feature) %>%
  summarise(total_count = sum(count)) %>%
  mutate(percentage = total_count / total_peaks * 100)

# Pie chart using ggplot
ggplot(feature_counts, aes(x = "", y = percentage, fill = feature)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  scale_fill_viridis(
    discrete = TRUE,
    name = "Genomic Feature",
    labels = paste0(feature_counts$feature, " (", round(feature_counts$percentage, 3), "%)")
  ) +
  labs(title = "Distribution of Narrow Open-Chromatin Peaks")

########################### Bar Plot ##############################

# Reading the file
tss_data <- read.table("1b-ii.tss_stats.out2", sep = "\t", header = FALSE)
colnames(tss_data)[c(2, 3, 5, 8)] <- c("Sample", "Stage", "Near_TSS", "Distal_TSS")
tss_data_clean <- tss_data[, c("Sample", "Stage", "Near_TSS", "Distal_TSS")]

# Reshape the data
tss_long <- pivot_longer(
  tss_data_clean,
  cols = c(Near_TSS, Distal_TSS),
  names_to = "Peak_Type",
  values_to = "Count"
)

# Bar plot 
bar_plot <- ggplot(tss_long, aes(x = Sample, y = Count, fill = Peak_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "TSS Distribution of Accessible Peaks",
       x = "Sample", y = "Number of Peaks", fill = "Peak Type") +
  scale_fill_manual(values = c("skyblue", "green")) +
  theme_void() +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11)
  )

# Reference table
dpf_table <- data.frame(
  Sample = c("1aAc", "1bAc", "2aAc", "2bAc", "3aAc", "3bAc"),
  Stage = c("3 dpf", "3 dpf", "7 dpf", "7 dpf", "12 dpf", "12 dpf")
)

# Table grob
table_plot <- tableGrob(dpf_table, rows = NULL)

# Combine plot and table
grid.arrange(bar_plot, table_plot, ncol = 2, widths = c(3.5, 1))