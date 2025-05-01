# Setting the working directory
setwd("~/UoL_Bioinformatics/Research Project/Analysis")

# Loading necessary libraries
library(tidyr)
library(ggplot2)
library(viridis)
library(dplyr)

########################### Pie chart 1 ############################

# Read the file
peak_data <- read.table("Peak_Annotation/1b-i.peakannot.counts2", header = FALSE, sep = "", stringsAsFactors = FALSE)
colnames(peak_data) <- c("species", "feature", "count")
head(peak_data)

# Calculate the total number of peaks to use for percentages
total_peaks <- sum(peak_data$count)

# Group by feature and calculate the sum of counts for each feature type
feature_counts <- peak_data %>%
  group_by(feature) %>%
  summarise(total_count = sum(count)) %>%
  mutate(percentage = total_count / total_peaks * 100)

# Create a pie chart using ggplot
ggplot(feature_counts, aes(x = "", y = percentage, fill = feature)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  scale_fill_viridis(discrete = TRUE) +  # Using viridis colors
  labs(title = "Distribution of Narrow Open-Chromatin Peaks")


########################### Bar Plot ##############################

# Reading the file
tss_data <- read.table("1b-ii.tss_stats.out2", sep = "\t", header = FALSE)
head(tss_data)

# Rename relevant columns
colnames(tss_data)[c(2, 3, 5, 8)] <- c("Sample", "Stage", "Near_TSS", "Distal_TSS")

# Keep only those useful columns
tss_data_clean <- tss_data[, c("Sample", "Stage", "Near_TSS", "Distal_TSS")]

# Reshape the data
tss_long <- pivot_longer(
  tss_data_clean,
  cols = c(Near_TSS, Distal_TSS),
  names_to = "Peak_Type",
  values_to = "Count"
)

# Bar plot
ggplot(tss_long, aes(x = Sample, y = Count, fill = Peak_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "TSS Distribution of Accessible Peaks",
       x = "Sample", y = "Number of Peaks", fill = "Peak Type") +
  scale_fill_manual(values = c("skyblue", "green")) +
  theme_void() +
  theme_minimal()

