# Setting the working directory
setwd("~/UoL_Bioinformatics/Research Project/Analysis/GAT")


# Load required libraries
library(ggplot2)
library(dplyr)

# Read the three files and add a "stage" column
data_3dpf <- read.table("gatnormed_Ac_3dpf_Peak_AllAnnot_allIDR_ChrBgd.tsv", header = TRUE, sep = "\t")
data_3dpf$stage <- "3dpf"

data_7dpf <- read.table("gatnormed_Ac_7dpf_Peak_AllAnnot_allIDR_ChrBgd.tsv", header = TRUE, sep = "\t")
data_7dpf$stage <- "7dpf"

data_12dpf <- read.table("gatnormed_Ac_12dpf_Peak_AllAnnot_allIDR_ChrBgd.tsv", header = TRUE, sep = "\t")
data_12dpf$stage <- "12dpf"

# Keep only significant results (qvalue < 0.05)
sig_data_3dpf <- data_3dpf[data_3dpf$qvalue < 0.05, ]
sig_data_7dpf <- data_7dpf[data_7dpf$qvalue < 0.05, ]
sig_data_12dpf <- data_12dpf[data_12dpf$qvalue < 0.05, ]

# Calculate -log10(qvalue)
sig_data_3dpf$log_qvalue <- -log10(sig_data_3dpf$qvalue)
sig_data_7dpf$log_qvalue <- -log10(sig_data_7dpf$qvalue)
sig_data_12dpf$log_qvalue <- -log10(sig_data_12dpf$qvalue)

# Dotplot

#3dpf
plot_3dpf <- ggplot(sig_data_3dpf, aes(x = fold, y = annotation, size = overlap_nsegments, color = qvalue)) +
  geom_point() +  # Draw circles
  scale_size(name = "Peak Count", range = c(2, 10)) +  # Size for peak count
  scale_color_gradient(low = "blue", high = "red", name = "qvalue", limits = c(0, 0.05)) +  # Linear color scale for qvalue
  theme_minimal() +  # Simple theme
  theme(
    axis.text = element_text(size = 10),  # Text size
    axis.title = element_text(size = 12),  # Title size
    plot.title = element_text(size = 12, face = "bold")  # Plot title
  ) +
  labs(
    x = "Fold Enrichment",  # X-axis label
    y = "Genomic Region",  # Y-axis label
    title = "Peak Enrichment at 3dpf in A. calliptera"
  )
show(plot_3dpf)

#7dpf
plot_7dpf <- ggplot(sig_data_7dpf, aes(x = fold, y = annotation, size = overlap_nsegments, color = qvalue)) +
  geom_point() +  # Draw circles
  scale_size(name = "Peak Count", range = c(2, 10)) +  # Size for peak count
  scale_color_gradient(low = "blue", high = "red", name = "qvalue", limits = c(0, 0.05)) +  # Linear color scale for qvalue
  theme_minimal() +  # Simple theme
  theme(
    axis.text = element_text(size = 10),  # Text size
    axis.title = element_text(size = 12),  # Title size
    plot.title = element_text(size = 12, face = "bold")  # Plot title
  ) +
  labs(
    x = "Fold Enrichment",  # X-axis label
    y = "Genomic Region",  # Y-axis label
    title = "Peak Enrichment at 7dpf in A. calliptera"
  )
show(plot_7dpf)

#12dpf
plot_12dpf <- ggplot(sig_data_12dpf, aes(x = fold, y = annotation, size = overlap_nsegments, color = qvalue)) +
  geom_point() +  # Draw circles
  scale_size(name = "Peak Count", range = c(2, 10)) +  # Size for peak count
  scale_color_gradient(low = "blue", high = "red", name = "qvalue", limits = c(0, 0.05)) +  # Linear color scale for qvalue
  theme_minimal() +  # Simple theme
  theme(
    axis.text = element_text(size = 10),  # Text size
    axis.title = element_text(size = 12),  # Title size
    plot.title = element_text(size = 12, face = "bold")  # Plot title
  ) +
  labs(
    x = "Fold Enrichment",  # X-axis label
    y = "Genomic Region",  # Y-axis label
    title = "Peak Enrichment at 12dpf in A. calliptera"
  )
show(plot_12dpf)