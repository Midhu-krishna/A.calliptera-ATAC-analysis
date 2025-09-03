# ------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------
library(ggplot2)

# ------------------------------------------------------------
# 2. Load differential peak results
# ------------------------------------------------------------
data_3vs7 <- read.table("differential_peaks_3vs7.tsv", header = TRUE)
data_3vs12 <- read.table("differential_peaks_3vs12.tsv", header = TRUE)
data_7vs12 <- read.table("differential_peaks_7vs12.tsv", header = TRUE)

# ------------------------------------------------------------
# 3. Function to generate volcano and MA plots
# ------------------------------------------------------------
plot_volcano_ma <- function(data, contrast_name) {
  # Mark significant peaks
  data$sig <- data$FDR < 0.05
  
  # ---- Volcano plot ----
  p_volcano <- ggplot(data, aes(x = logFC, y = -log10(FDR))) +
    geom_point(aes(color = sig), size = 1) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
    labs(title = paste("Volcano:", contrast_name),
         x = "log2 Fold Change", y = "-log10(FDR)") +
    theme_classic() +  # White background
    theme(legend.position = "none")
  
  ggsave(filename = paste0("volcano_", contrast_name, ".jpeg"),
         plot = p_volcano, width = 6, height = 4, units = "in", dpi = 300)
  
  # ---- MA plot ----
  p_ma <- ggplot(data, aes(x = logCPM, y = logFC)) +
    geom_point(aes(color = sig), size = 1) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
    labs(title = paste("MA:", contrast_name),
         x = "logCPM", y = "log2 Fold Change") +
    theme_classic() +  # White background
    theme(legend.position = "none")
  
  ggsave(filename = paste0("ma_", contrast_name, ".jpeg"),
         plot = p_ma, width = 6, height = 4, units = "in", dpi = 300)
}

# ------------------------------------------------------------
# 4. Generate plots
# ------------------------------------------------------------
plot_volcano_ma(data_3vs7, "3vs7")
plot_volcano_ma(data_3vs12, "3vs12")
plot_volcano_ma(data_7vs12, "7vs12")
