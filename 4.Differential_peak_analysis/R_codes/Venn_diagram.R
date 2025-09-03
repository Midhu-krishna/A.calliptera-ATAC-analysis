# ============================================================
# Venn Diagram of GO:BP Terms from Differential ATAC-seq Peaks
# ============================================================

library(ggplot2)
library(dplyr)
library(ggVennDiagram)
library(tidyr)
library(readr)

# ------------------------------------------------------------
# 1. Read the new CSV files (already filtered for BP)
# ------------------------------------------------------------
bp_3vs7_down  <- read_csv("genes_3vs7_down_BP.csv")  %>% mutate(comparison = "3vs7_down")
bp_3vs12_down <- read_csv("genes_3vs12_down_BP.csv") %>% mutate(comparison = "3vs12_down")
bp_7vs12_down <- read_csv("genes_7vs12_down_BP.csv") %>% mutate(comparison = "7vs12_down")
bp_7vs12_up   <- read_csv("genes_7vs12_up_BP.csv")   %>% mutate(comparison = "7vs12_up")

# ------------------------------------------------------------
# 2. 4-way Venn Diagram
# ------------------------------------------------------------
venn_data_4way <- list(
  `3vs7_down` = bp_3vs7_down,
  `3vs12_down` = bp_3vs12_down,
  `7vs12_down` = bp_7vs12_down,
  `7vs12_up`   = bp_7vs12_up
)

plot_venn_4way <- ggVennDiagram(venn_data_4way, label_alpha = 0) +
  scale_fill_gradient(low = "skyblue", high = "blue") +
  ggtitle("Overlap of GO:BP Terms Across Differential Peak Comparisons") +
  theme(plot.title = element_text(size = 14, face = "bold"))

show(plot_venn_4way)
ggsave("venn_4way_bp_diff_peaks.png", plot_venn_4way, width = 8, height = 8, dpi = 300, bg = "white")


# ------------------------------------------------------------
# 3. Combine all into one data frame
# ------------------------------------------------------------
bp_all <- bind_rows(bp_3vs7_down, bp_3vs12_down, bp_7vs12_down, bp_7vs12_up)

# ------------------------------------------------------------
# 4. Save combined file
# ------------------------------------------------------------
write_csv(bp_all, "combined_BP_significant_terms.csv")

