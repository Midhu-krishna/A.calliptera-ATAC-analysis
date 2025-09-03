# ============================================================
# GO:BP Term Overlap Analysis Across Developmental Stages
# ============================================================
# Description:
# This script generates Venn diagrams for GO Biological Process (GO:BP) terms 
# enriched in ATAC-seq differential peaks across 3dpf, 7dpf, and 12dpf stages 
# in Astatotilapia calliptera. Both 3-way and pairwise overlaps are visualized 
# using ggVennDiagram, with color gradients indicating the number of overlapping terms.
# ============================================================

setwd("~/UoL_Bioinformatics/LIFE_703/Analysis/GProfiler")

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggVennDiagram)
library(tidyr)

# ------------------------------------------------------------
# 1. Read G:Profiler CSV files for each stage
# ------------------------------------------------------------
data_3dpf <- read.csv("gProfiler_acalliptera_3dpf.csv", header = TRUE)
data_7dpf <- read.csv("gProfiler_acalliptera_7dpf.csv", header = TRUE)
data_12dpf <- read.csv("gProfiler_acalliptera_12dpf.csv", header = TRUE)

# ------------------------------------------------------------
# 2. Combine data with stage labels
# ------------------------------------------------------------
all_data <- bind_rows(
  mutate(data_3dpf, stage = "3dpf"),
  mutate(data_7dpf, stage = "7dpf"),
  mutate(data_12dpf, stage = "12dpf")
)

# ------------------------------------------------------------
# 3. Filter for significant GO:BP terms
# ------------------------------------------------------------
bp_data <- all_data %>%
  filter(source == "GO:BP", adjusted_p_value < 0.05)

# ------------------------------------------------------------
# 4. Extract GO:BP term names for each stage
# ------------------------------------------------------------
terms_3dpf <- bp_data %>% filter(stage == "3dpf") %>% pull(term_name)
terms_7dpf <- bp_data %>% filter(stage == "7dpf") %>% pull(term_name)
terms_12dpf <- bp_data %>% filter(stage == "12dpf") %>% pull(term_name)

# ------------------------------------------------------------
# 5. 3-Way Venn Diagram
# ------------------------------------------------------------
venn_data_3way <- list(
  `3dpf` = terms_3dpf,
  `7dpf` = terms_7dpf,
  `12dpf` = terms_12dpf
)

plot_venn_3way <- ggVennDiagram(venn_data_3way, label_alpha = 0) +
  scale_fill_gradient(low = "skyblue", high = "blue") +
  ggtitle("3-Way Overlap of GO:BP Terms Across 3dpf, 7dpf, and 12dpf") +
  theme(plot.title = element_text(size = 14, face = "bold"))

show(plot_venn_3way)
ggsave("venn_3way_go_bp.png", plot_venn_3way, width = 8, height = 8, dpi = 300, bg = "white")

# ------------------------------------------------------------
# 6. Pairwise Venn Diagrams
# ------------------------------------------------------------

# 3dpf vs 7dpf
venn_data_3v7 <- list(
  `3dpf` = terms_3dpf,
  `7dpf` = terms_7dpf
)
plot_venn_3v7 <- ggVennDiagram(venn_data_3v7, label_alpha = 0) +
  scale_fill_gradient(low = "skyblue", high = "blue") +
  ggtitle("Pairwise Overlap of GO:BP Terms: 3dpf vs 7dpf") +
  theme(plot.title = element_text(size = 14, face = "bold"))

show(plot_venn_3v7)
ggsave("venn_3v7_go_bp.png", plot_venn_3v7, width = 8, height = 8, dpi = 300, bg = "white")

# 3dpf vs 12dpf
venn_data_3v12 <- list(
  `3dpf` = terms_3dpf,
  `12dpf` = terms_12dpf
)
plot_venn_3v12 <- ggVennDiagram(venn_data_3v12, label_alpha = 0) +
  scale_fill_gradient(low = "skyblue", high = "blue") +
  ggtitle("Pairwise Overlap of GO:BP Terms: 3dpf vs 12dpf") +
  theme(plot.title = element_text(size = 14, face = "bold"))

show(plot_venn_3v12)
ggsave("venn_3v12_go_bp.png", plot_venn_3v12, width = 8, height = 8, dpi = 300, bg = "white")

# 7dpf vs 12dpf
venn_data_7v12 <- list(
  `7dpf` = terms_7dpf,
  `12dpf` = terms_12dpf
)
plot_venn_7v12 <- ggVennDiagram(venn_data_7v12, label_alpha = 0) +
  scale_fill_gradient(low = "skyblue", high = "blue") +
  ggtitle("Pairwise Overlap of GO:BP Terms: 7dpf vs 12dpf") +
  theme(plot.title = element_text(size = 14, face = "bold"))

show(plot_venn_7v12)
ggsave("venn_7v12_go_bp.png", plot_venn_7v12, width = 8, height = 8, dpi = 300, bg = "white")
