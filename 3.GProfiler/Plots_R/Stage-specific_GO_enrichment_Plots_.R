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

# Cross-Stage Comparison with Faceted Plot
top_terms <- all_data %>%
  filter(source == "GO:BP", adjusted_p_value < 0.05) %>%
  group_by(stage) %>%
  slice_min(adjusted_p_value, n = 30) %>%
  ungroup()

plot_fold_faceted <- ggplot(top_terms, aes(x = stage, y = intersection_size, group = term_name, color = term_name)) +
  geom_line() +
  geom_point() +
  facet_wrap(~term_name, scales = "free_y", ncol = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(x = "Stage",
       y = "Number of Overlapping Genes",
       title = "Overlap Size of Top GO Terms Across Stages",
       color = "GO Term")

show(plot_fold_faceted)
ggsave("overlap_size_faceted.png", plot_fold_faceted, width = 12, height = 9, dpi = 300, bg = "white")

# Venn Diagram
stages <- c("3dpf", "7dpf", "12dpf")
venn_data <- lapply(stages, function(s) {
  all_data %>% filter(stage == s) %>% pull(term_name)
})
names(venn_data) <- stages
plot_venn <- ggVennDiagram(venn_data, label_alpha = 0) +
  scale_fill_gradient(low = "skyblue", high = "blue")+
  ggtitle("Overlap of GO Terms Across Stages")

show(plot_venn)
ggsave("go_term_venn.png", plot_venn, width = 8, height = 8, dpi = 300, bg = "white")


#######################################################################################################



# Filter for stage-specific GO Biological Process terms
bp_data <- all_data %>%
  filter(source == "GO:BP") %>%
  mutate(stage = factor(stage, levels = c("3dpf", "7dpf", "12dpf")))

plot_bp_faceted <- ggplot(bp_data, aes(x = intersection_size, 
                                       y = reorder(term_name, -intersection_size), 
                                       size = intersection_size, 
                                       color = adjusted_p_value)) +
  geom_point(alpha = 0.6) +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "blue", name = "FDR (adj. p-value)", trans = "log10") +
  facet_wrap(~stage, scales = "free_y", ncol = 1, switch = "y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(x = "Number of Overlapping Genes",
       y = "GO Biological Process",
       title = "Stage-Specific GO Biological Process Enrichment in A. calliptera")

# Display and save
show(plot_bp_faceted)
ggsave("3dpf_7dpf_12dpf_go_bp_faceted.png", plot_bp_faceted, width = 12, height = 9, dpi = 300, bg = "white")

