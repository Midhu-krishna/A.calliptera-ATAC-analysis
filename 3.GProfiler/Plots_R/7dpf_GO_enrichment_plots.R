# Load libraries
library(ggplot2)

# Read the files 
data_7dpf <- read.csv("gProfiler_acalliptera_7dpf.csv", header = TRUE)

# Filter for significant terms (adjusted p-value < 0.05) and take top 20
sig_data_7dpf_MF <- data_7dpf[data_7dpf$source == "GO:MF" & data_7dpf$adjusted_p_value < 0.05, ]
sig_data_7dpf_MF <- sig_data_7dpf_MF[order(sig_data_7dpf_MF$adjusted_p_value), ][1:20, ]

sig_data_7dpf_BP <- data_7dpf[data_7dpf$source == "GO:BP" & data_7dpf$adjusted_p_value < 0.05, ]
sig_data_7dpf_BP <- sig_data_7dpf_BP[order(sig_data_7dpf_BP$adjusted_p_value), ][1:20, ]

sig_data_7dpf_CC <- data_7dpf[data_7dpf$source == "GO:CC" & data_7dpf$adjusted_p_value < 0.05, ]
sig_data_7dpf_CC <- sig_data_7dpf_CC[order(sig_data_7dpf_CC$adjusted_p_value), ][1:20, ]

# Use adjusted_p_value directly for color
sig_data_7dpf_MF$neg_log10_adj_p <- sig_data_7dpf_MF$adjusted_p_value
sig_data_7dpf_BP$neg_log10_adj_p <- sig_data_7dpf_BP$adjusted_p_value
sig_data_7dpf_CC$neg_log10_adj_p <- sig_data_7dpf_CC$adjusted_p_value

# Plot of Molecular Function
plot_7dpf_MF <- ggplot(sig_data_7dpf_MF, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                             size = intersection_size, color = neg_log10_adj_p)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "FDR") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Fold Enrichment",
    y = "GO Molecular Function",
    title = "GO Molecular Function Enrichment at 7dpf in A. calliptera"
  )

# Plot of Biological Processes
plot_7dpf_BP <- ggplot(sig_data_7dpf_BP, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                             size = intersection_size, color = neg_log10_adj_p)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "FDR") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Fold Enrichment",
    y = "GO Biological Processes",
    title = "GO Biological Processes at 7dpf in A. calliptera"
  )

# Plot of Cellular Components
plot_7dpf_CC <- ggplot(sig_data_7dpf_CC, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                             size = intersection_size, color = neg_log10_adj_p)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "FDR") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Fold Enrichment",
    y = "GO Cellular Component",
    title = "GO Cellular Component Enrichment at 7dpf in A. calliptera"
  )

# Display plots
show(plot_7dpf_MF)
show(plot_7dpf_BP)
show(plot_7dpf_CC)

# Save plots
ggsave("7dpf_go_mf_enrichment.jpeg", plot_7dpf_MF, width = 12, height = 9, dpi = 300, bg = "white")
ggsave("7dpf_go_bp_enrichment.jpeg", plot_7dpf_BP, width = 12, height = 9, dpi = 300, bg = "white")
ggsave("7dpf_go_cc_enrichment.jpeg", plot_7dpf_CC, width = 12, height = 9, dpi = 300, bg = "white")
