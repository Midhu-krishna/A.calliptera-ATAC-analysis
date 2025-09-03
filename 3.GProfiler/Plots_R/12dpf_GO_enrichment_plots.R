# Load libraries
library(ggplot2)

# Read the files 
data_12dpf <- read.csv("gProfiler_acalliptera_12dpf.csv", header = TRUE)

# Filter for significant terms (adjusted p-value < 0.05) and take top 20
sig_data_12dpf_MF <- data_12dpf[data_12dpf$source == "GO:MF" & data_12dpf$adjusted_p_value < 0.05, ]
sig_data_12dpf_MF <- sig_data_12dpf_MF[order(sig_data_12dpf_MF$adjusted_p_value), ][1:20, ]

sig_data_12dpf_BP <- data_12dpf[data_12dpf$source == "GO:BP" & data_12dpf$adjusted_p_value < 0.05, ]
sig_data_12dpf_BP <- sig_data_12dpf_BP[order(sig_data_12dpf_BP$adjusted_p_value), ][1:20, ]

sig_data_12dpf_CC <- data_12dpf[data_12dpf$source == "GO:CC" & data_12dpf$adjusted_p_value < 0.05, ]
sig_data_12dpf_CC <- sig_data_12dpf_CC[order(sig_data_12dpf_CC$adjusted_p_value), ][1:20, ]

# Plot of Molecular Function
plot_12dpf_MF <- ggplot(sig_data_12dpf_MF, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                               size = intersection_size, color = adjusted_p_value)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "FDR (adj p)", trans = "reverse") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Fold Enrichment",
    y = "GO Molecular Function",
    title = "GO Molecular Function Enrichment at 12dpf in A. calliptera"
  )

# Plot of Biological Processes
plot_12dpf_BP <- ggplot(sig_data_12dpf_BP, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                               size = intersection_size, color = adjusted_p_value)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "FDR (adj p)", trans = "reverse") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Fold Enrichment",
    y = "GO Biological Processes",
    title = "GO Biological Processes at 12dpf in A. calliptera"
  )

# Plot of Cellular Components
plot_12dpf_CC <- ggplot(sig_data_12dpf_CC, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                               size = intersection_size, color = adjusted_p_value)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "FDR (adj p)", trans = "reverse") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Fold Enrichment",
    y = "GO Cellular Component",
    title = "GO Cellular Component Enrichment at 12dpf in A. calliptera"
  )

# Display plots
show(plot_12dpf_MF)
show(plot_12dpf_BP)
show(plot_12dpf_CC)

# Save plots
ggsave("12dpf_go_mf_enrichment.jpeg", plot_12dpf_MF, width = 12, height = 9, dpi = 300, bg = "white")
ggsave("12dpf_go_bp_enrichment.jpeg", plot_12dpf_BP, width = 12, height = 9, dpi = 300, bg = "white")
ggsave("12dpf_go_cc_enrichment.jpeg", plot_12dpf_CC, width = 12, height = 9, dpi = 300, bg = "white")
