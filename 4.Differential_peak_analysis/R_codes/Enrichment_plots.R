# ============================================================
# G:Profiler GO Enrichment Plotting for Differential ATAC-seq Peaks 
# ============================================================

library(ggplot2)

# ============================================================
# 1. 3vs7 dpf comparison
# ------------------------------------------------------------
data_3vs7 <- read.csv("genes_3vs7_down_BP.csv", header = TRUE)
# Select top 20 by adjusted p-value
data_3vs7 <- data_3vs7[order(data_3vs7$adjusted_p_value), ][1:20, ]
data_3vs7$neg_log10_adj_p <- -log10(data_3vs7$adjusted_p_value)

plot_3vs7_BP <- ggplot(data_3vs7, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                      size = intersection_size, color = neg_log10_adj_p)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(x = "Number of Genes", y = "GO Biological Process", title = "GO Biological Process Enrichment 3vs7 dpf (Down, Top 20)")

print(plot_3vs7_BP)
ggsave("3vs7_BP_enrichment.jpeg", plot_3vs7_BP, width = 12, height = 9, dpi = 300, bg = "white")

# ============================================================
# 2. 3vs12 dpf comparison
# ------------------------------------------------------------
data_3vs12 <- read.csv("genes_3vs12_down_BP.csv", header = TRUE)
data_3vs12 <- data_3vs12[order(data_3vs12$adjusted_p_value), ][1:20, ]
data_3vs12$neg_log10_adj_p <- -log10(data_3vs12$adjusted_p_value)

plot_3vs12_BP <- ggplot(data_3vs12, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                        size = intersection_size, color = neg_log10_adj_p)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(x = "Number of Genes", y = "GO Biological Process", title = "GO Biological Process Enrichment 3vs12 dpf (Down, Top 20)")

print(plot_3vs12_BP)
ggsave("3vs12_BP_enrichment.jpeg", plot_3vs12_BP, width = 12, height = 9, dpi = 300, bg = "white")

# ============================================================
# 3. 7vs12 dpf comparison (Down)
# ------------------------------------------------------------
data_7vs12_down <- read.csv("genes_7vs12_down_BP.csv", header = TRUE)
data_7vs12_down$neg_log10_adj_p <- -log10(data_7vs12_down$adjusted_p_value)

plot_7vs12_down_BP <- ggplot(data_7vs12_down, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                                  size = intersection_size, color = neg_log10_adj_p)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(x = "Number of Genes", y = "GO Biological Process", title = "GO Biological Process Enrichment 7vs12 dpf (Down)")

print(plot_7vs12_down_BP)
ggsave("7vs12_down_BP_enrichment.jpeg", plot_7vs12_down_BP, width = 12, height = 9, dpi = 300, bg = "white")

# ============================================================
# 4. 7vs12 dpf comparison (Up)
# ------------------------------------------------------------
data_7vs12_up <- read.csv("genes_7vs12_up_BP.csv", header = TRUE)
data_7vs12_up$neg_log10_adj_p <- -log10(data_7vs12_up$adjusted_p_value)

plot_7vs12_up_BP <- ggplot(data_7vs12_up, aes(x = intersection_size, y = reorder(term_name, -intersection_size), 
                                              size = intersection_size, color = neg_log10_adj_p)) +
  geom_point() +
  scale_size(name = "Number of Genes", range = c(2, 10)) +
  scale_color_gradient(low = "red", high = "green", name = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(x = "Number of Genes", y = "GO Biological Process", title = "GO Biological Process Enrichment 7vs12 dpf (Up)")

print(plot_7vs12_up_BP)
ggsave("7vs12_up_BP_enrichment.jpeg", plot_7vs12_up_BP, width = 12, height = 9, dpi = 300, bg = "white")
