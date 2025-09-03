# ============================================================
# ATAC-seq Differential Accessibility Analysis (Ivanek protocol)
# ============================================================

# Setting the working directory
setwd("~/UoL_Bioinformatics/Research Project/Analysis/Differential_peak_analysis")

# Load libraries
library(csaw)
library(edgeR)
library(rtracklayer)
library(SummarizedExperiment)
library(GenomicRanges)

# ------------------------------------------------------------
# 1. Load BAM files
# ------------------------------------------------------------
bam_files <- c(
  "Files/1aAc_3dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam",
  "Files/1bAc_3dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam",
  "Files/2aAc_7dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam",
  "Files/2bAc_7dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam",
  "Files/3aAc_12dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam",
  "Files/3bAc_12dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam"
)

# Parameters for read counting
param <- readParam(max.frag = 1000, pe = "both")

# Count reads in 150 bp windows
data <- windowCounts(bam_files, bin = TRUE, width = 150, param = param)

# Filter out very low-count windows
keep <- filterWindowsGlobal(data)$filter > 1
filtered.data <- data[keep, ]

# Normalize counts
normfacs <- normOffsets(filtered.data)

# Define sample groups
cell_types <- c("dpf3", "dpf3", "dpf7", "dpf7", "dpf12", "dpf12")

# Convert to DGEList with coordinates
dgel <- asDGEList(filtered.data, group = cell_types, genes = as.data.frame(rowRanges(filtered.data)))
colnames(dgel) <- c("3dpf_1a", "3dpf_1b", "7dpf_2a", "7dpf_2b", "12dpf_3a", "12dpf_3b")

# ------------------------------------------------------------
# 2. Differential accessibility analysis
# ------------------------------------------------------------
moma <- model.matrix(~0 + group, data = dgel$samples)
colnames(moma) <- levels(factor(cell_types))

dgel <- estimateDisp(dgel, design = moma)
fit <- glmQLFit(dgel, moma)

# Define contrasts
contr <- makeContrasts(
  dpf3 - dpf7,
  dpf3 - dpf12,
  dpf7 - dpf12,
  levels = colnames(moma)
)

# Run tests
res_3vs7 <- glmQLFTest(fit, contrast = contr[,1])
res_3vs12 <- glmQLFTest(fit, contrast = contr[,2])
res_7vs12 <- glmQLFTest(fit, contrast = contr[,3])

# Extract results and sort by FDR
results_3vs7 <- topTags(res_3vs7, n = Inf, sort.by = "none")$table
results_3vs7 <- results_3vs7[order(results_3vs7$FDR), ]
write.table(results_3vs7, file = "differential_peaks_3vs7.tsv", sep = "\t", quote = FALSE)

results_3vs12 <- topTags(res_3vs12, n = Inf, sort.by = "none")$table
results_3vs12 <- results_3vs12[order(results_3vs12$FDR), ]
write.table(results_3vs12, file = "differential_peaks_3vs12.tsv", sep = "\t", quote = FALSE)

results_7vs12 <- topTags(res_7vs12, n = Inf, sort.by = "none")$table
results_7vs12 <- results_7vs12[order(results_7vs12$FDR), ]
write.table(results_7vs12, file = "differential_peaks_7vs12.tsv", sep = "\t", quote = FALSE)

# ------------------------------------------------------------
# 3. Filter for significance 
# ------------------------------------------------------------
results_3vs7_sig <- subset(results_3vs7, FDR < 0.05)
results_3vs12_sig <- subset(results_3vs12, FDR < 0.05)
results_7vs12_sig <- subset(results_7vs12, FDR < 0.05)

# ------------------------------------------------------------
# 4. Load fixed promoter BED file
# ------------------------------------------------------------
promoter_bed <- read.table("Files/Astatotilapia_calliptera.fAstCal1.2.101.5kb_promoters.bed", 
                           sep="\t", colClasses=c("character", "integer", "integer", 
                                                  "character", "character", "character"))
bed_fixed <- data.frame(
  chr = promoter_bed[,1],
  start = promoter_bed[,2],
  end = promoter_bed[,3],
  name = promoter_bed[,4],
  score = 0,
  strand = promoter_bed[,6]
)
write.table(bed_fixed, "Files/Astatotilapia_calliptera.fAstCal1.2.101.5kb_promoters_fixed.bed", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

promoters <- import.bed("Files/Astatotilapia_calliptera.fAstCal1.2.101.5kb_promoters_fixed.bed")

# ------------------------------------------------------------
# 5. Convert significant results to GRanges
# ------------------------------------------------------------
results_gr_3vs7 <- makeGRangesFromDataFrame(results_3vs7_sig, 
                                            seqnames.field = "seqnames", 
                                            start.field = "start", end.field = "end", 
                                            keep.extra.columns = TRUE)
results_gr_3vs12 <- makeGRangesFromDataFrame(results_3vs12_sig, 
                                             seqnames.field = "seqnames", 
                                             start.field = "start", end.field = "end", 
                                             keep.extra.columns = TRUE)
results_gr_7vs12 <- makeGRangesFromDataFrame(results_7vs12_sig, 
                                             seqnames.field = "seqnames", 
                                             start.field = "start", end.field = "end", 
                                             keep.extra.columns = TRUE)

# ------------------------------------------------------------
# 6. Find overlaps with promoters
# ------------------------------------------------------------
overlaps_3vs7 <- findOverlaps(results_gr_3vs7, promoters)
overlaps_3vs12 <- findOverlaps(results_gr_3vs12, promoters)
overlaps_7vs12 <- findOverlaps(results_gr_7vs12, promoters)

# ------------------------------------------------------------
# 7. Save overlapping promoter peaks and gene IDs
# ------------------------------------------------------------
promoter_peaks_3vs7 <- results_3vs7_sig[queryHits(overlaps_3vs7), ]
promoter_peaks_3vs7$gene_id <- mcols(promoters)$name[subjectHits(overlaps_3vs7)]

promoter_peaks_3vs12 <- results_3vs12_sig[queryHits(overlaps_3vs12), ]
promoter_peaks_3vs12$gene_id <- mcols(promoters)$name[subjectHits(overlaps_3vs12)]

promoter_peaks_7vs12 <- results_7vs12_sig[queryHits(overlaps_7vs12), ]
promoter_peaks_7vs12$gene_id <- mcols(promoters)$name[subjectHits(overlaps_7vs12)]

# Save tables
write.table(promoter_peaks_3vs7, file = "promoter_differential_peaks_3vs7.tsv", sep = "\t", quote = FALSE)
write.table(promoter_peaks_3vs12, file = "promoter_differential_peaks_3vs12.tsv", sep = "\t", quote = FALSE)
write.table(promoter_peaks_7vs12, file = "promoter_differential_peaks_7vs12.tsv", sep = "\t", quote = FALSE)

# Save unique gene lists (combined significant promoter peaks)
write.table(unique(promoter_peaks_3vs7$gene_id), file = "genes_3vs7.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique(promoter_peaks_3vs12$gene_id), file = "genes_3vs12.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique(promoter_peaks_7vs12$gene_id), file = "genes_7vs12.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ------------------------------------------------------------
# 7b. Split promoter peaks into activated (logFC > 0) and 
#     inactivated (logFC < 0) and save lists for GO
# ------------------------------------------------------------
# NOTE about direction: contrasts are dpf3 - dpf7, dpf3 - dpf12, dpf7 - dpf12.
# So logFC > 0 in 'dpf3 - dpf7' means higher at 3 dpf (relative to 7 dpf).

# 3 vs 7
promoter_3vs7_up   <- subset(promoter_peaks_3vs7, logFC > 0)
promoter_3vs7_down <- subset(promoter_peaks_3vs7, logFC < 0)
write.table(promoter_3vs7_up,   file = "promoter_differential_peaks_3vs7_up.tsv",   sep = "\t", quote = FALSE)
write.table(promoter_3vs7_down, file = "promoter_differential_peaks_3vs7_down.tsv", sep = "\t", quote = FALSE)
write.table(unique(promoter_3vs7_up$gene_id),   file = "genes_3vs7_up.txt",   row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique(promoter_3vs7_down$gene_id), file = "genes_3vs7_down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 3 vs 12
promoter_3vs12_up   <- subset(promoter_peaks_3vs12, logFC > 0)
promoter_3vs12_down <- subset(promoter_peaks_3vs12, logFC < 0)
write.table(promoter_3vs12_up,   file = "promoter_differential_peaks_3vs12_up.tsv",   sep = "\t", quote = FALSE)
write.table(promoter_3vs12_down, file = "promoter_differential_peaks_3vs12_down.tsv", sep = "\t", quote = FALSE)
write.table(unique(promoter_3vs12_up$gene_id),   file = "genes_3vs12_up.txt",   row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique(promoter_3vs12_down$gene_id), file = "genes_3vs12_down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 7 vs 12
promoter_7vs12_up   <- subset(promoter_peaks_7vs12, logFC > 0)
promoter_7vs12_down <- subset(promoter_peaks_7vs12, logFC < 0)
write.table(promoter_7vs12_up,   file = "promoter_differential_peaks_7vs12_up.tsv",   sep = "\t", quote = FALSE)
write.table(promoter_7vs12_down, file = "promoter_differential_peaks_7vs12_down.tsv", sep = "\t", quote = FALSE)
write.table(unique(promoter_7vs12_up$gene_id),   file = "genes_7vs12_up.txt",   row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique(promoter_7vs12_down$gene_id), file = "genes_7vs12_down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ------------------------------------------------------------
# 8. Quick summary 
# ------------------------------------------------------------
print("Significant peaks (FDR < 0.05):")
print(nrow(results_3vs7_sig))
print(nrow(results_3vs12_sig))
print(nrow(results_7vs12_sig))

print("Promoter-overlapping significant peaks (counts):")
print(paste0("3vs7 total = ", nrow(promoter_peaks_3vs7), 
             "; up = ", nrow(promoter_3vs7_up), "; down = ", nrow(promoter_3vs7_down)))
print(paste0("3vs12 total = ", nrow(promoter_peaks_3vs12), 
             "; up = ", nrow(promoter_3vs12_up), "; down = ", nrow(promoter_3vs12_down)))
print(paste0("7vs12 total = ", nrow(promoter_peaks_7vs12), 
             "; up = ", nrow(promoter_7vs12_up), "; down = ", nrow(promoter_7vs12_down)))

print("Unique overlapping peaks with promoters:")
print(length(unique(queryHits(overlaps_3vs7))))
print(length(unique(queryHits(overlaps_3vs12))))
print(length(unique(queryHits(overlaps_7vs12))))

print("Unique genes in promoter peaks (combined):")
print(length(unique(promoter_peaks_3vs7$gene_id)))
print(length(unique(promoter_peaks_3vs12$gene_id)))
print(length(unique(promoter_peaks_7vs12$gene_id)))

print("Unique genes in promoter peaks (directional):")
print(paste0("3vs7 up genes = ", length(unique(promoter_3vs7_up$gene_id)), 
             "; 3vs7 down genes = ", length(unique(promoter_3vs7_down$gene_id))))
print(paste0("3vs12 up genes = ", length(unique(promoter_3vs12_up$gene_id)), 
             "; 3vs12 down genes = ", length(unique(promoter_3vs12_down$gene_id))))
print(paste0("7vs12 up genes = ", length(unique(promoter_7vs12_up$gene_id)), 
             "; 7vs12 down genes = ", length(unique(promoter_7vs12_down$gene_id))))

# ------------------------------------------------------------
# 9. QC: MDS plot (promoter-focused)
# ------------------------------------------------------------
promoter_idx <- unique(c(queryHits(overlaps_3vs7),
                         queryHits(overlaps_3vs12),
                         queryHits(overlaps_7vs12)))
dgel_promoter <- dgel[promoter_idx, , keep.lib.sizes=FALSE]
dgel_promoter <- calcNormFactors(dgel_promoter)

png("mds_plot_promoters.png", width = 600, height = 400)
plotMDS(dgel_promoter, labels = colnames(dgel_promoter),
        col = c("blue", "blue", "green", "green", "red", "red"))
dev.off()

