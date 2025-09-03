# Load required package
library(rtracklayer)

# -----------------------------
# Step 1: Import GTF and extract gene_id → gene_name mapping
gtf_file <- "Astatotilapia_calliptera.fAstCal1.2.101.gtf"
gtf <- import(gtf_file)

gene_map <- unique(mcols(gtf)[, c("gene_id", "gene_name")])
gene_map <- as.data.frame(gene_map)

# -----------------------------
# Step 2: Function to map CSV IDs to gene names
# Handles both multi-line and single-line comma-separated files
map_ids_csv <- function(file, gene_map) {
  # Read file
  lines <- readLines(file)
  
  # Split by comma and remove extra spaces
  ids <- trimws(unlist(strsplit(lines, ",")))
  
  # Map to gene names
  mapped <- merge(data.frame(gene_id = ids),
                  gene_map,
                  by.x = "gene_id", by.y = "gene_id",
                  all.x = TRUE)
  return(mapped)
}

# -----------------------------
# Step 3: Map the three files
aa_met_process                 <- map_ids_csv("aa_met_process.txt", gene_map)
actin_cyto_org                 <- map_ids_csv("actin_cyto_org.txt", gene_map)
carbohydrate_meta_pro          <- map_ids_csv("carbohydrate_meta_pro.txt", gene_map)
carbohydrate_deri_cat_pro      <- map_ids_csv("carbohydrate_deri_cat_pro.txt", gene_map)
cell_migration                 <- map_ids_csv("cell_migration.txt", gene_map)
cyto_org                       <- map_ids_csv("cyto_org.txt", gene_map)
digestive_tract_dev            <- map_ids_csv("Digestive_tract_dev.txt", gene_map)
growth_inv_in_morpho           <- map_ids_csv("growth_inv_in_morpho.txt", gene_map)
loco_behavior                  <- map_ids_csv("loco_behavior.txt", gene_map)

# -----------------------------
# Step 4: Save mapped results
write.table(aa_met_process, "aa_met_process.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(actin_cyto_org, "actin_cyto_org.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(carbohydrate_meta_pro, "carbohydrate_meta_pro.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(carbohydrate_deri_cat_pro, "carbohydrate_deri_cat_pro.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(cell_migration, "cell_migration.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(cyto_org, "cyto_org.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(digestive_tract_dev, "digestive_tract_dev.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(growth_inv_in_morpho, "growth_inv_in_morpho.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(loco_behavior, "loco_behavior.tsv", sep="\t", row.names=FALSE, quote=FALSE)
