library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

base_dir <- "/projectnb/selneuro/Shankar/results/snRNA-Seq/Mice/"
samples <- c("SNS2-1", "SNS2-2", "SNS2-3", "SNS2-4", "SNS3_1", "SNS3_2", "SNS3_4", "SNS3_5", "SNS4_1", "SNS4_2", "SNS4_3", "SNS4_4", "SNS4_5")
orthologs_path <- "/projectnb/selneuro/Shankar/results/snRNA-Seq/Mice/orthologs_table.csv"

# Load orthologs mapping
orthologs <- read.csv(orthologs_path)
orthologs <- orthologs[!duplicated(orthologs$mouse_gene), ] 
orthologs <- orthologs[complete.cases(orthologs), ] 
gene_mapping <- setNames(orthologs$human_gene, orthologs$mouse_gene)

seurat_objects <- list()

for (sample in samples) {
  h5_file <- file.path(base_dir, sample, "outs/filtered_feature_bc_matrix.h5")
  
  # Read the h5 file
  data <- Read10X_h5(h5_file)
  
  # Map gene names, retaining unmapped genes
  updated_gene_names <- ifelse(rownames(data) %in% names(gene_mapping),
                               gene_mapping[rownames(data)],  # Map to human orthologs
                               rownames(data))               # Retain original names for unmapped genes
  rownames(data) <- updated_gene_names
  
  # Create Seurat object
  seurat_objects[[sample]] <- CreateSeuratObject(data, project = sample)
}

# Merge all Seurat objects into a single object
combined <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[-1]
)

# Calculate mitochondrial percentage
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-") # Adjusted for human mitochondrial gene prefix

# Subset to remove cells with >5% mitochondrial genes
combined_filtered <- subset(combined, subset = percent.mt < 5)

# --------------------------
#     SCTransform step
# --------------------------
combined_filtered <- SCTransform(
  combined_filtered,
  vars.to.regress = "nFeature_RNA", 
  verbose = FALSE
)

combined_filtered <- RunPCA(combined_filtered, npcs = 50, verbose = FALSE)
combined_filtered <- RunUMAP(combined_filtered, dims = 1:50, verbose = FALSE)
combined_filtered <- FindNeighbors(combined_filtered, dims = 1:50, verbose = FALSE)
combined_filtered <- FindClusters(combined_filtered, resolution = 0.1, verbose = FALSE)

# Save object if desired
saveRDS(combined_filtered, file = "/projectnb/selneuro/Shankar/results/snRNA-Seq/ISR_snrna/rds_objects/snsm_orthologs_r01_SCT.rds")
