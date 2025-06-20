library(Seurat)
library(sctransform)
library(glmGamPoi)
library(dplyr)
library(ggplot2)
library(patchwork)

#-------------------------------
# 1) Load & Prepare Data
#-------------------------------
mouse_exc <- readRDS("/projectnb/selneuro/Shankar/results/snRNA-Seq/Mice/mouse_orthologs_exc_cells_SCT.rds")
integrated_regressed <- readRDS("/projectnb/selneuro/Shankar/results/snRNA-Seq/SCT/integrated_regressed_SCT_exc_subset.rds")

DefaultAssay(mouse_exc) <- "RNA"
DefaultAssay(integrated_regressed) <- "RNA"

# Remove old SCT assay if it exists
mouse_exc[["SCT"]] <- NULL
integrated_regressed[["SCT"]] <- NULL
integrated_regressed[["integrated"]] <- NULL

#-------------------------------
# 2) Subset to common gene names
#   (IMPORTANT for cross-species)
#-------------------------------
common_genes <- intersect(rownames(mouse_exc), rownames(integrated_regressed))

mouse_exc <- subset(mouse_exc, features = common_genes)
integrated_regressed <- subset(integrated_regressed, features = common_genes)

#-------------------------------
# 3) Tag each object with sampleID
#-------------------------------
mouse_exc$sampleID <- mouse_exc$orig.ident
mouse_exc$orig.ident <- NULL

integrated_regressed$sampleID <- integrated_regressed$orig.ident
integrated_regressed$orig.ident <- NULL

mouse_exc$dataset <- "mouse_exc"
integrated_regressed$dataset <- "human_exc"

#-------------------------------
# 4) Merge everything into one
#-------------------------------
combined_obj <- merge(
  x = mouse_exc,
  y = integrated_regressed,
  project = "mouse_rosmap_merged"
)

#-------------------------------
# 5) Split by sample ID -> list
#-------------------------------
combined_list <- SplitObject(combined_obj, split.by = "sampleID")

#-------------------------------
# 6) Ensure All Objects Have the Same Genes
#-------------------------------
common_genes <- Reduce(intersect, lapply(combined_list, rownames))

combined_list <- lapply(combined_list, function(obj) {
  subset(obj, features = common_genes)
})

# Verify uniform gene count across all objects
print(lapply(combined_list, function(obj) length(rownames(obj))))

#-------------------------------
# 7) SCTransform each sub-object
#-------------------------------
combined_list <- lapply(
  combined_list,
  function(obj) {
    SCTransform(obj, vars.to.regress = "nFeature_RNA")
  }
)

#-------------------------------
# 8) Integration Workflow
#-------------------------------
features <- SelectIntegrationFeatures(combined_list)
combined_list <- PrepSCTIntegration(combined_list, anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list        = combined_list, 
  normalization.method = "SCT",
  anchor.features    = features
)

combined_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

#-------------------------------
# 9) Downstream Analysis
#-------------------------------
DefaultAssay(combined_integrated) <- "integrated"

combined_integrated <- RunPCA(combined_integrated, npcs = 50)
combined_integrated <- RunUMAP(combined_integrated, dims = 1:50)
combined_integrated <- FindNeighbors(combined_integrated, dims = 1:50)
combined_integrated <- FindClusters(combined_integrated, resolution = 1)

#-------------------------------
# 10) Save Output
#-------------------------------
saveRDS(combined_integrated, file = "/projectnb/selneuro/Shankar/results/snRNA-Seq/Mice/mouse_human_SCT.rds")

