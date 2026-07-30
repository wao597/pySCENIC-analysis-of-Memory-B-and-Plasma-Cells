# 02a_prepare_and_azimuth.R
#
# This script:
#   1. Reads CellRanger filtered_feature_bc_matrix.h5
#   2. Applies QC filters
#   3. Normalises
#   4. Runs Azimuth with bonemarrowref
#   5. Exports Azimuth labels + confidence scores to CSV
#      
#
# Usage:
#   Rscript 02a_prepare_and_azimuth.R <donor> <cellranger_dir> <outdir>


library(Seurat)
library(Azimuth)
library(Matrix)
options(bitmapType = "cairo")

args <- commandArgs(trailingOnly = TRUE)
donor     <- args[1]   # BM1/BM2 ...
cr_dir    <- args[2]   # cellranger output directory
outdir    <- args[3]   # output directory

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sample_id <- paste0("Manton", donor)
h5_path   <- file.path(cr_dir, sample_id, "outs", "filtered_feature_bc_matrix.h5")

cat(sprintf("\n=== %s: reading %s ===\n", donor, h5_path))

if (!file.exists(h5_path)) {
  stop(sprintf("ERROR: h5 file not found: %s", h5_path))
}

#Load CellRanger filtered matrix
counts <- Read10X_h5(h5_path)
obj <- CreateSeuratObject(
  counts    = counts,
  project   = sample_id,
  min.cells = 3,
  min.features = 200
)

#QC filtering
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
cat(sprintf("%s: %d cells before QC\n", donor, ncol(obj)))

obj <- subset(
  obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           percent.mt < 10
)
cat(sprintf("%s: %d cells after QC\n", donor, ncol(obj)))

#Normalisation
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

#Run Azimuth
cat(sprintf("%s: running Azimuth...\n", donor))
obj <- RunAzimuth(obj, reference = "bonemarrowref")

#Export labels & confidence scores to CSV 
meta_cols <- c(
  "predicted.celltype.l1",
  "predicted.celltype.l1.score",
  "predicted.celltype.l2",
  "predicted.celltype.l2.score"
)
available_cols <- meta_cols[meta_cols %in% colnames(obj@meta.data)]
labels_df <- obj@meta.data[, available_cols, drop = FALSE]
labels_df$donor     <- donor
labels_df$cell_id   <- rownames(labels_df)
labels_df$n_genes   <- obj$nFeature_RNA
labels_df$pct_mt    <- obj$percent.mt

out_csv <- file.path(outdir, paste0(donor, "_azimuth_labels.csv"))
write.csv(labels_df, out_csv, row.names = FALSE)
cat(sprintf("%s: labels written to %s\n", donor, out_csv))

#Cell type summary
cat(sprintf("\n%s: Azimuth l2 cell type counts:\n", donor))
print(sort(table(obj$predicted.celltype.l2), decreasing = TRUE))

#Save full RDS
rds_out <- file.path(outdir, paste0(donor, "_azimuth_out.rds"))
saveRDS(obj, rds_out)
cat(sprintf("%s: RDS saved to %s\n", donor, rds_out))

cat(sprintf("\n=== %s: complete ===\n", donor))
