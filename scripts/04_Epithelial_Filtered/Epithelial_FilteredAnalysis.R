# Stage 04 — filter_decision.csv로 저품질 cluster drop → SCT/Harmony/clustering 재계산
# epi_filtered_clustered.rds 가 있으면 compute 는 건너뛰고 figures만 재생성
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(harmony)
library(clustree)
library(future)

source("scripts/00_utils/scRNA_utils.R")

# Settings ----
N_CORES <- 30
plan("multicore", workers = N_CORES)
options(future.globals.maxSize = 128 * 1024^3)

IN_RDS       <- "Results/02_Epithelial_Initial/epi_clustered.rds"
DECISION_CSV <- "Results/03_Epithelial_FilterDecision/filter_decision.csv"
OUT_DIR      <- "Results/04_Epithelial_Filtered"
OUT_RDS      <- file.path(OUT_DIR, "epi_filtered_clustered.rds")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 2026-07-11: decontX 적용 후에는 모든 cluster가 적합하여 filter_decision.csv에서 모두 keep=TRUE → drop 없음.
RES_VEC  <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2)  # for clustree, UMAP_multires
RES_BASE <- "0.3"                                      # seurat_clusters resolution

# ============================================================
# Compute: cluster 필터 → SCT/Harmony → clustering → UMAP → cell cycle ----
# ============================================================
build_filtered_object <- function() {
    epi <- readRDS(IN_RDS)

    decision <- read.csv(DECISION_CSV)
    keep_clusters <- as.character(
        decision$cluster[decision$keep %in% c(TRUE, "TRUE", "true")]
    )
    drop_clusters <- setdiff(as.character(decision$cluster), keep_clusters)
    message("Keep ", length(keep_clusters), " clusters; drop ", length(drop_clusters))
    if (length(drop_clusters) > 0) {
        message("  dropped: ", paste(drop_clusters, collapse = ", "))
    }

    n_before <- ncol(epi)
    keep_cell <- as.character(epi$seurat_clusters) %in% keep_clusters
    epi <- subset(epi, cells = colnames(epi)[keep_cell])
    message("Cells: ", n_before, " → ", ncol(epi),
            " (", n_before - ncol(epi), " dropped)")

    # Reset assay/reduction state (drop stage-02 SCT/pca/umap, re-split RNA by patient)
    DefaultAssay(epi) <- "RNA"
    epi[["SCT"]]     <- NULL
    epi[["harmony"]] <- NULL
    epi[["pca"]]     <- NULL
    epi[["umap"]]    <- NULL
    epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)
    epi@meta.data <- dplyr::select(epi@meta.data,
        -matches("_snn_res\\."), -seurat_clusters)

    epi <- utils_sct_harmony(epi)
    utils_save_png(ElbowPlot(epi, ndims = 50), OUT_DIR, "ElbowPlot.png", w = 12, h = 8)

    epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:30)
    epi <- FindClusters(epi, resolution = RES_VEC)
    # FindClusters는 마지막 res가 seurat_clusters로 들어감 → RES_BASE(0.3)로 재설정
    epi$seurat_clusters <- epi@meta.data[[paste0("SCT_snn_res.", RES_BASE)]]
    Idents(epi) <- "seurat_clusters"

    epi <- RunUMAP(epi, reduction = "harmony", dims = 1:30,
                   n.neighbors = 40, min.dist = 0.3, spread = 1.5)
    # Cell cycle re-scoring (SCT 모델이 바뀌었으므로 재계산)
    epi <- CellCycleScoring(epi, s.features = cc.genes.updated.2019$s.genes,
                            g2m.features = cc.genes.updated.2019$g2m.genes, nbin = 12)

    saveRDS(epi, OUT_RDS)
    epi
}

if (file.exists(OUT_RDS)) {
    message("Loading ", OUT_RDS, " — skipping compute")
    epi <- readRDS(OUT_RDS)
} else {
    epi <- build_filtered_object()
}

# ============================================================
# Diagnostic figures ----
# ============================================================
# clustree: cluster stability tree across resolutions
p <- clustree(epi, prefix = "SCT_snn_res.")
utils_save_png(p, OUT_DIR, "clustree.png", w = 12, h = 14, dpi = 200)

# UMAPs at each resolution for side-by-side comparison
multires <- lapply(RES_VEC, function(r) {
    utils_umap_by(epi, paste0("SCT_snn_res.", r), label = TRUE, pt.size = 0.2) +
        ggtitle(paste0("res = ", r)) + NoLegend()
})
utils_save_png(wrap_plots(multires, ncol = 3), OUT_DIR, "UMAP_multires.png",
               w = 18, h = 14, dpi = 200)

cluster_cols <- utils_cb_cols(epi$seurat_clusters)
p <- DimPlot(epi, label = TRUE, pt.size = 0.3, cols = cluster_cols) +
    ggtitle(paste0("Filtered epithelial clusters (decontX, res ", RES_BASE, ")"))
utils_save_png(p, OUT_DIR, "UMAP_cluster.png", w = 10, h = 8)
p <- utils_umap_by(epi, "orig.ident", pt.size = 0.3) + ggtitle("By patient")
utils_save_png(p, OUT_DIR, "UMAP_patient.png", w = 10, h = 8)

p <- DimPlot(epi, split.by = "orig.ident", label = TRUE, pt.size = 0.3, cols = cluster_cols)
utils_save_png(p, OUT_DIR, "UMAP_split_by_patient.png", w = 24, h = 8)

p <- utils_umap_by(epi, "Phase", pt.size = 0.3)
utils_save_png(p, OUT_DIR, "CellCycle_UMAP.png", w = 10, h = 7)

# Post-filter QC verification ----
qc_vlns <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(feat) {
    VlnPlot(epi, features = feat, pt.size = 0, cols = cluster_cols) +
        NoLegend() + labs(x = "Cluster", title = feat)
})
utils_save_png(wrap_plots(qc_vlns, ncol = 1), OUT_DIR, "cluster_QC_violin.png",
               w = 12, h = 12)

lineage_markers <- list(
    Epithelial    = c("EPCAM", "KRT8", "KRT18", "KRT19"),
    Basal_Hillock = c("KRT5", "KRT14", "TP63", "KRT13", "KRT4"),
    Fibroblast    = c("COL1A1", "COL1A2", "PDGFRA"),
    Smooth_muscle = c("ACTA2", "MYH11", "RGS5", "TAGLN"),
    Endothelial   = c("PECAM1", "VWF"),
    Immune        = c("PTPRC", "CD68")
)
p <- DotPlot(epi, features = lineage_markers, cluster.idents = FALSE) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    labs(title = "Post-filter lineage purity check")
utils_save_png(p, OUT_DIR, "cluster_lineage_DotPlot.png", w = 18, h = 8)

message("Stage 04 results -> ", OUT_DIR, " | canonical res ", RES_BASE,
        " | ", nlevels(epi$seurat_clusters), " clusters, ", ncol(epi), " cells")
