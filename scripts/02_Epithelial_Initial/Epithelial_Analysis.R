# Stage 02 — Epi subset → SCT/Harmony → clustering (res 0.4)
# epi_clustered.rds 가 있으면 compute 는 건너뛰고 figures 만 재생성.

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

IN_RDS  <- "Results/01_Integrated/combined_CRPC.rds"
OUT_DIR <- "Results/02_Epithelial_Initial"
EPI_RDS <- file.path(OUT_DIR, "epi_clustered.rds")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

RES_VEC  <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2)  # clustree · UMAP_multires 용
RES_BASE <- "0.4"                                      # seurat_clusters 로 쓰는 해상도

# ============================================================
# Compute: Epi subset → SCT/Harmony → clustering → UMAP → cell cycle ----
# ============================================================
build_epi_object <- function() {
    epi <- subset(readRDS(IN_RDS), subset = celltype == "Epithelial")
    message("Epithelial subset: ", ncol(epi), " cells")

    # Reset assay state (drop inherited SCT, re-split RNA by patient for Harmony)
    DefaultAssay(epi) <- "RNA"
    epi[["SCT"]] <- NULL
    epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)
    epi@meta.data <- dplyr::select(epi@meta.data,
        -matches("_snn_res\\."), -seurat_clusters)

    epi <- utils_sct_harmony(epi)
    utils_save_png(ElbowPlot(epi, ndims = 50), OUT_DIR, "ElbowPlot.png", w = 12, h = 8)

    epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:30)
    epi <- FindClusters(epi, resolution = RES_VEC)
    # FindClusters 는 마지막 res 를 seurat_clusters 로 둠 → RES_BASE 로
    epi$seurat_clusters <- epi@meta.data[[paste0("SCT_snn_res.", RES_BASE)]]
    Idents(epi) <- "seurat_clusters"

    epi <- RunUMAP(epi, reduction = "harmony", dims = 1:30,
                   n.neighbors = 20, min.dist = 0.1, spread = 4.0)
    epi <- CellCycleScoring(epi, s.features = cc.genes.updated.2019$s.genes,
                            g2m.features = cc.genes.updated.2019$g2m.genes, nbin = 12)

    utils_save_all_markers(epi, file.path(OUT_DIR, "all_markers.csv"))
    saveRDS(epi, EPI_RDS)
    epi
}

if (file.exists(EPI_RDS)) {
    message("Loading ", EPI_RDS, " — skipping compute")
    epi <- readRDS(EPI_RDS)
} else {
    epi <- build_epi_object()
}

# ============================================================
# Figures ----
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
    ggtitle("Epithelial clusters")
utils_save_png(p, OUT_DIR, "UMAP_cluster.png", w = 10, h = 8)
p <- utils_umap_by(epi, "orig.ident", pt.size = 0.3) + ggtitle("By patient")
utils_save_png(p, OUT_DIR, "UMAP_patient.png", w = 10, h = 8)

# Patient-split UMAP
p <- DimPlot(epi, split.by = "orig.ident", label = TRUE, pt.size = 0.3, cols = cluster_cols)
utils_save_png(p, OUT_DIR, "UMAP_split_by_patient.png", w = 24, h = 8)

p <- utils_umap_by(epi, "Phase", pt.size = 0.3)
utils_save_png(p, OUT_DIR, "CellCycle_UMAP.png", w = 10, h = 7)

# Marker 패널 ----
LUMINAL_BASAL_CLUB <- c(
    "FOXA1", "HOXB13", "NKX3-1", "KLK3", "TOP2A", "MKI67",  # Luminal
    "KRT5", "KRT14", "TP63",                                # Basal
    "MMP7", "WFDC2"                                         # Club
)
DNPC_MARKERS <- c("CHD7", "MYC", "KMT2C", "KRT7", "SOX2", "SYP", "AR")

p <- FeaturePlot(epi, features = LUMINAL_BASAL_CLUB, ncol = 4, pt.size = 0.1)
utils_save_png(p, OUT_DIR, "Luminal_Basal_Club_FeaturePlot.png", w = 20, h = 16)
p <- DotPlot(epi, features = LUMINAL_BASAL_CLUB) + RotatedAxis()
utils_save_png(p, OUT_DIR, "Luminal_Basal_Club_DotPlot.png", w = 14, h = 8)
p <- FeaturePlot(epi, features = DNPC_MARKERS, ncol = 4, pt.size = 0.1)
utils_save_png(p, OUT_DIR, "DNPC_FeaturePlot.png", w = 20, h = 16)
p <- DotPlot(epi, features = DNPC_MARKERS) + RotatedAxis()
utils_save_png(p, OUT_DIR, "DNPC_DotPlot.png", w = 12, h = 8)

marker_dotplot <- function(features, title) {
    DotPlot(epi, features = features, cluster.idents = FALSE) +
        RotatedAxis() +
        scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
        labs(title = title)
}

# CRPC subtype panel: ARPC / DNPC-FGF / DNPC-KRT7 / NEPC
SUBTYPE_MARKERS <- list(
    ARPC        = c("AR", "KLK3", "FOLH1"),
    `DNPC-FGF`  = c("FGF8", "FGFR1", "CCL2", "DKK1"),
    `DNPC-KRT7` = c("KRT7", "SOX2", "FOXA2"),
    NEPC        = c("CHGA", "SYP", "PROX1")
)
p <- marker_dotplot(SUBTYPE_MARKERS,
                    "CRPC subtype markers (ARPC / DNPC-FGF / DNPC-KRT7 / NEPC)")
utils_save_png(p, OUT_DIR, "Subtype_DotPlot.png", w = 14, h = 8)

# Hillock / club-hillock verification panel
HILLOCK_MARKERS <- list(
    Hillock_core  = c("KRT13", "KRT4", "AQP3", "LY6D", "PSCA", "LYPD3"),
    Basal_axis    = c("KRT5", "KRT14", "TP63"),
    Squamous_diff = c("SERPINB5", "SERPINB1", "ZNF750", "S100A2"),
    Club          = c("MMP7", "WFDC2", "PIGR", "LCN2"),
    Luminal_ref   = c("KRT8", "KRT18", "AR", "KLK3")
)
p <- marker_dotplot(HILLOCK_MARKERS, "Hillock / club-hillock verification")
utils_save_png(p, OUT_DIR, "Hillock_DotPlot.png", w = 16, h = 8)

message("Stage 02 results -> ", OUT_DIR)
