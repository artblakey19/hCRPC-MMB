# Stage 07-02 epithelial subclustering → 세부 annotation
# 07-01 의 celltype == "Epithelial" 만 사용
# <PATIENT>_epi.rds 가 있으면 compute 는 건너뛰고 figures / annotation만 생성

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(clustree)
library(future)

source("scripts/00_utils/scRNA_utils.R")

# Settings ----
PATIENT <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(PATIENT %in% c("CRPC1", "CRPC2", "CRPC3"))

N_CORES <- 30
plan("multicore", workers = N_CORES)
options(future.globals.maxSize = 128 * 1024^3)

IN_RDS  <- file.path("Results/07_PerPatient", PATIENT, "01_Major", paste0(PATIENT, "_major.rds"))
INT_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- file.path("Results/07_PerPatient", PATIENT, "02_Epithelial")
OBJ_RDS <- file.path(OUT_DIR, paste0(PATIENT, "_epi.rds"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

RES_VEC <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)  # clustree, UMAP_multires 용
RES_EPI <- "0.4"                                      # seurat_clusters 해상도

# ============================================================
# Compute: Epi subset → SCT → clustering → UMAP → cell cycle ----
# ============================================================
# integrated anlysis에서 상피로 annotation을 불러와 메타데이터 추가
# integrated anlysis에서 상피로 분류되지 않았으나 환자별 분석에서 상피로 분류되었다면 NA 로 처리
integrated_annotation <- function(cells) {
    int_epi <- readRDS(INT_RDS)
    as.character(int_epi$annotation)[match(cells, colnames(int_epi))]
}

build_epi_object <- function() {
    epi <- subset(readRDS(IN_RDS), subset = celltype == "Epithelial")
    message(PATIENT, " epithelial: ", ncol(epi), " cells")

    epi$major_cluster <- as.character(epi$seurat_clusters)
    epi$integrated_epi_annotation <- integrated_annotation(colnames(epi))

    DefaultAssay(epi) <- "RNA"
    epi[["SCT"]] <- NULL
    epi@meta.data <- dplyr::select(epi@meta.data,
        -matches("_snn_res\\."), -seurat_clusters, -celltype)

    epi <- utils_sct_pca(epi)
    utils_save_png(ElbowPlot(epi, ndims = 50), OUT_DIR, "ElbowPlot.png", w = 10, h = 6)

    epi <- FindNeighbors(epi, reduction = "pca", dims = 1:30)
    epi <- FindClusters(epi, resolution = RES_VEC, random.seed = 42)
    epi$seurat_clusters <- epi@meta.data[[paste0("SCT_snn_res.", RES_EPI)]]
    Idents(epi) <- "seurat_clusters"
    epi <- RunUMAP(epi, reduction = "pca", dims = 1:30,
                   n.neighbors = 20, min.dist = 0.1, spread = 4.0)
    epi <- CellCycleScoring(epi, s.features = cc.genes.updated.2019$s.genes,
                            g2m.features = cc.genes.updated.2019$g2m.genes, nbin = 12)

    utils_save_all_markers(epi, file.path(OUT_DIR, "all_markers.csv"))
    saveRDS(epi, OBJ_RDS)
    epi
}

if (file.exists(OBJ_RDS)) {
    message("Loading cached ", OBJ_RDS, " — skipping compute")
    epi <- readRDS(OBJ_RDS)
} else {
    epi <- build_epi_object()
}

# ============================================================
# Figures ----
# ============================================================
p <- clustree(epi, prefix = "SCT_snn_res.")
utils_save_png(p, OUT_DIR, "clustree.png", w = 12, h = 14, dpi = 200)

multires <- lapply(RES_VEC, function(r) {
    utils_umap_by(epi, paste0("SCT_snn_res.", r), label = TRUE, pt.size = 0.2) +
        ggtitle(paste0("res = ", r)) + NoLegend()
})
utils_save_png(wrap_plots(multires, ncol = 4), OUT_DIR, "UMAP_multires.png",
               w = 24, h = 12, dpi = 200)

p <- utils_umap_by(epi, "seurat_clusters", label = TRUE, pt.size = 0.3) +
    ggtitle(paste0(PATIENT, " alone — epithelial clusters (res ", RES_EPI, ")"))
utils_save_png(p, OUT_DIR, "UMAP_cluster.png", w = 10, h = 8)

p <- utils_umap_by(epi, "integrated_epi_annotation", label = TRUE, repel = TRUE, pt.size = 0.3) +
    ggtitle("Integrated epithelial annotation on patient-alone UMAP")
utils_save_png(p, OUT_DIR, "UMAP_integrated_annotation.png", w = 11, h = 8)

p <- utils_umap_by(epi, "Phase", pt.size = 0.3)
utils_save_png(p, OUT_DIR, "UMAP_CellCycle.png", w = 10, h = 7)

p <- VlnPlot(epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination"),
             group.by = "seurat_clusters", pt.size = 0, ncol = 2,
             cols = utils_cb_cols(epi$seurat_clusters))
utils_save_png(p, OUT_DIR, "QC_violin_by_cluster.png", w = 16, h = 10)

# Epithelial subtype marker ----
EPI_MARKERS <- list(
    Luminal_AR  = c("AR", "KLK3", "KLK2", "NKX3-1", "FOLH1", "ACP3", "TMPRSS2", "HOXB13", "FOXA1"),
    Luminal_pan = c("KRT8", "KRT18", "CD24"),
    Basal       = c("KRT5", "KRT14", "KRT15", "TP63", "COL17A1"),
    Club        = c("MMP7", "WFDC2", "PIGR", "SCGB3A1", "SCGB1A1", "LCN2", "LTF"),
    Hillock     = c("KRT13", "KRT4", "AQP3", "LY6D", "PSCA", "S100A2", "KRT6A"),
    OE_CFTR     = c("CFTR", "SLC4A4", "MUC5B", "OLFM4"),
    Ionocyte    = c("FOXI1", "ATP6V1B1", "ASCL3"),
    NE          = c("CHGA", "CHGB", "SYP", "ENO2", "ASCL1", "INSM1", "NCAM1"),
    DNPC        = c("KRT7", "SOX2", "FOXA2", "FGFR1", "FGF8", "DKK1", "CHD7", "MYC"),
    Tumor_nonAR = c("PCA3", "SCHLAP1", "GOLM1", "AMACR", "ERG"),
    Stress_AP1  = c("FOS", "JUN", "JUNB", "ATF3", "HSPA1A", "DDIT3"),
    Cycling     = c("MKI67", "TOP2A"),
    Squamous    = c("KRT1", "KRT10", "SPRR1B", "ZNF750"),
    Urothelial  = c("UPK1B", "UPK3A"),
    Contam      = c("PTPRC", "DCN", "PECAM1")
)

EPI_MARKERS <- lapply(EPI_MARKERS, intersect, rownames(epi))

marker_dotplot <- function(group, title = NULL) {
    DotPlot(epi, features = EPI_MARKERS, group.by = group, cluster.idents = FALSE) +
        RotatedAxis() +
        scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
        theme(axis.text.x = element_text(size = 7)) +
        labs(title = title)
}
p <- marker_dotplot("seurat_clusters", paste0(PATIENT, " alone — epithelial subtype markers by cluster"))
utils_save_png(p, OUT_DIR, "DotPlot_epi_markers_by_cluster.png", w = 30, h = 8)

FP_GENES <- c("AR", "KLK3", "NKX3-1", "KRT8", "KRT5", "KRT14", "TP63", "MMP7", "PIGR",
              "KRT13", "KRT4", "CFTR", "FOXI1", "SYP", "CHGA", "KRT7", "SOX2", "PCA3",
              "AMACR", "FOS", "MKI67", "TOP2A", "S100A2", "UPK1B")
p <- FeaturePlot(epi, features = FP_GENES, ncol = 6, pt.size = 0.1)
utils_save_png(p, OUT_DIR, "FeaturePlot_epi_markers.png", w = 30, h = 20)

# 다른 annotation과 비교----
utils_write_crosstab(epi$seurat_clusters, addNA(epi$integrated_epi_annotation),
                     OUT_DIR, "cluster_vs_integrated_annotation.csv")
utils_write_crosstab(epi$seurat_clusters, epi$major_cluster, OUT_DIR, "cluster_vs_major_cluster.csv")

marker_genes <- unlist(EPI_MARKERS)
avg <- AverageExpression(epi, features = marker_genes, group.by = "seurat_clusters",
                         assay = "SCT", layer = "data")$SCT
write.csv(round(as.matrix(avg), 3), file.path(OUT_DIR, "epi_marker_avg_by_cluster.csv"))

expressed <- GetAssayData(epi, assay = "SCT", layer = "data")[marker_genes, ] > 0
pct <- sapply(split(colnames(epi), epi$seurat_clusters),
              function(cells) Matrix::rowMeans(expressed[, cells]))
write.csv(round(pct, 3), file.path(OUT_DIR, "epi_marker_pct_by_cluster.csv"))

qc_tab <- epi@meta.data |>
    group_by(cluster = seurat_clusters) |>
    summarise(n = n(),
              median_nCount = median(nCount_RNA),
              median_nFeature = median(nFeature_RNA),
              median_mt = round(median(percent.mt), 2),
              median_decontX = round(median(decontX_contamination), 3),
              frac_cycling = round(mean(Phase != "G1"), 3),
              .groups = "drop")
write.csv(qc_tab, file.path(OUT_DIR, "cluster_QC_summary.csv"), row.names = FALSE)

# ============================================================
# Annotation (res 0.4) ----
# ============================================================
ANNOTATION_MAPS <- list(
    CRPC1 = c(
        `0`  = "Club MMP7/SCGB3A1",              # WFDC2/OLFM4/FOLR1, integrated Club
        `1`  = "OE CFTR (low-count)",            # SLC4A4/PDE8B, integrated OE
        `2`  = "Hillock KRT13/PSCA",             # LY6D/UPK1B/S100P, integrated Hillock 1
        `3`  = "Basal TP63 (low-count)",         # integrated BE 6
        `4`  = "Basal KRT5/KRT15 (AP-1 high)",   # COL17A1/CRYAB/SFRP1, integrated BE 1 + BE 5
        `5`  = "Basal-squamous DSG3/KRT13",      # DSC3/OLFM4/SNAI2, integrated BE 1 + Hillock 2
        `6`  = "Basal (low-count)",              # 핵 lncRNA, mt ~10%, integrated BE 4
        `7`  = "Basal DKK1/SHH",                 # ALDH3A1/S100A2/AKR1C3, integrated BE 2
        `8`  = "Luminal TMPRSS2 (low-count)",    # SHROOM3/PDE4C, integrated OE + BE 4
        `9`  = "Low-quality PPARG (mt high)",    # mt ~12%, SCHLAP1, integrated BE 4
        `10` = "Club-Hillock OLFM4/MMP7",        # GNG4/CALML3/GABRP/KRT13, integrated BE 5 + Club
        `11` = "Luminal AR+ (KLK3)"              # KLK2/ACP3/MSMB/FOLH1, integrated LE(ARPC)
    ),
    CRPC2 = c(
        `0`  = "Basal DKK1/CCL2",                # KRT5/KRT14/KRT15, SFRP1, integrated BE 1/2/3
        `1`  = "Basal (low-count)",              # 핵 lncRNA (MAGI2/ZBTB20), mt ~11%, integrated BE 4/BE 6
        `2`  = "Basal stress (FOS/JUN)",         # ATF3/DDIT3, integrated BE 5
        `3`  = "Hillock KRT13/OLFM4",            # GABRP/FOXC1/AREG, basal→hillock 이행 (integrated BE 1/2 + Hillock 1/2)
        `4`  = "Hillock squamous (KRT1/SPRR1B)", # KRT13/KRT6A/S100A8, integrated Hillock 2
        `5`  = "Luminal TMPRSS2/PPARG (low-count)", # BMPR1B, integrated OE + LE(ARPC)
        `6`  = "Club SCGB1A1/TFF1",              # SCGB3A1/PIGR/MMP7, integrated Club + LE(ARPC)
        `7`  = "Ionocyte FOXI1",                 # ATP6V0A4/ATP6V0D2/TMEM213, integrated Ionocyte
        `8`  = "Squamous DSG3/KRT6A",            # DSG1/SERPINB13/ESR1, TP63+, mt ~16%
        `9`  = "Luminal AR+ (KLK3)",             # TMPRSS2/FOLH1/ACP3/KLK2 + PCA3/SCHLAP1/AMACR, KRT8 낮음, integrated LE(ARPC)
        `10` = "Ionocyte-like KIT+"              # TFCP2L1/ATP6V0D2, FOXI1 낮음
    ),
    CRPC3 = c(
        `0` = "Basal DKK1/TNC",                  # KRT5/KRT15, CCL2/S100A2, integrated BE 1/2/4
        `1` = "Hillock KRT13/KRT6A",             # LY6D/OLFM4/SERPINB3, integrated Hillock 1/2 + BE 1
        `2` = "Basal stress (FOS/JUN)",          # ATF3/HSPA1A, integrated BE 5
        `3` = "Basal KRT14/POSTN",               # VCAN/DKK2, integrated BE 3
        `4` = "Club MMP7/PIGR",                  # WFDC2/SCGB3A1/LTF, KRT5-, integrated Club
        `5` = "Hillock squamous (KRT13/ESR1)",   # ZNF750/SPRR1B/OVOL1, integrated Hillock 1
        `6` = "Basal cycling (TOP2A)",           # MKI67, integrated BE 1
        `7` = "Luminal AR+ (KLK3)"               # KLK2/FOLH1/ACP3/MSMB, integrated LE(ARPC)
    )
)

annotation_map <- ANNOTATION_MAPS[[PATIENT]]
stopifnot(all(levels(epi$seurat_clusters) %in% names(annotation_map)))
epi$annotation <- factor(unname(annotation_map[as.character(epi$seurat_clusters)]),
                         levels = unique(annotation_map))
Idents(epi) <- "annotation"

p <- utils_umap_by(epi, "annotation", label = TRUE, repel = TRUE, pt.size = 0.3) +
    ggtitle(paste0(PATIENT, " alone — epithelial annotation"))
utils_save_png(p, OUT_DIR, "UMAP_annotation.png", w = 11, h = 8)

p1 <- utils_umap_by(epi, "annotation", label = TRUE, repel = TRUE, pt.size = 0.2) +
    ggtitle("Patient-alone annotation")
p2 <- utils_umap_by(epi, "integrated_epi_annotation", label = TRUE, repel = TRUE, pt.size = 0.2) +
    ggtitle("Integrated annotation")
utils_save_png(p1 + p2, OUT_DIR, "UMAP_annotation_vs_integrated.png", w = 22, h = 8)

utils_save_png(marker_dotplot("annotation"), OUT_DIR, "DotPlot_epi_markers_by_annotation.png",
               w = 30, h = 7)

utils_write_crosstab(epi$annotation, addNA(epi$integrated_epi_annotation),
                     OUT_DIR, "annotation_vs_integrated_annotation.csv")
write.csv(as.data.frame(table(annotation = epi$annotation)),
          file.path(OUT_DIR, "annotation_counts.csv"), row.names = FALSE)
write.csv(
    data.frame(
        cell = colnames(epi),
        cluster = as.character(epi$seurat_clusters),
        annotation = as.character(epi$annotation),
        integrated_epi_annotation = epi$integrated_epi_annotation
    ),
    file.path(OUT_DIR, "annotation_per_cell.csv"),
    row.names = FALSE
)

utils_save_all_markers(epi, file.path(OUT_DIR, "all_markers_by_annotation.csv"))
saveRDS(epi, OBJ_RDS)
message("annotation saved -> ", OBJ_RDS)
print(table(epi$annotation))
