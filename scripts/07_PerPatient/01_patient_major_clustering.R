# 07-01 - integration 없이 환자별 clustering, major annotation
# decontX, scDblFinder, QC 는 stage 01 에서 sample 별로 처리함 -> combined_CRPC.rds 의 RNA count 를 환자별로 불러와 사용
# <PATIENT>_major.rds 가 있으면 compute 는 건너뛰고 figures / celltype 만 다시 생성

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(clustree)
library(SingleCellExperiment)
library(BiocParallel)
library(SingleR)
library(celldex)
library(future)

source("scripts/00_utils/scRNA_utils.R")

# Settings ----
PATIENT <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(PATIENT %in% c("CRPC1", "CRPC2", "CRPC3"))

N_CORES <- 30
plan("multicore", workers = N_CORES)
options(future.globals.maxSize = 128 * 1024^3)

IN_RDS  <- "Results/01_Integrated/combined_CRPC.rds"
OUT_DIR <- file.path("Results/07_PerPatient", PATIENT, "01_Major")
OBJ_RDS <- file.path(OUT_DIR, paste0(PATIENT, "_major.rds"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

RES_VEC   <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)  # clustree, UMAP_multires 용
RES_MAJOR <- "0.5"                                      # seurat_clusters 해상도

# ============================================================
# Compute: patient subset → SCT → clustering → UMAP → SingleR ----
# ============================================================
# integrated anlysis에서는 count 와 QC 메타데이터만 가져옴
subset_patient <- function() {
    combined <- readRDS(IN_RDS)
    cells  <- colnames(combined)[combined$orig.ident == PATIENT]
    counts <- GetAssayData(combined, assay = "RNA", layer = "counts")[, cells]
    counts <- counts[Matrix::rowSums(counts) > 0, ]
    md <- combined@meta.data[cells, c("orig.ident", "percent.mt", "decontX_contamination",
                                      "scDblFinder.score", "celltype")]
    md <- dplyr::rename(md, integrated_celltype = celltype)
    CreateSeuratObject(counts = counts, meta.data = md, project = PATIENT)
}

build_major_object <- function() {
    obj <- subset_patient()
    message(PATIENT, ": ", ncol(obj), " cells x ", nrow(obj), " genes")

    obj <- utils_sct_pca(obj)
    utils_save_png(ElbowPlot(obj, ndims = 50), OUT_DIR, "ElbowPlot.png", w = 10, h = 6)

    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
    obj <- FindClusters(obj, resolution = RES_VEC, random.seed = 42)
    obj$seurat_clusters <- obj@meta.data[[paste0("SCT_snn_res.", RES_MAJOR)]]
    Idents(obj) <- "seurat_clusters"
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

    obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes,
                            g2m.features = cc.genes.updated.2019$g2m.genes, nbin = 12)

    # SingleR (HPCA main)
    ref  <- celldex::HumanPrimaryCellAtlasData()
    pred <- SingleR(test = as.SingleCellExperiment(obj, assay = "SCT"), ref = ref,
                    labels = ref$label.main, BPPARAM = MulticoreParam(workers = N_CORES))
    obj$SingleR        <- pred$labels
    obj$SingleR_pruned <- pred$pruned.labels

    utils_save_all_markers(obj, file.path(OUT_DIR, "all_markers.csv"))
    saveRDS(obj, OBJ_RDS)
    obj
}

if (file.exists(OBJ_RDS)) {
    message("Loading cached ", OBJ_RDS, " — skipping compute")
    obj <- readRDS(OBJ_RDS)
} else {
    obj <- build_major_object()
}

# ============================================================
# Figures ----
# ============================================================
p <- clustree(obj, prefix = "SCT_snn_res.")
utils_save_png(p, OUT_DIR, "clustree.png", w = 12, h = 14, dpi = 200)

multires <- lapply(RES_VEC, function(r) {
    utils_umap_by(obj, paste0("SCT_snn_res.", r), label = TRUE, pt.size = 0.2) +
        ggtitle(paste0("res = ", r)) + NoLegend()
})
utils_save_png(wrap_plots(multires, ncol = 4), OUT_DIR, "UMAP_multires.png",
               w = 24, h = 12, dpi = 200)

p <- utils_umap_by(obj, "seurat_clusters", label = TRUE, pt.size = 0.3) +
    ggtitle(paste0(PATIENT, " alone — clusters (res ", RES_MAJOR, ")"))
utils_save_png(p, OUT_DIR, "UMAP_cluster.png", w = 10, h = 8)

p <- utils_umap_by(obj, "integrated_celltype", label = TRUE, repel = TRUE, pt.size = 0.3) +
    ggtitle("Integrated (Harmony) celltype projected on patient-alone UMAP")
utils_save_png(p, OUT_DIR, "UMAP_integrated_celltype.png", w = 11, h = 8)

p <- utils_umap_by(obj, "SingleR_pruned", label = TRUE, repel = TRUE, pt.size = 0.3) +
    ggtitle("SingleR (HPCA main, patient alone)")
utils_save_png(p, OUT_DIR, "UMAP_SingleR.png", w = 12, h = 8)

p <- utils_umap_by(obj, "Phase", pt.size = 0.3)
utils_save_png(p, OUT_DIR, "UMAP_CellCycle.png", w = 10, h = 7)

p <- FeaturePlot(obj, features = "decontX_contamination", pt.size = 0.2) +
    viridis::scale_color_viridis()
utils_save_png(p, OUT_DIR, "UMAP_decontX.png", w = 10, h = 8)

p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination"),
             group.by = "seurat_clusters", pt.size = 0, ncol = 2,
             cols = utils_cb_cols(obj$seurat_clusters))
utils_save_png(p, OUT_DIR, "QC_violin_by_cluster.png", w = 16, h = 10)

# Lineage marker ----
LINEAGE_MARKERS <- list(
    Epithelial  = c("EPCAM", "KRT8", "KRT18", "KRT5", "KRT15", "AR", "KLK3"),
    NE          = c("CHGA", "SYP", "ENO2"),
    Fibroblast  = c("DCN", "LUM", "COL1A1", "PDGFRA"),
    Pericyte    = c("RGS5", "NOTCH3", "HIGD1B"),
    SMC         = c("MYH11", "ACTG2", "DES", "CNN1"),
    Schwann     = c("PLP1", "S100B", "SOX10"),
    Endothelial = c("PECAM1", "VWF", "CLDN5"),
    Lymphatic   = c("PROX1", "LYVE1"),
    T           = c("CD3D", "CD3E", "CD8A", "CD4", "IL7R"),
    NK          = c("NKG7", "GNLY", "KLRD1"),
    B           = c("MS4A1", "CD79A"),
    Plasma      = c("JCHAIN", "MZB1", "IGKC"),
    Myeloid     = c("LYZ", "CD68", "C1QA", "CD14", "FCGR3A"),
    DC          = c("CLEC9A", "CD1C", "LAMP3"),
    Mast        = c("TPSAB1", "CPA3", "KIT"),
    Cycling     = c("MKI67", "TOP2A")
)

LINEAGE_MARKERS <- lapply(LINEAGE_MARKERS, intersect, rownames(obj))

marker_dotplot <- function(group, title = NULL) {
    DotPlot(obj, features = LINEAGE_MARKERS, group.by = group, cluster.idents = FALSE) +
        RotatedAxis() +
        scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
        theme(axis.text.x = element_text(size = 8)) +
        labs(title = title)
}
p <- marker_dotplot("seurat_clusters", paste0(PATIENT, " alone — lineage markers by cluster"))
utils_save_png(p, OUT_DIR, "DotPlot_lineage_by_cluster.png", w = 26, h = 8)

FP_GENES <- c("EPCAM", "KRT8", "KRT5", "AR", "KLK3", "PTPRC", "CD3E", "LYZ",
              "PECAM1", "DCN", "RGS5", "MYH11", "TPSAB1", "JCHAIN", "SYP", "MKI67")
p <- FeaturePlot(obj, features = FP_GENES, ncol = 4, pt.size = 0.1)
utils_save_png(p, OUT_DIR, "FeaturePlot_lineage.png", w = 20, h = 20)

# SingleR, Integrated anlysis의 annotation과 비교 ----
utils_write_crosstab(obj$seurat_clusters, obj$integrated_celltype,
                     OUT_DIR, "cluster_vs_integrated_celltype.csv")
utils_write_crosstab(obj$seurat_clusters, obj$SingleR, OUT_DIR, "cluster_vs_SingleR.csv")

avg <- AverageExpression(obj, features = unlist(LINEAGE_MARKERS), group.by = "seurat_clusters",
                         assay = "SCT", layer = "data")$SCT
write.csv(round(as.matrix(avg), 3), file.path(OUT_DIR, "lineage_marker_avg_by_cluster.csv"))

qc_tab <- obj@meta.data |>
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
# Celltype assignment (res 0.5) ----
# ============================================================
# cluster별로 대표적인 마커 확인해 major cell type 배정
CELLTYPE_MAPS <- list(
    CRPC1 = c(
        `0`  = "Epithelial",
        `1`  = "Epithelial",
        `2`  = "Epithelial",
        `3`  = "Fibroblast",
        `4`  = "T/NK cells",
        `5`  = "Epithelial",
        `6`  = "Endothelial",
        `7`  = "Epithelial",
        `8`  = "Epithelial",
        `9`  = "Epithelial",
        `10` = "Epithelial",
        `11` = "Epithelial",
        `12` = "T/NK cells",
        `13` = "Pericyte",
        `14` = "Pericyte",
        `15` = "T/NK cells",
        `16` = "Phagocytes",
        `17` = "Smooth muscle cells",
        `18` = "Mast cells",
        `19` = "Phagocytes",
        `20` = "Erythrocytes"
    ),
    CRPC2 = c(
        `0`  = "Epithelial",
        `1`  = "Fibroblast",
        `2`  = "Epithelial",
        `3`  = "Epithelial",
        `4`  = "Epithelial",
        `5`  = "Endothelial",
        `6`  = "T/NK cells",
        `7`  = "Fibroblast",
        `8`  = "Epithelial",
        `9`  = "Epithelial",
        `10` = "Pericyte",
        `11` = "Smooth muscle cells",
        `12` = "Mast cells",
        `13` = "Epithelial",
        `14` = "Phagocytes",
        `15` = "Epithelial",
        `16` = "Schwann cell",
        `17` = "Epithelial",
        `18` = "Epithelial"
    ),
    CRPC3 = c(
        `0`  = "Epithelial",
        `1`  = "Epithelial",
        `2`  = "Epithelial",
        `3`  = "Epithelial",
        `4`  = "Fibroblast",
        `5`  = "Epithelial",
        `6`  = "Epithelial",
        `7`  = "Fibroblast",
        `8`  = "T/NK cells",          
        `9`  = "Epithelial",
        `10` = "Endothelial",
        `11` = "Low-quality",         # 대표 마커 없음, low count
        `12` = "Phagocytes",
        `13` = "Pericyte",
        `14` = "Phagocytes",
        `15` = "Epithelial"
    )
)
CELLTYPE_LEVELS <- c("Epithelial", "Fibroblast", "Pericyte", "Smooth muscle cells", "Schwann cell",
                     "Endothelial", "T/NK cells", "Phagocytes", "Mast cells", "Erythrocytes",
                     "Low-quality")

celltype_map <- CELLTYPE_MAPS[[PATIENT]]
stopifnot(all(levels(obj$seurat_clusters) %in% names(celltype_map)))
obj$celltype <- droplevels(factor(unname(celltype_map[as.character(obj$seurat_clusters)]),
                                  levels = CELLTYPE_LEVELS))
Idents(obj) <- "celltype"

p <- utils_umap_by(obj, "celltype", label = TRUE, repel = TRUE, pt.size = 0.3) +
    ggtitle(paste0(PATIENT, " alone — cell types"))
utils_save_png(p, OUT_DIR, "UMAP_celltype.png", w = 11, h = 8)

p1 <- utils_umap_by(obj, "celltype", label = TRUE, repel = TRUE, pt.size = 0.2) +
    ggtitle("Patient-alone annotation")
p2 <- utils_umap_by(obj, "integrated_celltype", label = TRUE, repel = TRUE, pt.size = 0.2) +
    ggtitle("Integrated annotation")
utils_save_png(p1 + p2, OUT_DIR, "UMAP_celltype_vs_integrated.png", w = 22, h = 8)

utils_save_png(marker_dotplot("celltype"), OUT_DIR, "DotPlot_lineage_by_celltype.png", w = 26, h = 7)

utils_write_crosstab(obj$celltype, obj$integrated_celltype,
                     OUT_DIR, "celltype_vs_integrated_celltype.csv")
write.csv(as.data.frame(table(celltype = obj$celltype)),
          file.path(OUT_DIR, "celltype_counts.csv"), row.names = FALSE)
write.csv(
    data.frame(
        cell = colnames(obj),
        cluster = as.character(obj$seurat_clusters),
        celltype = as.character(obj$celltype),
        integrated_celltype = as.character(obj$integrated_celltype)
    ),
    file.path(OUT_DIR, "celltype_per_cell.csv"),
    row.names = FALSE
)

saveRDS(obj, OBJ_RDS)  # 07-02 입력
message("celltype saved -> ", OBJ_RDS)
print(table(obj$celltype))
