# Stage 01 — CRPC1-3 integration: decontX → scDblFinder → QC → SCT/Harmony → SingleR → celltype
# combined_CRPC.rds 가 있으면 compute 는 건너뛰고 figures / celltype 만 다시 생성

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(SingleR)
library(celldex)
library(decontX)
library(future)

source("scripts/00_utils/scRNA_utils.R")

# Settings ----
N_CORES <- 30
plan("multicore", workers = N_CORES)
options(future.globals.maxSize = 128 * 1024^3)
# scDblFinder / SingleR (RNGseed: scDblFinder 재현성)
bpp <- MulticoreParam(workers = N_CORES, RNGseed = 42)

SAMPLE_DIRS <- c(
    CRPC1 = "Raw_data/CRPC1/filtered_feature_bc_matrix/",
    CRPC2 = "Raw_data/CRPC2/filtered_feature_bc_matrix/",
    CRPC3 = "Raw_data/CRPC3/filtered_feature_bc_matrix/"
)
CELL_ID_PREFIX <- c("P1", "P2", "P3")  # merge 시 barcode 충돌 방지

# QC cutoffs
MIN_FEATURES <- 200
MAX_FEATURES <- 9000
MAX_MT       <- 20

OUT_DIR <- "Results/01_Integrated"
QC_DIR  <- file.path(OUT_DIR, "QC")
DX_DIR  <- file.path(OUT_DIR, "DecontX")
STR_DIR <- file.path(OUT_DIR, "Stromal")
INT_RDS <- file.path(OUT_DIR, "combined_CRPC.rds")
for (d in c(QC_DIR, DX_DIR, STR_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Compute: decontX → scDblFinder → QC → SCT/Harmony → SingleR ----
# ============================================================

# decontX per sample ----
# ambient leakage 체크용 lineage markers
LEAK_MARKERS <- list(
    Luminal_secretory = c("KLK3", "KLK2", "MSMB", "ACP3"),
    Epithelial        = c("EPCAM", "KRT8", "KRT18"),
    Tcell             = c("CD3D", "CD3E", "TRAC"),
    Myeloid           = c("LYZ", "C1QA", "CD68"),
    Plasma_Ig         = c("IGKC", "IGHG1", "JCHAIN"),
    Endothelial       = c("VWF", "CLDN5"),
    Fibroblast        = c("DCN", "LUM", "COL1A1"),
    SmoothMuscle      = c("ACTA2", "MYH11", "DES"),
    Erythrocyte       = c("HBB", "HBA1", "HBA2")
)

save_decontX_diagnostics <- function(sce, project) {
    p <- plotDecontXContamination(sce) +
        ggtitle(paste0(project, " — decontX contamination"))
    utils_save_png(p, DX_DIR, paste0("contamination_UMAP_", project, ".png"), w = 8, h = 6)

    p <- plotDecontXMarkerPercentage(scater::logNormCounts(sce, assay.type = "counts"),
                                     markers = LEAK_MARKERS, groupClusters = NULL) +
        ggtitle(paste0(project, " — marker %: raw vs decontaminated"))
    utils_save_png(p, DX_DIR, paste0("marker_percentage_", project, ".png"), w = 14, h = 9)
}

# Read10X → decontX → rounded decontXcounts 로 RNA 교체
run_decontX_sample <- function(data_dir, project) {
    seu <- CreateSeuratObject(counts = Read10X(data.dir = data_dir), project = project,
                              min.cells = 3, min.features = 200)

    counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
    sce <- decontX(SingleCellExperiment(assays = list(counts = counts)))
    save_decontX_diagnostics(sce, project)

    decont <- as(round(decontXcounts(sce)), "CsparseMatrix")
    seu[["RNA"]] <- CreateAssay5Object(counts = decont)
    # counts 가 바뀌었으니 QC metrics 재계산
    seu$decontX_contamination <- sce$decontX_contamination
    seu$nCount_RNA   <- colSums(decont)
    seu$nFeature_RNA <- colSums(decont > 0)
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

    message(sprintf(
        "  %s: %d cells, contamination mean=%.3f median=%.3f",
        project, ncol(seu),
        mean(seu$decontX_contamination), median(seu$decontX_contamination)))
    seu
}

# contamination summary (pre-QC)
save_contamination_summary <- function(samples) {
    contam_df <- dplyr::bind_rows(lapply(samples, function(x) {
        data.frame(orig.ident = x$orig.ident,
                   decontX_contamination = x$decontX_contamination)
    }))
    contam_summary <- contam_df |>
        dplyr::group_by(orig.ident) |>
        dplyr::summarise(n = dplyr::n(),
                         mean = mean(decontX_contamination),
                         median = median(decontX_contamination),
                         q90 = quantile(decontX_contamination, 0.90),
                         frac_gt_0.2 = mean(decontX_contamination > 0.2),
                         .groups = "drop")
    write.csv(contam_summary, file.path(DX_DIR, "contamination_summary_by_sample.csv"),
              row.names = FALSE)

    p <- ggplot(contam_df, aes(orig.ident, decontX_contamination, fill = orig.ident)) +
        geom_violin(scale = "width", trim = TRUE) +
        geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "white") +
        scale_fill_manual(values = utils_cb_cols(contam_df$orig.ident)) +
        labs(x = NULL, y = "decontX contamination fraction",
             title = "Estimated ambient RNA contamination per sample") +
        theme_bw() + theme(legend.position = "none")
    utils_save_png(p, DX_DIR, "contamination_violin_by_sample.png", w = 8, h = 6)
}

# scDblFinder per sample (merge 전) ----
run_scDblFinder <- function(seu) {
    sce <- as.SingleCellExperiment(seu)
    sce <- scDblFinder(sce, BPPARAM = bpp)
    seu$scDblFinder.class <- sce$scDblFinder.class
    seu$scDblFinder.score <- sce$scDblFinder.score
    subset(seu, subset = scDblFinder.class == "singlet")
}

# QC plots (doublet 제거 후, filter 전) ----
save_qc_plots <- function(obj) {
    patient_cols <- utils_cb_cols(obj$orig.ident)

    p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "orig.ident", ncol = 3, pt.size = 0, cols = patient_cols)
    utils_save_png(p, QC_DIR, "mt_by_patient_VlnPlot.png", w = 12, h = 8)

    p <- VlnPlot(obj, features = "decontX_contamination",
                 group.by = "orig.ident", pt.size = 0, cols = patient_cols) + NoLegend()
    utils_save_png(p, QC_DIR, "decontX_contamination_VlnPlot.png", w = 8, h = 6)

    cutoffs <- list(
        nFeature_RNA = c(MIN_FEATURES, MAX_FEATURES),
        nCount_RNA   = NULL,
        percent.mt   = MAX_MT
    )
    vlns <- lapply(names(cutoffs), function(feat) {
        p <- VlnPlot(obj, features = feat, group.by = "orig.ident",
                     pt.size = 0, cols = patient_cols) + NoLegend()
        for (y in cutoffs[[feat]]) {
            p <- p + geom_hline(yintercept = y, linetype = "dashed", color = "red")
        }
        p
    })
    utils_save_png(wrap_plots(vlns, ncol = 3), QC_DIR, "QC_VlnPlot_with_thresholds.png", w = 14, h = 6)

    per_patient <- SplitObject(obj, split.by = "orig.ident")
    scatter_rows <- lapply(names(per_patient), function(s) {
        x  <- per_patient[[s]]
        p1 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle(s)
        p2 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(s)
        p1 + p2
    })
    utils_save_png(wrap_plots(scatter_rows, ncol = 1), QC_DIR, "FeatureScatter_by_patient.png",
             w = 12, h = 15)
}

# SingleR (HPCA main) ----
run_singleR <- function(obj) {
    ref  <- celldex::HumanPrimaryCellAtlasData()
    pred <- SingleR(
        test = as.SingleCellExperiment(obj, assay = "SCT"),
        ref = ref,
        labels = ref$label.main,
        BPPARAM = bpp
    )
    obj$SingleR        <- pred$labels
    obj$SingleR_pruned <- pred$pruned.labels

    write.csv(data.frame(cell = rownames(pred),
                         label = pred$labels,
                         pruned_label = pred$pruned.labels),
              file.path(OUT_DIR, "SingleR_cell_labels.csv"), row.names = FALSE)
    png(file.path(OUT_DIR, "SingleR_score_heatmap.png"), width = 1200, height = 900, bg = "white")
    plotScoreHeatmap(pred)
    dev.off()
    png(file.path(OUT_DIR, "SingleR_delta_distribution.png"), width = 1200, height = 900, bg = "white")
    print(plotDeltaDistribution(pred, ncol = 4))  # print 없으면 빈 png
    dev.off()
    obj
}

# full compute → checkpoint (celltype 없음) ----
build_integrated_object <- function() {
    samples <- Map(run_decontX_sample, SAMPLE_DIRS, names(SAMPLE_DIRS))
    save_contamination_summary(samples)

    set.seed(42)
    samples <- lapply(samples, run_scDblFinder)

    merged <- merge(samples[[1]], y = unname(samples[-1]),
                    add.cell.ids = CELL_ID_PREFIX, project = "CRPC")
    save_qc_plots(merged)

    obj <- subset(merged, subset = nFeature_RNA > MIN_FEATURES &
                                   nFeature_RNA < MAX_FEATURES &
                                   percent.mt < MAX_MT)

    obj <- utils_sct_harmony(obj)
    utils_save_png(ElbowPlot(obj, ndims = 50), QC_DIR, "ElbowPlot.png", w = 12, h = 15)

    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
    obj <- FindClusters(obj, resolution = 0.5)
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)

    obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes,
                            g2m.features = cc.genes.updated.2019$g2m.genes)
    obj <- run_singleR(obj)

    utils_save_all_markers(obj, file.path(OUT_DIR, "all_markers.csv"))
    saveRDS(obj, INT_RDS)
    obj
}

if (file.exists(INT_RDS)) {
    message("Loading cached ", INT_RDS, " — skipping compute")
    combined_CRPC <- readRDS(INT_RDS)
} else {
    combined_CRPC <- build_integrated_object()
}

# ============================================================
# Figures ----
# ============================================================
p <- utils_umap_by(combined_CRPC, "Phase", pt.size = 0.3)
utils_save_png(p, OUT_DIR, "CellCycle_Phase_UMAP.png", w = 10, h = 7)

# names = genes, values = panel titles
LINEAGE_MARKERS <- c(
    EPCAM  = "Epithelial",
    KRT8   = "Epithelial (luminal)",
    AR     = "Luminal / AR pathway",
    PTPRC  = "Immune",
    VIM    = "Mesenchymal",
    PECAM1 = "Endothelial",
    ACTA2  = "Smooth muscle (ACTA2)",
    MYH11  = "Smooth muscle (MYH11)",
    SYP    = "Neuroendocrine",
    COL1A1 = "Fibroblast",
    DCN    = "Fibroblast / stromal",
    TPSAB1 = "Mast cell"
)
fp <- FeaturePlot(combined_CRPC, features = names(LINEAGE_MARKERS), pt.size = 0.1,
                  combine = FALSE)
fp <- Map(function(p, gene, lineage) p + ggtitle(bquote(.(lineage) ~ "—" ~ italic(.(gene)))),
          fp, names(LINEAGE_MARKERS), LINEAGE_MARKERS)
utils_save_png(wrap_plots(fp, ncol = 4), OUT_DIR, "FeaturePlot_markers.png", w = 20, h = 19)

# cluster 번호 확인용 (celltype map)
p <- utils_umap_by(combined_CRPC, "seurat_clusters", label = TRUE)
utils_save_png(p, OUT_DIR, "Cluster_UMAP_integrated.png", w = 15, h = 15)
p <- utils_umap_by(combined_CRPC, "seurat_clusters", split.by = "orig.ident", label = TRUE)
utils_save_png(p, OUT_DIR, "Cluster_UMAP_by_patient.png", w = 24, h = 15)

p <- FeaturePlot(combined_CRPC, features = "decontX_contamination", pt.size = 0.2) +
    viridis::scale_color_viridis(option = "viridis") +
    ggtitle("decontX contamination (integrated UMAP)")
utils_save_png(p, OUT_DIR, "decontX_contamination_UMAP.png", w = 10, h = 8)

p <- utils_umap_by(combined_CRPC, "SingleR_pruned", label = TRUE, repel = TRUE)
utils_save_png(p, OUT_DIR, "SingleR_UMAP.png", w = 15, h = 15)
p <- utils_umap_by(combined_CRPC, "SingleR", split.by = "orig.ident", label = TRUE, repel = TRUE)
utils_save_png(p, OUT_DIR, "SingleR_UMAP_by_patient.png", w = 24, h = 15)

# ============================================================
# Celltype assignment (res 0.5, 19 clusters) ----
# ============================================================
# 실행환경(아마 N_CORES)이 달라지면 cluster 번호가 바뀔 수 있음.
# 직접 실행하기 보다는 combined_CRPC.rds를 로드해 사용할 것을 권장.
CELLTYPE_MAP <- c(
    `0`  = "Epithelial",          `1`  = "Fibroblast",
    `2`  = "Epithelial",          `3`  = "Epithelial",
    `4`  = "Epithelial",          `5`  = "Epithelial",
    `6`  = "T/NK cells",          `7`  = "Epithelial",
    `8`  = "Endothelial",         `9`  = "Epithelial",
    `10` = "Epithelial",          `11` = "Smooth muscle cells",
    `12` = "Smooth muscle cells", `13` = "Phagocytes",
    `14` = "T/NK cells",          `15` = "Epithelial",
    `16` = "Epithelial",          `17` = "Mast cells",
    `18` = "Epithelial"
)
combined_CRPC$celltype <- factor(unname(CELLTYPE_MAP[as.character(combined_CRPC$seurat_clusters)]))
Idents(combined_CRPC) <- "celltype"

# ============================================================
# Stromal refinement ----
# ============================================================
# cl 11/12에 마커 유전자로 보아 pericyte / vSMC / fibroblast / Schwann이 혼합된 것으로 보임.
# stromal 세포만 subset → 재클러스터링 → subcluster annotation → combined_CRPC celltype 교체
STR_RES_VEC <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2)  # UMAP_multires.png 용
STR_RES     <- 0.3                                         # annotation 해상도

recluster_stromal <- function(obj) {
    str_obj <- subset(obj, subset = celltype %in% c("Fibroblast", "Smooth muscle cells"))
    message("Stromal subset: ", ncol(str_obj), " cells")
    str_obj$parent_celltype <- as.character(str_obj$celltype)
    str_obj$parent_cluster  <- as.character(str_obj$seurat_clusters)

    # Reset assay state (drop inherited SCT, re-split RNA by patient for Harmony)
    DefaultAssay(str_obj) <- "RNA"
    str_obj[["SCT"]] <- NULL
    str_obj[["RNA"]] <- split(str_obj[["RNA"]], f = str_obj$orig.ident)
    str_obj@meta.data <- dplyr::select(str_obj@meta.data,
        -matches("_snn_res\\."), -seurat_clusters)

    str_obj <- utils_sct_harmony(str_obj)
    str_obj <- FindNeighbors(str_obj, reduction = "harmony", dims = 1:30)
    str_obj <- FindClusters(str_obj, resolution = STR_RES_VEC)
    str_obj$stromal_subcluster <- str_obj@meta.data[[paste0("SCT_snn_res.", STR_RES)]]
    RunUMAP(str_obj, reduction = "harmony", dims = 1:30,
            n.neighbors = 20, min.dist = 0.1, spread = 4.0)
}

save_stromal_diagnostics <- function(str_obj) {
    p <- utils_umap_by(str_obj, "stromal_subcluster", label = TRUE, pt.size = 0.4) +
        ggtitle(paste0("Stromal subclusters (res ", STR_RES, ")"))
    utils_save_png(p, STR_DIR, "UMAP_subcluster.png", w = 10, h = 8)
    p <- utils_umap_by(str_obj, "celltype", label = TRUE, pt.size = 0.4) +
        ggtitle("Stromal cell types")
    utils_save_png(p, STR_DIR, "UMAP_celltype.png", w = 10, h = 8)
    p <- utils_umap_by(str_obj, "celltype", split.by = "orig.ident", label = TRUE, pt.size = 0.4)
    utils_save_png(p, STR_DIR, "UMAP_celltype_by_patient.png", w = 24, h = 8)
    p <- utils_umap_by(str_obj, "parent_cluster", label = TRUE, pt.size = 0.4) +
        ggtitle("Parent res-0.5 cluster (1 = Fibroblast, 11/12 = coarse 'SMC')")
    utils_save_png(p, STR_DIR, "UMAP_parent_cluster.png", w = 10, h = 8)

    multires <- lapply(STR_RES_VEC, function(r) {
        utils_umap_by(str_obj, paste0("SCT_snn_res.", r), label = TRUE, pt.size = 0.3) +
            ggtitle(paste0("res = ", r)) + NoLegend()
    })
    utils_save_png(wrap_plots(multires, ncol = 3), STR_DIR, "UMAP_multires.png",
             w = 18, h = 14, dpi = 200)

    p <- VlnPlot(str_obj, group.by = "stromal_subcluster", pt.size = 0, ncol = 2,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination"),
                 cols = utils_cb_cols(str_obj$stromal_subcluster))
    utils_save_png(p, STR_DIR, "QC_violin_by_subcluster.png", w = 16, h = 10)

    stromal_markers <- list(
        Fibroblast    = c("DCN", "LUM", "PDGFRA", "COL1A1", "FBLN1", "APOD"),
        Pericyte      = c("RGS5", "NOTCH3", "HIGD1B", "NDUFA4L2", "KCNJ8", "GJA4"),
        `Smooth muscle cells` = c("MYH11", "ACTG2", "DES", "MYLK", "CNN1", "ACTA2"),
        Schwann       = c("PLP1", "S100B", "MPZ", "SOX10", "CDH19"),
        Contamination = c("PECAM1", "PTPRC", "EPCAM", "KRT8")
    )
    marker_dotplot <- function(group, title) {
        DotPlot(str_obj, features = stromal_markers, group.by = group,
                cluster.idents = FALSE) +
            RotatedAxis() +
            scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
            labs(title = title)
    }
    utils_save_png(marker_dotplot("stromal_subcluster", "Stromal markers by subcluster"),
             STR_DIR, "DotPlot_by_subcluster.png", w = 20, h = 8)
    utils_save_png(marker_dotplot("celltype", "Stromal markers by cell type"),
             STR_DIR, "DotPlot_by_celltype.png", w = 18, h = 6)

    utils_write_crosstab(str_obj$parent_cluster, str_obj$stromal_subcluster,
                   STR_DIR, "parent_vs_subcluster.csv")
    utils_write_crosstab(str_obj$stromal_subcluster, str_obj$celltype,
                   STR_DIR, "subcluster_to_celltype.csv")
    utils_write_crosstab(str_obj$celltype, str_obj$orig.ident,
                   STR_DIR, "celltype_by_patient.csv")
    qc_tab <- str_obj@meta.data |>
        group_by(subcluster = stromal_subcluster) |>
        summarise(n = n(),
                  median_nCount = median(nCount_RNA),
                  median_nFeature = median(nFeature_RNA),
                  median_mt = round(median(percent.mt), 2),
                  median_decontX = round(median(decontX_contamination), 3),
                  .groups = "drop")
    write.csv(qc_tab, file.path(STR_DIR, "subcluster_QC_summary.csv"), row.names = FALSE)

    Idents(str_obj) <- "stromal_subcluster"
    utils_save_all_markers(str_obj, file.path(STR_DIR, "all_markers_by_subcluster.csv"))
    invisible(NULL)  # autoprint 방지
}

str_obj <- recluster_stromal(combined_CRPC)

# res 0.3 subcluster → label.
#   c6  pericyte (HIGD1B/NDUFA4L2/GJA4/ADGRF5)
#   c4  vascular SMC (PLN/RERGL/CASQ2, MYH11+ 99.8%, NOTCH3+ 65%) → Pericyte 로 묶음
#   c9  분화 SM (DES/ACTG2/MYLK/CNN1),  c8  MYH11+ 저복잡도 SM
#   c10 Schwann (PLP1/S100B/SOX10/CDH19)
#   c0 APOD+ / c1 ASPN·OGN·C7 / c3 PENK·PTGDS / c7 CA3·APOE fibroblast
# rerun 시 subcluster 번호가 바뀔 수 있음 → Stromal/ 그림 보고 map 을 확인·수정할 것
STROMAL_MAP <- c(
    `0` = "Fibroblast", `1` = "Fibroblast", `2` = "Fibroblast",
    `3` = "Fibroblast", `7` = "Fibroblast",
    `4` = "Pericyte",   `5` = "Pericyte",   `6` = "Pericyte",
    `8` = "Smooth muscle cells", `9` = "Smooth muscle cells",
    `10` = "Schwann cell"
)
STROMAL_LEVELS  <- c("Fibroblast", "Pericyte", "Smooth muscle cells", "Schwann cell")
CELLTYPE_LEVELS <- c("Epithelial", "Fibroblast", "Pericyte", "Smooth muscle cells",
                     "Schwann cell", "Endothelial", "T/NK cells", "Phagocytes", "Mast cells")

str_obj$celltype <- factor(unname(STROMAL_MAP[as.character(str_obj$stromal_subcluster)]),
                           levels = STROMAL_LEVELS)
Idents(str_obj) <- "celltype"

# stromal cells 만 celltype 교체
ct <- as.character(combined_CRPC$celltype)
ct[match(colnames(str_obj), colnames(combined_CRPC))] <- as.character(str_obj$celltype)
combined_CRPC$celltype <- factor(ct, levels = CELLTYPE_LEVELS)
Idents(combined_CRPC) <- "celltype"
message("Stromal refinement applied:")
print(table(str_obj$parent_celltype, str_obj$celltype))

save_stromal_diagnostics(str_obj)
saveRDS(str_obj, file.path(STR_DIR, "stromal_annotated.rds"))
rm(str_obj)
invisible(gc())

# ============================================================
# Final figures + save ----
# ============================================================
p <- utils_umap_by(combined_CRPC, "celltype", label = TRUE)
utils_save_png(p, OUT_DIR, "Labelled_UMAP_integrated.png", w = 15, h = 15)
p <- utils_umap_by(combined_CRPC, "celltype", split.by = "orig.ident", label = TRUE)
utils_save_png(p, OUT_DIR, "Labelled_UMAP_by_patient.png", w = 24, h = 15)

comp <- as.data.frame(table(celltype = combined_CRPC$celltype,
                            patient = combined_CRPC$orig.ident))
p <- ggplot(comp, aes(x = patient, y = Freq, fill = celltype)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = utils_cb_cols(combined_CRPC$celltype)) +
    labs(y = "fraction of cells", x = NULL, title = "Cell type composition by patient") +
    theme_classic(base_size = 14)
utils_save_png(p, OUT_DIR, "Celltype_composition_by_patient.png", w = 8, h = 6)
utils_write_crosstab(combined_CRPC$celltype, combined_CRPC$orig.ident,
               OUT_DIR, "celltype_by_patient.csv")

saveRDS(combined_CRPC, INT_RDS)  # stage 02 입력
message("celltype saved -> ", INT_RDS)
print(table(combined_CRPC$celltype))
