# Stage 03 — 상피 cluster 별 QC, lineage purity 진단 후 filter_decision.csv 작성
# filter_decision.csv에서 keep = FALSE로 설정된 cluster는 이후 단계에 포함되지 않음.
# 저품질 cluster를 제거하기 위한 스크립트였으나 decontX 적용 후에는 모두 적합해 전부 keep=TRUE로 설정함.

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

source("scripts/00_utils/scRNA_utils.R")

# Settings ----
IN_RDS       <- "Results/02_Epithelial_Initial/epi_clustered.rds"
OUT_DIR      <- "Results/03_Epithelial_FilterDecision"
DECISION_CSV <- file.path(OUT_DIR, "filter_decision.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

epi <- readRDS(IN_RDS)
cluster_cols <- utils_cb_cols(epi$seurat_clusters)

# Cluster-level QC metrics ----
qc_vlns <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(feat) {
    VlnPlot(epi, features = feat, pt.size = 0, cols = cluster_cols) +
        NoLegend() + labs(x = "Cluster", title = feat)
})
utils_save_png(wrap_plots(qc_vlns, ncol = 1), OUT_DIR, "cluster_QC_violin.png",
               w = 12, h = 12)

# Lineage purity — epi vs stromal/immune/endothelial ----
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
    labs(title = "Lineage purity — epi/stromal/immune/endo markers per cluster")
utils_save_png(p, OUT_DIR, "cluster_lineage_DotPlot.png", w = 18, h = 8)

key_markers <- c(
    "EPCAM", "KRT8", "KRT5", "TP63",
    "COL1A1", "ACTA2", "RGS5",
    "PECAM1", "PTPRC"
)
p <- FeaturePlot(epi, features = key_markers, ncol = 4, pt.size = 0.1)
utils_save_png(p, OUT_DIR, "cluster_lineage_FeaturePlot.png", w = 20, h = 14)

# Cluster × patient cross-table ----
patient_n <- dplyr::count(epi@meta.data, seurat_clusters, orig.ident, name = "n")

patient_table <- patient_n %>%
    tidyr::pivot_wider(names_from = orig.ident, values_from = n, values_fill = 0)
patient_table$total <- rowSums(dplyr::select(patient_table, -seurat_clusters))
write.csv(patient_table, file.path(OUT_DIR, "cluster_patient_table.csv"), row.names = FALSE)

# filter_decision.csv (템플릿 생성) ----
# 이미 있으면 덮어쓰지 않음
if (file.exists(DECISION_CSV)) {
    message(DECISION_CSV, " already exists; leaving it unchanged")
} else {
    filter_decision <- data.frame(
        cluster = levels(epi$seurat_clusters),
        keep    = TRUE,
        note    = NA_character_
    )
    write.csv(filter_decision, DECISION_CSV, row.names = FALSE)
}

message("Stage 03 figures/tables -> ", OUT_DIR)
