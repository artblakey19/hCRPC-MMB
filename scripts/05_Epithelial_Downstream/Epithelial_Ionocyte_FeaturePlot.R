# Stage 05 — Ionocyte 후보 마커 FeaturePlot + violin (Ionocyte = cl11 라벨 확인용)

suppressMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
})
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Annotation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MARKERS <- c("FOXI1", "ATP6V1G3", "ATP6V1B1", "KRT7", "KIT")

epi <- readRDS(IN_RDS)
# RNA 는 counts 만 → SCT data 사용
DefaultAssay(epi) <- "SCT"

# FeaturePlot ----
fp_list <- lapply(MARKERS, function(g) {
    FeaturePlot(epi, features = g, order = TRUE, pt.size = 0.25) +
        ggtitle(g) +
        theme(plot.title = element_text(face = "bold"))
})
utils_save_png(wrap_plots(fp_list, ncol = 3), OUT_DIR, "FeaturePlot_ionocyte_markers.png",
               w = 15, h = 10, dpi = 200)

# VlnPlot ----
# stack = TRUE 는 상수(전부 0) 유전자 있으면 에러 → violin 에서만 제외
# ATP6V1B1: SCT data 전부 0 (FeaturePlot 은 회색 패널)
sct_data <- GetAssayData(epi, assay = "SCT", layer = "data")[MARKERS, ]
vln_markers <- MARKERS[apply(sct_data, 1, function(v) length(unique(v)) > 1)]
message("Excluded from VlnPlot: ", paste(setdiff(MARKERS, vln_markers), collapse = ", "))

vp <- VlnPlot(epi, features = vln_markers, group.by = "annotation",
              pt.size = 0, stack = TRUE, flip = TRUE) +
    NoLegend() +
    ggtitle("Ionocyte candidate markers by annotation")
utils_save_png(vp, OUT_DIR, "VlnPlot_ionocyte_markers.png", w = 9, h = 7, dpi = 200)
