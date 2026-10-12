# Stage 05 — epithelial markers by annotation label (FindAllMarkers → top-N table / DotPlot / heatmap)

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
})
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/ClusterMarkers"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_N_TABLE <- 25
TOP_N_PLOT  <- 5

# FindAllMarkers ----
epi <- readRDS(IN_RDS)
Idents(epi) <- "annotation"  # level 순서 = 표·그림 순서

# SCT assay + PrepSCTFindMarkers, matching utils_save_all_markers() convention.
DefaultAssay(epi) <- "SCT"
epi <- PrepSCTFindMarkers(epi)
all_markers <- FindAllMarkers(
    epi,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.5
)
write.csv(all_markers, file.path(OUT_DIR, "all_markers_by_celltype.csv"),
          row.names = FALSE)

# cell type 별 avg_log2FC 상위 n 개 (p_adj < 0.05)
top_markers <- function(n) {
    all_markers %>%
        filter(p_val_adj < 0.05) %>%
        group_by(cluster) %>%
        slice_max(order_by = avg_log2FC, n = n, with_ties = FALSE) %>%
        ungroup()
}
write.csv(top_markers(TOP_N_TABLE), file.path(OUT_DIR, "top_markers_by_celltype.csv"),
          row.names = FALSE)

# Figures: top-N markers per cell type ----
plot_features <- unique(top_markers(TOP_N_PLOT)$gene)

p_dot <- DotPlot(epi, features = plot_features) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    labs(title = paste0("Top ", TOP_N_PLOT, " markers per epithelial cell type"),
         x = NULL, y = NULL)
utils_save_png(p_dot, OUT_DIR, "top_markers_DotPlot.png",
               w = max(16, length(plot_features) * 0.28), h = 8, limitsize = FALSE)

# 그룹당 300 세포 downsample (큰 cell type 이 heatmap 덮지 않게)
epi_hm <- subset(epi, downsample = 300)
epi_hm <- ScaleData(epi_hm, features = plot_features, verbose = FALSE)
p_hm <- DoHeatmap(epi_hm, features = plot_features,
                  group.colors = utils_cb_cols(Idents(epi_hm))) +
    scale_fill_gradient2(low = "#1F77B4", mid = "grey95", high = "#D7261E") +
    theme(axis.text.y = element_text(size = 6))
utils_save_png(p_hm, OUT_DIR, "top_markers_Heatmap.png",
               w = 14, h = max(10, length(plot_features) * 0.14), limitsize = FALSE)

message("Markers by cell type → ", OUT_DIR)
message("  all_markers_by_celltype.csv  (", nrow(all_markers), " rows)")
message("  top_markers_by_celltype.csv  (top ", TOP_N_TABLE, "/cell type)")
