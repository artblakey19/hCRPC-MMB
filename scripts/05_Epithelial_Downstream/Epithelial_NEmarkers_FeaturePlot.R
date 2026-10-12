# Stage 05 — NE 마커 patient-split FeaturePlot (DNPC NE-negativity 확인; SCT)

suppressMessages({
    library(Seurat); library(ggplot2); library(patchwork)
})
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
EPI_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Annotation"
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 표준 neuroendocrine/NEPC 패널
NE_GENES <- c("CHGA", "CHGB", "SYP", "NCAM1", "ASCL1")

epi <- readRDS(EPI_RDS)
DefaultAssay(epi) <- "SCT"

# 발현 세포 수 요약 (콘솔) ----
expr_all <- FetchData(epi, vars = c(NE_GENES, "orig.ident"), layer = "data")
message("=== Cells expressing NE markers (>0) / ", ncol(epi), " ===")
for (g in NE_GENES) {
    tb <- tapply(expr_all[[g]] > 0, expr_all$orig.ident, sum)
    message(sprintf("  %-6s total=%4d (%.2f%%)  CRPC1=%d CRPC2=%d CRPC3=%d",
                    g, sum(expr_all[[g]] > 0), 100 * mean(expr_all[[g]] > 0),
                    tb[["CRPC1"]], tb[["CRPC2"]], tb[["CRPC3"]]))
}

# patient-split FeaturePlot ----
# 모든 패널 공통 UMAP 프레임
emb <- Embeddings(epi, "umap")
xr  <- range(emb[, 1]) + 0.03 * c(-1, 1) * diff(range(emb[, 1]))
yr  <- range(emb[, 2]) + 0.03 * c(-1, 1) * diff(range(emb[, 2]))

base_theme <- theme(plot.title = element_text(size = 12, face = "bold"),
                    axis.title = element_text(size = 8), axis.text = element_text(size = 7))

p_anno <- utils_umap_by(epi, "annotation", label = TRUE, repel = TRUE,
                        pt.size = 0.25, label.size = 3) +
    ggtitle("epithelial annotation") +
    coord_equal(xlim = xr, ylim = yr) + NoLegend() + base_theme

# 유전자별 3환자 공통 색스케일
gene_panel <- function(g, s, lim) {
    cells <- colnames(epi)[epi$orig.ident == s]
    FeaturePlot(epi, features = g, cells = cells, order = TRUE, pt.size = 0.3) +
        scale_color_viridis_c(option = "viridis", limits = lim, name = "expr") +
        ggtitle(sprintf("%s — %s", g, s)) +
        coord_equal(xlim = xr, ylim = yr) + base_theme
}

# 행 배치: [annotation | spacer | spacer] 다음 유전자마다 [CRPC1 | CRPC2 | CRPC3]
panels <- list(p_anno, patchwork::plot_spacer(), patchwork::plot_spacer())
for (g in NE_GENES) {
    lim <- c(0, max(expr_all[[g]]))
    panels <- c(panels, lapply(SAMPLES, function(s) gene_panel(g, s, lim)))
}

combined <- wrap_plots(panels, ncol = 3) +
    plot_annotation(
        title    = "Neuroendocrine markers — patient-split (DNPC NE-negativity check)",
        subtitle = "Rows = NE markers (shared color scale across 3 patients per gene), columns = CRPC1/2/3 · top-left = cluster reference",
        caption  = "assay=SCT (data); INSM1 not shown (excluded from SCT, near-zero expression). Not to be used as evidence in the malignancy argument.",
        theme    = theme(plot.title    = element_text(size = 18, face = "bold"),
                         plot.subtitle = element_text(size = 11),
                         plot.caption  = element_text(size = 9, colour = "grey35"))
    )

utils_save_png(combined, OUT_DIR, "FeaturePlot_NEmarkers_byPatient.png",
               w = 15, h = 4.3 * (length(NE_GENES) + 1), dpi = 200, limitsize = FALSE)
