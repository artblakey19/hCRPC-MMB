# Stage 05 — notable genes (lineage plasticity / NEPC / tumor suppressor / 대사) FeaturePlot + violin, 환자별 split 포함

suppressMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
})
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
IN_RDS   <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR  <- "Results/05_Epithelial_Downstream/Annotation"
PATIENTS <- c("CRPC1", "CRPC2", "CRPC3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 관련 연구에서 많이 언급되는 genes로 구성
GENES <- c(
    "CHD7", "MYC", "KMT2C", "KRT7", "SOX2",
    "FOXA1",      # ARPC pioneer TF
    "SYP", "AR", "TP53", "RB1", "PTEN", "REST", "YAP1",
    "TP63",       # dNp63 isoform 은 구분 불가, 유전자 수준
    "AMACR",      # PCa marker
    "KRT13",      # Hillock
    "SERPINB1",   # Hillock
    "KRT15",      # basal epithelial marker
    "KRT5",       # basal epithelial marker
    "GGCT",       # gamma-glutamyl cyclotransferase
    "PRPS1",      # phosphoribosyl pyrophosphate synthetase 1
    "PNP",        # purine nucleoside phosphorylase
    "CCL2",       # C-C motif chemokine ligand 2 / MCP-1
    "KLF5"        # Kruppel-like factor 5; basal/squamous lineage TF
)

epi <- readRDS(IN_RDS)
DefaultAssay(epi) <- "SCT"

# 전부 HGNC 심볼 그대로 (requested = resolved)
status <- data.frame(requested = GENES, resolved = GENES,
                     in_object = GENES %in% rownames(epi))
write.table(status, file.path(OUT_DIR, "NotableGenes_gene_status.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 전체 상피 ----
fp_list <- lapply(GENES, function(g) {
    FeaturePlot(epi, features = g, order = TRUE, pt.size = 0.25) +
        ggtitle(g) +
        theme(plot.title = element_text(face = "bold"))
})
utils_save_png(wrap_plots(fp_list, ncol = 4), OUT_DIR, "FeaturePlot_NotableGenes.png",
               w = 20, h = 18, dpi = 200)

vp <- VlnPlot(epi, features = GENES, group.by = "annotation",
              pt.size = 0, stack = TRUE, flip = TRUE) +
    NoLegend() +
    ggtitle("Notable genes by annotation")
utils_save_png(vp, OUT_DIR, "VlnPlot_NotableGenes.png", w = 10, h = 11, dpi = 200)

# 환자별 split ----
# 색 상한: 유전자별 전체 상피 q99(발현 세포 기준)로 고정
expr_mat <- GetAssayData(epi, assay = "SCT", layer = "data")[GENES, ]
gene_q99 <- apply(expr_mat, 1, function(v) quantile(v[v > 0], 0.99))

fp_split_rows <- lapply(GENES, function(g) {
    row_plots <- lapply(PATIENTS, function(p) {
        cells_p <- colnames(epi)[epi$orig.ident == p]
        FeaturePlot(epi, features = g, cells = cells_p,
                    order = TRUE, pt.size = 0.25,
                    min.cutoff = 0, max.cutoff = gene_q99[[g]]) +
            ggtitle(paste0(g, " | ", p)) +
            theme(plot.title = element_text(face = "bold", size = 10),
                  legend.position = "right")
    })
    wrap_plots(row_plots, nrow = 1)
})
utils_save_png(wrap_plots(fp_split_rows, ncol = 1), OUT_DIR,
               "FeaturePlot_NotableGenes_byPatient.png",
               w = 5 * length(PATIENTS), h = 4 * length(GENES),
               dpi = 180, limitsize = FALSE)

vp_split <- VlnPlot(epi, features = GENES, group.by = "annotation",
                    split.by = "orig.ident",
                    pt.size = 0, stack = TRUE, flip = TRUE) +
    ggtitle("Notable genes by annotation, split by patient")
utils_save_png(vp_split, OUT_DIR, "VlnPlot_NotableGenes_byPatient.png",
               w = 12, h = 12, dpi = 200)

message("Notable genes figures → ", OUT_DIR)
