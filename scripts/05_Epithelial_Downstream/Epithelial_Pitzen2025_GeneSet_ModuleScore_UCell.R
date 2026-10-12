#!/usr/bin/env Rscript
# Pitzen 2025 gene set score — AddModuleScore / UCell → 상피 UMAP
# sets: Pitzen 2025 PNAS 122:e2415308122 SI Table S1 (BPECT, Henry 2018), S2 (CRPCST), S4 (AR/KLF5 activity, BCH TFs, cell cycle)

suppressPackageStartupMessages({
    library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(UCell)
})
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
N_CORES <- 30  # UCell

EPI_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
GS_TSV  <- "Resources/Pitzen2025/Pitzen2025_gene_sets.tsv"
OUT_DIR <- "Results/05_Epithelial_Downstream/Pitzen2025_GeneSet_Scores"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SET_GROUPS <- list(
    BPECT     = c("Luminal", "Basal", "Club", "Hillock", "NE"),
    OtherSets = c("CRPC_AR", "CRPC_WNT", "CRPC_SCL", "CRPC_NE",
                  "AR_activity", "KLF5_activity", "BCH_TFs", "CellCycle"))
SET_LABEL <- c(Luminal = "Luminal", Basal = "Basal", Club = "Club", Hillock = "Hillock", NE = "Neuroendocrine",
               CRPC_AR = "CRPC-AR", CRPC_WNT = "CRPC-WNT", CRPC_SCL = "CRPC-SCL", CRPC_NE = "CRPC-NE",
               AR_activity = "AR activity", KLF5_activity = "KLF5 activity", BCH_TFs = "BCH TFs",
               CellCycle = "Cell cycle")
GROUP_TITLE <- c(BPECT     = "Pitzen 2025 BPECT gene sets (SI Table S1)",
                 OtherSets = "Pitzen 2025 CRPC subtype / activity gene sets (SI Tables S2, S4)")
METHODS <- c(MS = "AddModuleScore", UC = "UCell")

# Load ----
epi <- readRDS(EPI_RDS)
DefaultAssay(epi) <- "RNA"  # counts layer 만 있어 normalize
epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)

# Gene sets (RNA 에 없는 유전자 제외) ----
gs <- read.delim(GS_TSV, comment.char = "#")
gs <- gs[gs$set %in% unlist(SET_GROUPS), ]
gs$in_RNA <- gs$gene_current %in% rownames(epi[["RNA"]])
write.table(gs, file.path(OUT_DIR, "gene_set_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sets <- lapply(split(gs$gene_current[gs$in_RNA], gs$set[gs$in_RNA]), unique)[unlist(SET_GROUPS)]
coverage <- gs |> group_by(set) |>
    summarise(n_paper = n(), n_used = sum(in_RNA),
              missing = paste(gene_paper[!in_RNA], collapse = ","), .groups = "drop")
cat("\n=== gene set coverage ===\n"); print(as.data.frame(coverage), row.names = FALSE)
stopifnot(all(lengths(sets) >= 5))

# AddModuleScore (RNA data, seed 42) ----
for (s in names(sets)) {
    nm  <- paste0("MS_", s)
    epi <- AddModuleScore(epi, features = list(sets[[s]]), name = nm, assay = "RNA", seed = 42)
    epi@meta.data[[nm]] <- epi@meta.data[[paste0(nm, "1")]]  # AddModuleScore 가 붙인 "1" 제거
    epi@meta.data[[paste0(nm, "1")]] <- NULL
}

# UCell (RNA counts, rank 기반) ----
epi <- AddModuleScore_UCell(epi, features = sets, assay = "RNA", slot = "counts", ncores = N_CORES)
for (s in names(sets)) {
    epi@meta.data[[paste0("UC_", s)]] <- epi@meta.data[[paste0(s, "_UCell")]]
    epi@meta.data[[paste0(s, "_UCell")]] <- NULL
}

# Score tables ----
score_cols <- c(paste0("MS_", names(sets)), paste0("UC_", names(sets)))
write.csv(data.frame(cell = colnames(epi), patient = epi$orig.ident, annotation = epi$annotation,
                     epi@meta.data[, score_cols]),
          file.path(OUT_DIR, "per_cell_scores.csv"), row.names = FALSE)
score_summary <- epi@meta.data |> group_by(plot_group = annotation) |>
    summarise(n = n(), across(all_of(score_cols), ~ round(mean(.x), 4)), .groups = "drop")
write.csv(score_summary, file.path(OUT_DIR, "score_mean_byAnnotation.csv"), row.names = FALSE)

# UMAP: annotation + scores ----
# 패널 공통 좌표 범위 (UMAP range +3%)
emb <- Embeddings(epi, "umap")
xr  <- range(emb[, 1]) + 0.03 * c(-1, 1) * diff(range(emb[, 1]))
yr  <- range(emb[, 2]) + 0.03 * c(-1, 1) * diff(range(emb[, 2]))

base_theme <- theme(plot.title  = element_text(size = 11, face = "bold"),
                    axis.title  = element_text(size = 7),
                    axis.text   = element_text(size = 6),
                    legend.text = element_text(size = 6),
                    legend.key.width = unit(0.25, "cm"))

p_anno <- utils_umap_by(epi, "annotation", label = TRUE, repel = TRUE,
                        pt.size = 0.25, label.size = 3) +
    ggtitle("epithelial annotation") +
    coord_equal(xlim = xr, ylim = yr) + NoLegend() + base_theme

score_umap <- function(col, title) {
    FeaturePlot(epi, features = col, order = TRUE, pt.size = 0.25) +
        scale_color_viridis_c(option = "viridis", name = "score") +
        ggtitle(title) + coord_equal(xlim = xr, ylim = yr) + base_theme
}

for (grp in names(SET_GROUPS)) {
    ss <- SET_GROUPS[[grp]]
    for (m in names(METHODS)) {
        ps <- lapply(ss, function(s) score_umap(paste0(m, "_", s),
                                                sprintf("%s (%d genes)", SET_LABEL[[s]], length(sets[[s]]))))
        comb <- wrap_plots(c(list(p_anno), ps), ncol = 3) +
            plot_annotation(
                title    = sprintf("%s — %s", GROUP_TITLE[[grp]], METHODS[[m]]),
                subtitle = "Pitzen et al. 2025 PNAS 122:e2415308122",
                caption  = if (m == "MS") "AddModuleScore: assay=RNA(data, decontX), ctrl=100, seed=42"
                           else "UCell: assay=RNA(counts, decontX), maxRank=1500",
                theme    = theme(plot.title    = element_text(size = 16, face = "bold"),
                                 plot.subtitle = element_text(size = 9),
                                 plot.caption  = element_text(size = 8, colour = "grey35")))
        utils_save_png(comb, OUT_DIR, sprintf("UMAP_%s_%s.png", grp, METHODS[[m]]),
                       w = 13, h = 4.4 * ceiling((length(ss) + 1) / 3),
                       dpi = 200, limitsize = FALSE)
    }
}

cat("\nDone ->", OUT_DIR, "\n")
