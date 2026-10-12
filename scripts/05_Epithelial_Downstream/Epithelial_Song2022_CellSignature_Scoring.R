# Song 2022 cell-type signature score — Normal vs PCa signature, AddModuleScore
# Song2022_signature_scores.rds 가 있으면 scoring 은 건너뛰고 figures만 재생성 (FORCE_RECOMPUTE 로 재계산)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(msigdbr)
library(patchwork)
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
EPI_RDS      <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
SIG_XLSX     <- "Resources/Song2022_SuppData2.xlsx"
ANNOT_CSV    <- "Results/05_Epithelial_Downstream/Annotation/annotation_per_cell.csv"
ANNOT_LEVELS <- "Results/05_Epithelial_Downstream/Annotation/label_levels.txt"

OUT_DIR <- "Results/05_Epithelial_Downstream/Song2022_signature_scores"
OUT_RDS <- "Results/05_Epithelial_Downstream/Song2022_signature_scores.rds"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FORCE_RECOMPUTE <- FALSE

# Song 2022 Nat Commun Supp Data 2 — sheet 와 상피 signature 열 (기질·면역 열 제외)
SONG_SIGNATURES <- list(
    Normal = list(sheet = "Normal Signature",
                  cols  = c("BE", "Club", "Hillock", "LE", "NE")),
    PCa    = list(sheet = "PCa signature",
                  cols  = c("BE", "Club", "LE", "ERGpos_Tumor", "ERGneg_Tumor"))
)

# Signature score ----
# score 열: <Normal|PCa>_<col>_Score, Hallmark_AR_Score
if (file.exists(OUT_RDS) && !FORCE_RECOMPUTE) {
    message("Loading cached object: ", OUT_RDS)
    epi <- readRDS(OUT_RDS)
} else {
    epi <- readRDS(EPI_RDS)
    DefaultAssay(epi) <- "RNA"  # counts layer 만 있어 normalize
    epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)

    sigs <- list()
    for (src in names(SONG_SIGNATURES)) {
        df <- read.xlsx(SIG_XLSX, sheet = SONG_SIGNATURES[[src]]$sheet)
        for (col in SONG_SIGNATURES[[src]]$cols) {
            sigs[[paste(src, col, sep = "_")]] <- unique(na.omit(df[[col]]))
        }
    }
    # Hallmark AR score. Normal/PCa 소속이 아니라 공통 reference score로 둠.
    hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
    sigs[["Hallmark_AR"]] <- unique(hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_ANDROGEN_RESPONSE"])

    for (s in names(sigs)) {
        genes <- intersect(sigs[[s]], rownames(epi))
        message(sprintf("  %s: %d / %d genes present", s, length(genes), length(sigs[[s]])))
        if (length(genes) < 5) {
            message("Skip ", s, " — fewer than 5 valid genes")
            next
        }
        sc  <- paste0(s, "_Score")
        epi <- AddModuleScore(epi, features = list(genes), name = sc, assay = "RNA")
        epi@meta.data[[sc]] <- epi@meta.data[[paste0(sc, "1")]]  # AddModuleScore 가 붙인 "1" 제거
        epi@meta.data[[paste0(sc, "1")]] <- NULL
    }
    saveRDS(epi, OUT_RDS)  # score 만 (annotation 은 매번 CSV 로)
}

# Annotation (Epithelial_Annotation.R 의 per-cell label) ----
ann <- read.csv(ANNOT_CSV)
epi$annotation <- ann$annotation[match(colnames(epi), ann$cell)]
epi$plot_group <- factor(epi$annotation, levels = readLines(ANNOT_LEVELS))

score_cols <- list(
    Normal   = grep("^Normal_.*_Score$",   colnames(epi@meta.data), value = TRUE),
    PCa      = grep("^PCa_.*_Score$",      colnames(epi@meta.data), value = TRUE),
    Hallmark = grep("^Hallmark_AR_Score$", colnames(epi@meta.data), value = TRUE)
)
for (src in names(score_cols)) {
    message(src, " score columns: ", paste(score_cols[[src]], collapse = ", "))
}

# Per-score plots → {Normal,PCa,Hallmark}/ ----
group_cols  <- utils_cb_cols(epi$plot_group)
x_labels_45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
n_samples   <- length(unique(epi$orig.ident))

# "Normal_BE_Score" → "Normal - BE"
score_title <- function(sc) sub("^([^_]+)_(.*)_Score$", "\\1 - \\2", sc)

plot_score <- function(sc, out_dir) {
    title <- score_title(sc)
    lims  <- range(epi@meta.data[[sc]], na.rm = TRUE)
    df    <- data.frame(group = epi$plot_group, score = epi@meta.data[[sc]])

    p <- FeaturePlot(epi, features = sc, pt.size = 0.3, order = TRUE) +
        scale_color_viridis_c(option = "viridis", limits = lims) +
        ggtitle(title)
    utils_save_png(p, out_dir, paste0("UMAP_", sc, ".png"),
                   w = 20, h = 20, units = "cm", dpi = 200)

    # split 패널 공통 color range
    p <- FeaturePlot(epi, features = sc, split.by = "orig.ident", pt.size = 0.3, order = TRUE) &
        scale_color_viridis_c(option = "viridis", limits = lims)
    p <- p + plot_annotation(title = paste(title, "split by sample", sep = " - "))
    utils_save_png(p, out_dir, paste0("UMAP_", sc, "_splitBySample.png"),
                   w = 8 * n_samples, h = 8, dpi = 200)

    # 파란 막대(shape 95) = group 평균
    p <- ggplot(df, aes(group, score, fill = group)) +
        geom_boxplot(outlier.size = 0.1) +
        scale_fill_manual(values = group_cols) +
        stat_summary(fun = mean, geom = "point", size = 12, colour = "blue", shape = 95) +
        ggtitle(title) + NoLegend() + x_labels_45
    utils_save_png(p, out_dir, paste0("Boxplot_", sc, ".png"),
                   w = 20, h = 20, units = "cm", dpi = 200)

    p <- ggplot(df, aes(group, score, fill = group)) +
        geom_violin(scale = "width") +
        scale_fill_manual(values = group_cols) +
        ggtitle(title) + NoLegend() + x_labels_45
    utils_save_png(p, out_dir, paste0("Violin_", sc, ".png"),
                   w = 20, h = 15, units = "cm", dpi = 200)
}

for (src in names(score_cols)) {
    out <- file.path(OUT_DIR, src)
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    for (sc in score_cols[[src]]) plot_score(sc, out)
}

# Normal vs PCa combined violin ----
long_df <- epi@meta.data %>%
    select(group = plot_group, all_of(c(score_cols$Normal, score_cols$PCa))) %>%
    pivot_longer(-group, names_to = "score_name", values_to = "score") %>%
    mutate(source    = sub("^([^_]+)_.*_Score$", "\\1", score_name),
           signature = sub("^[^_]+_(.*)_Score$", "\\1", score_name))

p <- ggplot(long_df, aes(group, score, fill = group)) +
    geom_violin(scale = "width") +
    facet_grid(source ~ signature, scales = "free_y") +
    scale_fill_manual(values = group_cols) +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text  = element_text(size = 10)) +
    ggtitle("Normal vs PCa signature scores")
utils_save_png(p, OUT_DIR, "Violin_Normal_vs_PCa_all_signatures.png",
               w = 34, h = 18, units = "cm", dpi = 200)

message("Done. Outputs in: ", OUT_DIR)
