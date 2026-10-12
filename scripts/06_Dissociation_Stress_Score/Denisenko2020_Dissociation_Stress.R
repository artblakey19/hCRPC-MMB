# Denisenko 2020 dissociation stress score (17 genes, AddModuleScore) 확인
# BE 5 top 마커(FOS/JUN/EGR1/IER2/IER5L/DDIT3/GADD45G)가 논문에서 지적한 stress signal과 겹침

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
})
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/06_Dissociation_Stress_Score/Denisenko2020"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Denisenko et al. 2020 Genome Biol (doi:10.1186/s13059-020-02048-6) Methods 의 17 genes
# "Systematic assessment of tissue dissociation and storage biases in single-cell and single-nucleus RNA-seq workflows"
DENISENKO2020 <- c("FOS","FOSB","JUN","JUNB","JUND","ATF3","EGR1",
                   "HSPA1A","HSPA1B","HSP90AB1","HSPA8","HSPB1",
                   "IER3","IER2","BTG1","BTG2","DUSP1")

# Load ----
epi <- readRDS(IN_RDS)
message("Loaded: ", IN_RDS, "  (", ncol(epi), " cells)")

DefaultAssay(epi) <- "RNA"  # counts layer 만 있어 normalize
epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)

# Module score ----
genes <- intersect(DENISENKO2020, rownames(epi[["RNA"]]))
miss  <- setdiff(DENISENKO2020, genes)
message(sprintf("  Denisenko2020 %d/%d genes present%s", length(genes), length(DENISENKO2020),
                if (length(miss)) paste0(" (missing: ", paste(miss, collapse = ","), ")") else ""))

# AddModuleScore
epi <- AddModuleScore(epi, features = list(genes), name = "Denisenko2020_sc", assay = "RNA")
epi$Denisenko2020_sc  <- epi$Denisenko2020_sc1  # AddModuleScore 가 붙인 "1" 제거
epi$Denisenko2020_sc1 <- NULL

# annotation 별 요약 ----
summ <- epi@meta.data %>%
    group_by(annotation) %>%
    summarise(n = n(),
              Denisenko2020_sc_mean   = round(mean(Denisenko2020_sc), 4),
              Denisenko2020_sc_median = round(median(Denisenko2020_sc), 4),
              .groups = "drop") %>%
    arrange(desc(Denisenko2020_sc_mean))
write.csv(summ, file.path(OUT_DIR, "Denisenko2020_per_cluster_scores.csv"), row.names = FALSE)
cat("\n=== Per-cluster stress scores (sorted by Denisenko2020 mean) ===\n")
print(as.data.frame(summ), row.names = FALSE)

# Plots ----
p <- VlnPlot(epi, features = "Denisenko2020_sc", group.by = "annotation",
             cols = utils_cb_cols(epi$annotation), pt.size = 0) +
    ggtitle("Denisenko 2020 dissociation stress score (17 genes)") +
    theme(legend.position = "none", axis.title.x = element_blank())
utils_save_png(p, OUT_DIR, "VlnPlot_Denisenko2020_byCluster.png", w = 11, h = 6, dpi = 200)

p <- FeaturePlot(epi, features = "Denisenko2020_sc", pt.size = 0.3, order = TRUE) +
    scale_color_viridis_c(option = "viridis") +
    ggtitle("Denisenko 2020 dissociation stress score")
utils_save_png(p, OUT_DIR, "UMAP_Denisenko2020_score.png", w = 10, h = 8, dpi = 200)

cat("\nDONE — outputs in ", OUT_DIR, "\n")
