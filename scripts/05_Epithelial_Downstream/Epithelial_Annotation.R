# Stage 05 — filtered epi res-0.3 cluster(12개) → cell-type label
# epi_annotated.rds, label_levels.txt(표시 순서) → 이후 스크립트 input

suppressMessages({
    library(Seurat)
    library(ggplot2)
})
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
IN_RDS  <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Annotation"
OUT_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Cluster → label map ----
#   BE 5 = cl3 (AP-1/heat-stress), OE = cl6 (CFTR+), LE(ARPC) = cl9 (KLK2/FOLH1),
#   Ionocyte = cl11 (FOXI1/ATP6V), Club = cl2 (MMP7/PIGR).
cluster_to_label <- c(
    "0"  = "BE 1",
    "1"  = "BE 4",
    "2"  = "Club",
    "3"  = "BE 5",
    "4"  = "Hillock 1",
    "5"  = "BE 2",
    "6"  = "OE",
    "7"  = "Hillock 2",
    "8"  = "BE 6",
    "9"  = "LE(ARPC)",
    "10" = "BE 3",
    "11" = "Ionocyte"
)

# 표시 순서
label_levels <- c(
    "LE(ARPC)", "Club", "Hillock 1", "Hillock 2",
    "BE 1", "BE 2", "BE 3",
    "BE 4", "BE 5", "BE 6",
    "OE",
    "Ionocyte"
)

# Annotate ----
epi <- readRDS(IN_RDS)
epi$annotation <- factor(unname(cluster_to_label[as.character(epi$seurat_clusters)]),
                         levels = label_levels)
Idents(epi) <- "annotation"

# Figures & tables ----
p <- utils_umap_by(epi, "annotation", label = TRUE, repel = TRUE, pt.size = 0.3) +
    ggtitle("Filtered epithelial — annotation")
utils_save_png(p, OUT_DIR, "UMAP_annotation.png", w = 10, h = 8)

p <- utils_umap_by(epi, "annotation", split.by = "orig.ident",
                   label = TRUE, repel = TRUE, pt.size = 0.3)
utils_save_png(p, OUT_DIR, "UMAP_annotation_split.png", w = 24, h = 8)

tab <- as.data.frame(table(cluster = epi$seurat_clusters,
                           annotation = epi$annotation))
tab <- tab[tab$Freq > 0, ]
write.csv(tab, file.path(OUT_DIR, "cluster_to_label.csv"), row.names = FALSE)

write.csv(
    data.frame(
        cell = colnames(epi),
        seurat_cluster = as.character(epi$seurat_clusters),
        annotation = as.character(epi$annotation)
    ),
    file.path(OUT_DIR, "annotation_per_cell.csv"),
    row.names = FALSE
)

writeLines(label_levels, file.path(OUT_DIR, "label_levels.txt"))

saveRDS(epi, OUT_RDS)
message("Annotated object: ", OUT_RDS)
message("Figures + mapping: ", OUT_DIR)
