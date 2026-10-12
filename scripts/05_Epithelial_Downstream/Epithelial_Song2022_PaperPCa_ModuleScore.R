# Song 2022 paper PCa gene sets (MSigDB C2:CGP Liu/Tomlins/Wallace/Setlur) — AddModuleScore → UMAP 그림 생성
# epi_song2022_paperpca_modulescore.rds 있으면 그림만 재생성
# 원 논문의 코드인 Seqwell_Combined_E.R를 참고함.

library(Seurat)
library(msigdbr)
library(ggplot2)
library(patchwork)   # plot_annotation / `&` for split.by FeaturePlot panels
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
IN_RDS  <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Song2022_PaperPCa_ModuleScore"
OUT_RDS <- "Results/05_Epithelial_Downstream/epi_song2022_paperpca_modulescore.rds"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Paper named PCa gene sets (Seqwell_Combined_E.R line 896-906에서 사용한 4종)
PAPER_PCA_SETS <- c(
    "LIU_PROSTATE_CANCER_UP",
    "TOMLINS_PROSTATE_CANCER_UP",
    "WALLACE_PROSTATE_CANCER_UP",
    "SETLUR_PROSTATE_CANCER_TMPRSS2_ERG_FUSION_UP"
)

# Module score ----
if (file.exists(OUT_RDS)) {
    message("Loading cached ", OUT_RDS, " — skipping module scoring, regenerating figures")
    epi <- readRDS(OUT_RDS)
} else {
    epi <- readRDS(IN_RDS)
    message("Loaded: ", IN_RDS)

    DefaultAssay(epi) <- "RNA"  # counts layer 만 있어 normalize
    epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)

    c2cgp_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
    pca_sets <- split(c2cgp_df$gene_symbol, c2cgp_df$gs_name)[PAPER_PCA_SETS]

    # Seqwell_Combined_E.R:600에서 사용한 방식 — assay="RNA" + ctrl=default(100).
    for (set_name in PAPER_PCA_SETS) {
        sig <- paste0(set_name, "_sig")
        epi <- AddModuleScore(epi, features = list(pca_sets[[set_name]]), name = sig, assay = "RNA")
        epi@meta.data[[sig]] <- epi@meta.data[[paste0(sig, "1")]]  # AddModuleScore 가 붙인 "1" 제거
        epi@meta.data[[paste0(sig, "1")]] <- NULL
    }
    saveRDS(epi, OUT_RDS)
}

# UMAP: 전체 / 샘플별 ----
n_samples <- length(unique(epi$orig.ident))
for (set_name in PAPER_PCA_SETS) {
    sc_col <- paste0(set_name, "_sig")

    p <- FeaturePlot(epi, features = sc_col, pt.size = 0.3, order = TRUE) +
        scale_color_viridis_c(option = "viridis") + ggtitle(set_name)
    utils_save_png(p, OUT_DIR, paste0("UMAP_", set_name, ".png"), w = 10, h = 8, dpi = 200)

    rng <- range(epi@meta.data[[sc_col]], na.rm = TRUE)
    p <- FeaturePlot(epi, features = sc_col, split.by = "orig.ident",
                     pt.size = 0.3, order = TRUE) &
        scale_color_viridis_c(option = "viridis", limits = rng)
    p <- p + plot_annotation(title = set_name)
    utils_save_png(p, OUT_DIR, paste0("UMAP_", set_name, "_splitBySample.png"),
                   w = 8 * n_samples, h = 8, dpi = 200)
}

message("Done. Outputs in: ", OUT_DIR)
message("Annotated object: ", OUT_RDS)
