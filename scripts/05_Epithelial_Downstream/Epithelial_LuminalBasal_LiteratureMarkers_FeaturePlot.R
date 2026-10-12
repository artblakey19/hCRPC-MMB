# Stage 05 — Song 2022 외 문헌의 luminal / basal / club / hillock 마커 FeaturePlot · DotPlot · AddModuleScore (SCT)
# lineage identity 확인용 — malignancy 판정 근거 아님
# 마커 출처 (human prostate scRNA-seq, 본문/그림에 유전자명 명시)
#   Henry2018  Cell Rep 2018;25:3530 (PMID 30566875) — BE / LE / Club(OE1) / Hillock(OE2) 라벨 원출처
#   Joseph2020 Prostate 2020;80:872 (PMID 32497356) — urethral luminal(club/hillock) progenitor 마커 dot plot
#   Hirz2023   Nat Commun 2023;14:663 (PMID 36750562) — PCa 조직에서 4 아형 key-marker 재현
#   Pitzen2025 PNAS 2025;122:e2415308122 (PMID 39913208) — BPECT set, CRPC 에서 basal/club/hillock 동시 상승 (readout PI3)
#   Canonical  교과서적 basal/luminal 마커 (위 논문 공통 인용)

suppressMessages({
    library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})
source("scripts/00_utils/scRNA_utils.R")

# Settings ----
EPI_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/LuminalBasal_LiteratureMarkers"
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Marker table (cell_type, source, gene) ----
# ACPP → ACP3 (현재 HGNC 심볼)
marker_tbl <- rbind(
    # Basal epithelia
    data.frame(cell_type = "Basal",   source = "Henry2018",
               gene = c("KRT5", "KRT14", "TP63", "NOTCH4", "LTBP2", "DKK1")),
    data.frame(cell_type = "Basal",   source = "Canonical",
               gene = c("KRT15", "DST", "COL17A1", "ITGA6", "ITGB4",
                        "NGFR", "PDPN", "BNC1", "MIR205HG", "CDH3")),

    # Luminal epithelia
    data.frame(cell_type = "Luminal", source = "Henry2018",
               gene = c("KLK3", "ACP3", "MSMB", "GP2", "NEFH", "NPY", "DHRS7")),
    data.frame(cell_type = "Luminal", source = "Canonical",
               gene = c("KRT8", "KRT18", "AR", "NKX3-1", "KLK2", "TMPRSS2",
                        "FOLH1", "TGM4", "DPP4", "CD38", "PLA2G2A", "RDH11")),

    # Club (Henry OE1)
    data.frame(cell_type = "Club",    source = "Henry2018",
               gene = c("SCGB1A1", "SCGB3A1", "LY6D", "PSCA")),
    data.frame(cell_type = "Club",    source = "Joseph2020",
               gene = c("TACSTD2", "KRT4", "PIGR", "MMP7", "LTF", "CP")),
    data.frame(cell_type = "Club",    source = "Hirz2023",
               gene = c("WFDC2", "LCN2", "OLFM4", "CRABP2")),
    data.frame(cell_type = "Club",    source = "Pitzen2025",
               gene = "PI3"),

    # Hillock (Henry OE2)
    data.frame(cell_type = "Hillock", source = "Henry2018",
               gene = c("KRT13", "AKR1C1", "AKR1C2", "KRT19")),
    data.frame(cell_type = "Hillock", source = "Joseph2020",
               gene = c("KRT7", "RARRES1", "LYPD3", "AQP3")),
    data.frame(cell_type = "Hillock", source = "Hirz2023",
               gene = c("S100A2", "SERPINB3", "SERPINB4", "APOBEC3A", "CLDN4"))
)

CELL_TYPES <- c("Basal", "Luminal", "Club", "Hillock")

# core 마커 (환자별 split 용)
CORE_GENES <- c("KRT5", "TP63", "KRT14",          # basal
                "KLK3", "NKX3-1", "AR",           # luminal
                "SCGB1A1", "PIGR",                # club
                "KRT13", "KRT19")                 # hillock

# Load ----
epi <- readRDS(EPI_RDS)
DefaultAssay(epi) <- "SCT"

# 객체 내 marker 유무 기록
marker_tbl$in_SCT <- marker_tbl$gene %in% rownames(epi[["SCT"]])
marker_tbl$in_RNA <- marker_tbl$gene %in% rownames(epi[["RNA"]])
write.table(marker_tbl, file.path(OUT_DIR, "marker_gene_status.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 공통 UMAP 범위 / 테마 / annotation 참조 패널 ----
emb <- Embeddings(epi, "umap")
xr  <- range(emb[, 1]) + 0.03 * c(-1, 1) * diff(range(emb[, 1]))
yr  <- range(emb[, 2]) + 0.03 * c(-1, 1) * diff(range(emb[, 2]))

base_theme <- theme(plot.title  = element_text(size = 11, face = "bold"),
                    axis.title  = element_text(size = 7),
                    axis.text   = element_text(size = 6),
                    legend.text = element_text(size = 6),
                    legend.key.width = unit(0.25, "cm"))

panel_titles <- function(title, subtitle, caption, title_size = 16, subtitle_size = 9) {
    plot_annotation(
        title = title, subtitle = subtitle, caption = caption,
        theme = theme(plot.title    = element_text(size = title_size, face = "bold"),
                      plot.subtitle = element_text(size = subtitle_size),
                      plot.caption  = element_text(size = 8, colour = "grey35"))
    )
}

p_anno <- utils_umap_by(epi, "annotation", label = TRUE, repel = TRUE,
                        pt.size = 0.25, label.size = 3) +
    ggtitle("epithelial annotation") +
    coord_equal(xlim = xr, ylim = yr) + NoLegend() + base_theme
utils_save_png(p_anno, OUT_DIR, "UMAP_annotation_reference.png",
               w = 16, h = 16, units = "cm", dpi = 200)

gene_umap <- function(g, src) {
    FeaturePlot(epi, features = g, order = TRUE, pt.size = 0.25) +
        scale_color_viridis_c(option = "viridis", name = "expr") +
        ggtitle(sprintf("%s  [%s]", g, src)) +
        coord_equal(xlim = xr, ylim = yr) + base_theme
}

# cell type 별 FeaturePlot grid ----
for (ct in CELL_TYPES) {
    ct_tbl <- marker_tbl[marker_tbl$cell_type == ct, ]
    panels <- c(list(p_anno), Map(gene_umap, ct_tbl$gene, ct_tbl$source))
    n_col  <- 4
    n_row  <- ceiling(length(panels) / n_col)

    comb <- wrap_plots(panels, ncol = n_col) +
        panel_titles(
            title    = sprintf("%s markers — literature panel (excluding Song 2022)", ct),
            subtitle = paste0("Top-left = annotation reference · [source] shown per panel · ",
                              "Sources: Henry 2018 Cell Rep / Joseph 2020 Prostate / ",
                              "Hirz 2023 Nat Commun / Pitzen 2025 PNAS"),
            caption  = "assay=SCT(data), color=viridis. For lineage identity check — not evidence for malignancy."
        )
    utils_save_png(comb, OUT_DIR, sprintf("FeaturePlot_%s.png", ct),
                   w = 4.2 * n_col, h = 4.4 * n_row, dpi = 200, limitsize = FALSE)
}

# core 마커 patient-split FeaturePlot ----
# 유전자별 3환자 공통 색스케일
expr_core <- FetchData(epi, vars = CORE_GENES, layer = "data")

gene_panel_patient <- function(g, s, lim) {
    cells <- colnames(epi)[epi$orig.ident == s]
    FeaturePlot(epi, features = g, cells = cells, order = TRUE, pt.size = 0.3) +
        scale_color_viridis_c(option = "viridis", limits = lim, name = "expr") +
        ggtitle(sprintf("%s — %s", g, s)) +
        coord_equal(xlim = xr, ylim = yr) + base_theme
}

# 행 배치: [annotation | spacer | spacer] 다음 유전자마다 [CRPC1 | CRPC2 | CRPC3]
panels <- list(p_anno, patchwork::plot_spacer(), patchwork::plot_spacer())
for (g in CORE_GENES) {
    lim <- c(0, max(expr_core[[g]]))
    panels <- c(panels, lapply(SAMPLES, function(s) gene_panel_patient(g, s, lim)))
}

comb <- wrap_plots(panels, ncol = 3) +
    panel_titles(
        title    = "Core luminal / basal / club / hillock markers — patient-split",
        subtitle = "Rows = markers (shared color scale across 3 patients per gene), columns = CRPC1/2/3 · top-left = annotation reference",
        caption  = "assay=SCT(data). Sources: Henry 2018 / Joseph 2020 / Hirz 2023 / Pitzen 2025.",
        title_size = 17, subtitle_size = 10
    )
utils_save_png(comb, OUT_DIR, "FeaturePlot_CoreMarkers_byPatient.png",
               w = 14, h = 4.4 * (length(CORE_GENES) + 1), dpi = 200, limitsize = FALSE)

# DotPlot: 전체 마커 × annotation ----
dot_genes <- marker_tbl$gene

p_dot <- DotPlot(epi, features = dot_genes, group.by = "annotation",
                 assay = "SCT", dot.scale = 5) +
    scale_color_viridis_c(option = "viridis", name = "avg expr\n(scaled)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = element_text(size = 9)) +
    labs(x = NULL, y = NULL,
         title = "Literature luminal/basal/club/hillock markers by epithelial cluster",
         subtitle = paste("Gene order =", paste(CELL_TYPES, collapse = " → "),
                          "· Sources: Henry 2018 / Joseph 2020 / Hirz 2023 / Pitzen 2025"))

utils_save_png(p_dot, OUT_DIR, "DotPlot_AllMarkers_byCluster.png",
               w = 0.28 * length(dot_genes) + 6, h = 10, dpi = 200, limitsize = FALSE)

# cell type 별 AddModuleScore ----
# 개별 유전자 dropout 에 덜 흔들리는 요약
score_cols   <- paste0("Lit_", CELL_TYPES, "_Score")
score_titles <- paste(CELL_TYPES, "signature score")

for (i in seq_along(CELL_TYPES)) {
    nm    <- score_cols[i]
    genes <- marker_tbl$gene[marker_tbl$cell_type == CELL_TYPES[i]]
    epi <- AddModuleScore(epi, features = list(genes), name = nm, assay = "SCT", ctrl = 50)
    # 접미사 "1" 제거
    epi@meta.data[[nm]] <- epi@meta.data[[paste0(nm, "1")]]
    epi@meta.data[[paste0(nm, "1")]] <- NULL
}

score_umaps <- Map(function(sc, title) {
    FeaturePlot(epi, features = sc, order = TRUE, pt.size = 0.25) +
        scale_color_viridis_c(option = "viridis", name = "score") +
        ggtitle(title) +
        coord_equal(xlim = xr, ylim = yr) + base_theme
}, score_cols, score_titles)

comb <- wrap_plots(c(list(p_anno), score_umaps), ncol = 3) +
    panel_titles(
        title    = "Literature cell-type signature scores (AddModuleScore)",
        subtitle = "Union of Henry 2018 / Joseph 2020 / Hirz 2023 / Pitzen 2025 markers, per cell type",
        caption  = "assay=SCT(data), ctrl=50. Summary score less sensitive to single-gene dropout."
    )
utils_save_png(comb, OUT_DIR, "FeaturePlot_SignatureScores.png",
               w = 13, h = 4.4 * ceiling((length(score_cols) + 1) / 3),
               dpi = 200, limitsize = FALSE)

vln <- Map(function(sc, title) {
    df <- data.frame(group = epi$annotation, score = epi@meta.data[[sc]])
    ggplot(df, aes(group, score, fill = group)) +
        geom_violin(scale = "width") +
        geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.6) +
        scale_fill_manual(values = utils_cb_cols(df$group)) +
        ggtitle(title) +
        labs(x = NULL, y = "module score") + NoLegend() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title  = element_text(size = 12, face = "bold"))
}, score_cols, score_titles)
utils_save_png(wrap_plots(vln, ncol = 2), OUT_DIR, "Violin_SignatureScores_byCluster.png",
               w = 26, h = 9 * ceiling(length(vln) / 2), units = "cm",
               dpi = 200, limitsize = FALSE)

# cluster 별 평균 점수 표
score_summary <- epi@meta.data %>%
    group_by(plot_group = annotation) %>%
    summarise(n = n(), across(all_of(score_cols), ~ round(mean(.x), 4)), .groups = "drop")

write.csv(score_summary, file.path(OUT_DIR, "signature_score_mean_byCluster.csv"),
          row.names = FALSE)
print(as.data.frame(score_summary))

message("Done. Outputs in: ", OUT_DIR)
