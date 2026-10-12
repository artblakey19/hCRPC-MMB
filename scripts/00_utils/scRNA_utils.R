# scRNA_utils.R — 단계별 스크립트가 공통으로 source 하는 함수

library(Seurat)

# Colorblind-friendly palette ----
#' Okabe-Ito palette, n > 8 이면 interpolate
utils_cb_palette <- function(n) {
  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999"
  )
  if (n <= length(okabe_ito)) {
    okabe_ito[seq_len(n)]
  } else {
    grDevices::colorRampPalette(okabe_ito)(n)
  }
}

#' x 의 그룹 수만큼 utils_cb_palette (NA 제외)
utils_cb_cols <- function(x) utils_cb_palette(length(unique(na.omit(x))))

# Plot helpers ----
#' group 열로 칠한 UMAP, ... → DimPlot
utils_umap_by <- function(obj, group, ...) {
  DimPlot(obj, reduction = "umap", group.by = group,
          cols = utils_cb_cols(obj@meta.data[[group]]), ...)
}

#' 흰 배경 ggsave, ... → ggsave (dpi 등)
utils_save_png <- function(plot, dir, file, w, h, ...) {
  ggplot2::ggsave(file.path(dir, file), plot = plot, width = w, height = h,
                  bg = "white", ...)
}

# Tables ----
utils_write_crosstab <- function(rows, cols, dir, file) {
  write.csv(as.data.frame.matrix(table(rows, cols)), file.path(dir, file))
}

# Clustering ----
# SCTransform 결과(cluster 구성에도 영향)는 future worker 수에 따라 달라짐 (seed 와 무관).
# sctransform 이 유전자를 worker 수만큼 나눠 glmGamPoi로 연산하는데 같이 계산하는 유전자 구성에 따라 결과가 약간씩 다름.
# 저장된 결과와 label map 은 30 workers 기준 → N_CORES 와 무관하게 SCT 는 SCT_WORKERS 로 실행.
# RStudio, Windows는 fork 가 안 돼 worker를 하나만 사용해서 결과를 재현할 수 없음, worker 수가 맞지 않으면 raise error.
SCT_WORKERS <- 30

#' SCTransform → PCA → Harmony → JoinLayers (RNA 는 환자별 split 상태여야 함)
utils_sct_harmony <- function(obj) {
  with(future::plan("multicore", workers = SCT_WORKERS), local = TRUE)
  # worker 수가 맞지 않으면 error
  stopifnot(future::nbrOfWorkers() == SCT_WORKERS)
  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  set.seed(42)  # Harmony k-means init RNG
  obj <- IntegrateLayers(
    object = obj,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    normalization.method = "SCT",
    verbose = FALSE
  )
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj
}

# Markers ----
#' FindAllMarkers (only.pos, min.pct 0.25) → CSV
utils_save_all_markers <- function(
  seurat_obj, output_csv, logfc_threshold = 0.5
) {
  if (DefaultAssay(seurat_obj) == "SCT") {
    seurat_obj <- PrepSCTFindMarkers(seurat_obj)
  }

  all_markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25, logfc.threshold = logfc_threshold
  )
  write.csv(all_markers, file = output_csv, row.names = FALSE)
  return(all_markers)
}
