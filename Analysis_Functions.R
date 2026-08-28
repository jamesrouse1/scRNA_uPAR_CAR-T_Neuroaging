#!/usr/bin/env Rscript

# Reproducible single-cell analyses used in the updated manuscript.
# Mouse datasets are rebuilt from filtered 10x matrices, filtered with the
# dataset-specific thresholds in Methods, and integrated with SCT or RPCA.
# Published human annotations and embeddings are retained. Observed RNA is
# LogNormalized before module scoring and differential expression. ALRA is
# restricted to PLAUR expression panels. Module scores are calculated before
# cell-type subsetting, and Figures 2D and 2G use signature-gene averages.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(data.table)
  library(readxl)
  library(openxlsx)
})

set.seed(1)
options(future.globals.maxSize = 16 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
if (exists("mem.maxVSize")) mem.maxVSize(130000)

script_arg <- commandArgs()[grepl("^--file=", commandArgs())]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else normalizePath(".")
default_project_root <- dirname(dirname(normalizePath(script_file, mustWork = FALSE)))
project_root <- normalizePath(Sys.getenv("ALEX_MANUSCRIPT_PROJECT_ROOT", default_project_root), mustWork = FALSE)
alex_root <- normalizePath(Sys.getenv("ALEX_ROOT", dirname(project_root)), mustWork = FALSE)
output_root <- file.path(project_root, "outputs")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args)) tolower(args[[1]]) else "all"
allowed <- c("all", "sorted_brain", "post_treatment_brain", "bone_marrow",
             "human_brain_atlas", "freshmg", "psychad", "ainciburu", "validate")
if (!dataset %in% allowed) {
  stop("Dataset must be one of: ", paste(allowed, collapse = ", "))
}

required_packages <- c("Seurat", "Matrix", "dplyr", "tidyr", "ggplot2",
                       "readr", "data.table", "readxl", "openxlsx", "msigdbr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R package(s): ", paste(missing_packages, collapse = ", "))

save_pdf <- function(plot, path, width = 8, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = width, height = height, family = "Helvetica",
                 useDingbats = FALSE, onefile = FALSE, compress = TRUE,
                 bg = "white", version = "1.4")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}

read_10x_any <- function(path) {
  x <- if (dir.exists(path)) Read10X(path, gene.column = 2, unique.features = TRUE) else Read10X_h5(path)
  if (is.list(x)) {
    if (!"Gene Expression" %in% names(x)) stop("Gene Expression matrix not found in ", path)
    x <- x[["Gene Expression"]]
  }
  x
}

make_objects <- function(sample_table) {
  objects <- lapply(seq_len(nrow(sample_table)), function(i) {
    sample <- sample_table$sample[[i]]
    x <- CreateSeuratObject(read_10x_any(sample_table$path[[i]]), project = sample)
    x$sample <- sample
    for (field in setdiff(names(sample_table), c("sample", "path"))) {
      x[[field]] <- sample_table[[field]][[i]]
    }
    x
  })
  names(objects) <- sample_table$sample
  merge(objects[[1]], y = objects[-1], add.cell.ids = sample_table$sample)
}

qc_filter <- function(obj, mt_max, count_min, count_max, feature_min, feature_max, outdir) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  before <- as.data.frame(table(obj$sample), stringsAsFactors = FALSE)
  names(before) <- c("sample", "cells_before")
  p_before <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      group.by = "sample", ncol = 1, pt.size = 0, layer = "counts")
  save_pdf(p_before, file.path(outdir, "01_qc_before_filtering.pdf"), 13, 12)
  keep <- obj$percent.mt < mt_max & obj$nCount_RNA > count_min &
    obj$nCount_RNA < count_max & obj$nFeature_RNA > feature_min &
    obj$nFeature_RNA < feature_max
  obj <- subset(obj, cells = colnames(obj)[keep])
  after <- as.data.frame(table(obj$sample), stringsAsFactors = FALSE)
  names(after) <- c("sample", "cells_after")
  summary <- full_join(before, after, by = "sample") %>%
    mutate(across(starts_with("cells_"), ~replace_na(.x, 0L)),
           cells_removed = cells_before - cells_after,
           percent_retained = 100 * cells_after / cells_before)
  write_csv(summary, file.path(outdir, "qc_cell_counts.csv"))
  p_after <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     group.by = "sample", ncol = 1, pt.size = 0, layer = "counts")
  save_pdf(p_after, file.path(outdir, "02_qc_after_filtering.pdf"), 13, 12)
  obj
}

integrate_sorted_brain <- function(obj) {
  pieces <- SplitObject(obj, split.by = "sample")
  pieces <- lapply(pieces, SCTransform, verbose = FALSE)
  features <- SelectIntegrationFeatures(pieces, nfeatures = 2000)
  pieces <- PrepSCTIntegration(pieces, anchor.features = features, verbose = FALSE)
  anchors <- FindIntegrationAnchors(pieces, normalization.method = "SCT",
                                    anchor.features = features, verbose = FALSE)
  obj <- IntegrateData(anchors, normalization.method = "SCT", verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 1, verbose = FALSE)
  RunUMAP(obj, dims = 1:5, n.neighbors = 30, min.dist = 0.05,
          seed.use = 1, verbose = FALSE)
}

integrate_rpca <- function(obj, dims, neighbors, min_dist, resolution = 0.5) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- IntegrateLayers(obj, method = RPCAIntegration,
                         normalization.method = "SCT", verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims, reduction = "integrated.dr", verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  RunUMAP(obj, dims = dims, reduction = "integrated.dr", n.neighbors = neighbors,
          min.dist = min_dist, seed.use = 1, verbose = FALSE)
}

attach_labels <- function(obj, file, label_columns, target_columns = label_columns) {
  labels <- fread(file, data.table = FALSE)
  idx <- match(colnames(obj), labels$cell_id)
  if (anyNA(idx)) stop(sum(is.na(idx)), " cells are absent from ", file)
  for (i in seq_along(label_columns)) {
    obj[[target_columns[[i]]]] <- labels[[label_columns[[i]]]][idx]
  }
  obj
}

attach_labels_from_rds <- function(obj, file, source_columns, target_columns = source_columns) {
  annotated <- readRDS(file)
  idx <- match(colnames(obj), colnames(annotated))
  if (anyNA(idx)) stop(sum(is.na(idx)), " cells are absent from ", file)
  for (i in seq_along(source_columns)) {
    obj[[target_columns[[i]]]] <- annotated[[source_columns[[i]], drop = TRUE]][idx]
  }
  obj
}

join_and_normalize <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                verbose = FALSE)
}

read_signature <- function(file, species = c("mouse", "human")) {
  species <- match.arg(species)
  genes <- trimws(readLines(file, warn = FALSE))
  unique(genes[nzchar(genes)])
}

signature_files <- list(
  myeloid_mouse = file.path(alex_root, "manuscript_signature_validation_20260820/tables/new_signature_71_mouse.txt"),
  myeloid_human = file.path(alex_root, "manuscript_signature_validation_20260820/tables/new_signature_human_orthologs.txt")
)

read_bm_signature <- function(column) {
  x <- read_excel(file.path(alex_root, "BM signatures.xlsx"), col_names = FALSE)
  unique(na.omit(trimws(as.character(x[[column]][-(1:2)]))))
}

signatures_mouse <- list(
  Myeloid_Inflammatory = read_signature(signature_files$myeloid_mouse, "mouse"),
  Hallmark_Inflammatory_Response = msigdbr::msigdbr(species = "Mus musculus", collection = "H") %>%
    filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% pull(gene_symbol) %>% unique(),
  Chambers_Inflammaging = read_bm_signature(4),
  Kovtonyuk_Myeloid_Bias = read_bm_signature(7)
)

score_object <- function(obj, signatures, prefix = "") {
  for (nm in names(signatures)) {
    genes <- intersect(signatures[[nm]], rownames(obj))
    if (length(genes) < 3) stop("Too few genes found for ", nm)
    set.seed(1)
    obj <- AddModuleScore(obj, features = list(genes), name = paste0(prefix, nm),
                          assay = "RNA", nbin = 24, ctrl = 100, seed = 1,
                          search = FALSE, slot = "data")
    avg_name <- paste0(prefix, nm, "_Average")
    obj[[avg_name]] <- Matrix::colMeans(
      GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
    )
  }
  obj
}

split_violin_geom <- ggproto(
  "GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x))
    group <- data[1, "group"]
    data2 <- transform(data, x = if (group %% 2 == 1) xminv else xmaxv)
    data2 <- data2[order(if (group %% 2 == 1) data2$y else -data2$y), ]
    data2 <- rbind(data2[1, ], data2, data2[nrow(data2), ], data2[1, ])
    data2[c(1, nrow(data2) - 1, nrow(data2)), "x"] <- round(data2[1, "x"])
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(data2, ...))
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ..., trim = TRUE,
                              scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = split_violin_geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, na.rm = na.rm, ...))
}

format_p <- function(p) {
  if (!is.finite(p)) "NA" else format.pval(p, digits = 3, eps = 1e-300)
}

make_split_violin <- function(df, value, group, levels, title, ylab, path,
                              colors = c("#7CA1CC", "#FF4902")) {
  d <- df %>% filter(.data[[group]] %in% levels, is.finite(.data[[value]]))
  d$plot_group <- factor(d[[group]], levels = levels)
  d$comparison <- factor(paste(levels, collapse = " vs "))
  d$side <- ifelse(d$plot_group == levels[[1]], 1L, 2L)
  pvalue <- wilcox.test(d[[value]] ~ d$plot_group, exact = FALSE)$p.value
  ymax <- max(d[[value]], na.rm = TRUE)
  yrange <- diff(range(d[[value]], na.rm = TRUE))
  p <- ggplot(d, aes(comparison, .data[[value]], fill = plot_group, group = side)) +
    geom_split_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.25) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35, color = "black") +
    annotate("text", x = 1, y = ymax + max(0.08 * yrange, 0.03),
             label = paste0("p = ", format_p(pvalue)), size = 3.4) +
    scale_fill_manual(values = setNames(colors, levels)) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.17))) +
    theme_classic(base_size = 11) +
    labs(title = title, x = NULL, y = ylab, fill = NULL)
  save_pdf(p, path, 7.5, 5.5)
  data.frame(group_1 = levels[[1]], n_1 = sum(d$plot_group == levels[[1]]),
             group_2 = levels[[2]], n_2 = sum(d$plot_group == levels[[2]]),
             wilcox_p = pvalue)
}

run_deg <- function(obj, group_col, ident1, ident2, outdir, prefix) {
  Idents(obj) <- factor(obj[[group_col, drop = TRUE]])
  deg <- FindMarkers(obj, ident.1 = ident1, ident.2 = ident2, assay = "RNA",
                     test.use = "wilcox", min.pct = 0.10, logfc.threshold = 0,
                     only.pos = FALSE)
  deg$gene <- rownames(deg)
  fc <- intersect(c("avg_log2FC", "avg_logFC"), names(deg))[[1]]
  keep <- deg$p_val < 0.05 & abs(deg[[fc]]) > 0.5
  filtered <- deg[keep, , drop = FALSE]
  filtered$direction <- ifelse(filtered[[fc]] > 0, "UP", "DOWN")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write_csv(deg, file.path(outdir, paste0(prefix, "_all_genes.csv")))
  write_csv(filtered, file.path(outdir, paste0(prefix, "_filtered_for_enrichr.csv")))
  filtered
}

run_enrichr <- function(deg, outdir, prefix) {
  if (!requireNamespace("enrichR", quietly = TRUE)) stop("The enrichR package is required")
  databases <- c("GO_Biological_Process_2025", "MSigDB_Hallmark_2020", "Reactome_2022")
  results <- list()
  for (direction in c("UP", "DOWN")) {
    genes <- unique(deg$gene[deg$direction == direction])
    if (!length(genes)) next
    er <- enrichR::enrichr(genes, databases)
    for (db in names(er)) results[[paste(direction, db, sep = "_")]] <- er[[db]]
  }
  if (length(results)) write.xlsx(results, file.path(outdir, paste0(prefix, "_enrichr.xlsx")), overwrite = TRUE)
  invisible(results)
}

plot_cell_umap <- function(obj, field, path, title) {
  p <- DimPlot(obj, reduction = "umap", group.by = field, label = TRUE,
               repel = TRUE, raster = FALSE, pt.size = 0.12) +
    ggtitle(title) + theme_classic(base_size = 11)
  save_pdf(p, path, 11, 8)
}

run_sorted_brain <- function() {
  outdir <- file.path(output_root, "01_sorted_brain")
  samples <- data.frame(
    sample = c("F_neg", "F_plus", "M_neg", "M_plus"),
    path = file.path(alex_root, "dataset1", c(
      "Amor_CAV17_Female_UPARneg/filtered_feature_bc_matrix",
      "Amor_CAV17_Female_UPARpos/filtered_feature_bc_matrix",
      "Amor_CAV16_Male_UPARneg/filtered_feature_bc_matrix",
      "Amor_CAV16_Male_UPAR_pos/filtered_feature_bc_matrix")),
    treatment = c("uPAR-", "uPAR+", "uPAR-", "uPAR+"),
    sex = c("Female", "Female", "Male", "Male")
  )
  obj <- qc_filter(make_objects(samples), 15, 500, 20000, 250, 5000, outdir)
  obj <- integrate_sorted_brain(obj)
  obj <- attach_labels_from_rds(
    obj, file.path(alex_root, "dataset1/seurat_objs/annotated_cell_types_d1.RDS"),
    c("cell_type_grouped", "cell_type_grouped"), c("cell_type", "cell_type_grouped")
  )
  obj <- score_object(join_and_normalize(obj), signatures_mouse)
  plot_cell_umap(obj, "cell_type_grouped", file.path(outdir, "cell_type_umap.pdf"), "Sorted aged brain")
  bam <- obj@meta.data %>% mutate(cell = rownames(.), treatment = factor(treatment, c("uPAR-", "uPAR+"))) %>%
    filter(grepl("BAM", cell_type, ignore.case = TRUE))
  stats <- make_split_violin(
    bam, "Myeloid_Inflammatory_Average", "treatment", c("uPAR-", "uPAR+"),
    "Myeloid inflammatory signature in BAMs", "Mean LogNormalized expression",
    file.path(outdir, "Figure_2D_myeloid_inflammatory_BAM.pdf")
  )
  write_csv(stats, file.path(outdir, "Figure_2D_statistics.csv"))
  deg <- run_deg(obj, "treatment", "uPAR+", "uPAR-", file.path(outdir, "DEG"), "uPARplus_vs_uPARminus")
  run_enrichr(deg, file.path(outdir, "DEG"), "uPARplus_vs_uPARminus")
  saveRDS(obj, file.path(outdir, "sorted_brain_processed.rds"), compress = FALSE)
}

run_post_treatment_brain <- function() {
  outdir <- file.path(output_root, "02_post_treatment_brain")
  ids <- c("Y_uPAR_F", "Y_UT_F", "O_UT_F", "Y_uPAR_M", "Y_UT_M", "O_uPAR_F", "O_uPAR_M", "O_UT_M")
  samples <- data.frame(
    sample = ids,
    path = file.path(alex_root, "dataset2", ids, "filtered_feature_bc_matrix.h5"),
    age = ifelse(grepl("^Y_", ids), "Young", "Old"),
    treatment = ifelse(grepl("uPAR", ids), "uPAR", "UT"),
    sex = ifelse(grepl("_F$", ids), "Female", "Male")
  )
  obj <- qc_filter(make_objects(samples), 5, 500, 20000, 250, 5000, outdir)
  obj <- integrate_rpca(obj, 1:6, 30, 0.3)
  obj <- attach_labels_from_rds(
    obj, file.path(alex_root, "dataset2/objects/d2_annotated_final.RDS"),
    c("cell_type_grouped", "cell_type_grouped"), c("cell_type", "cell_type_grouped")
  )
  obj$analysis_group <- factor(paste(substr(obj$age, 1, 1), obj$treatment),
                               levels = c("Y UT", "Y uPAR", "O UT", "O uPAR"))
  obj <- score_object(join_and_normalize(obj), signatures_mouse)
  plot_cell_umap(obj, "cell_type_grouped", file.path(outdir, "cell_type_umap.pdf"), "Post-treatment brain")
  bam <- obj@meta.data %>% mutate(cell = rownames(.)) %>% filter(grepl("BAM", cell_type, ignore.case = TRUE))
  stats <- make_split_violin(
    bam, "Myeloid_Inflammatory_Average", "analysis_group", c("O UT", "O uPAR"),
    "Myeloid inflammatory signature in BAMs", "Mean LogNormalized expression",
    file.path(outdir, "Figure_2G_myeloid_inflammatory_BAM.pdf"), c("#08306B", "#FF4902")
  )
  write_csv(stats, file.path(outdir, "Figure_2G_statistics.csv"))
  deg <- run_deg(obj, "analysis_group", "O uPAR", "O UT", file.path(outdir, "DEG"), "Old_uPAR_vs_Old_UT")
  run_enrichr(deg, file.path(outdir, "DEG"), "Old_uPAR_vs_Old_UT")
  saveRDS(obj, file.path(outdir, "post_treatment_brain_processed.rds"), compress = FALSE)
}

run_bone_marrow <- function() {
  outdir <- file.path(output_root, "03_bone_marrow")
  ids <- c("Amor_AH05_Y_UT_1_M", "Amor_AH05_Y_UT_2_M", "Amor_AH05_Y_uPAR_1_M", "Amor_AH05_Y_uPAR_2_M",
           "Amor_AH06_O_UT_1_M", "Amor_AH06_O_UT_2_M", "Amor_AH06_O_uPAR_1_M", "Amor_AH06_O_uPAR_2_M")
  samples <- data.frame(
    sample = ids,
    path = file.path(alex_root, "dataset3/counts/bone_marrow", ids, "filtered_feature_bc_matrix"),
    age = ifelse(grepl("_Y_", ids), "Young", "Old"),
    treatment = ifelse(grepl("uPAR", ids), "uPAR", "UT"),
    replicate = sub(".*_([12])_M$", "\\1", ids)
  )
  obj <- qc_filter(make_objects(samples), 10, 500, 60000, 250, 6000, outdir)
  obj <- integrate_rpca(obj, 1:20, 150, 0.5)
  obj <- attach_labels(
    obj,
    file.path(alex_root, "dataset3/bone_marrow_updated_march/results/baccin_reference_transfer/baccin_reference_label_transfer_per_cell.csv"),
    "baccin_reference_label", "cell_type"
  )
  obj$analysis_group <- factor(paste(substr(obj$age, 1, 1), obj$treatment),
                               levels = c("Y UT", "Y uPAR", "O UT", "O uPAR"))
  obj <- score_object(join_and_normalize(obj), signatures_mouse)
  plot_cell_umap(obj, "cell_type", file.path(outdir, "cell_type_umap.pdf"), "Progenitor-enriched bone marrow")
  monocytes <- obj@meta.data %>% mutate(cell = rownames(.)) %>% filter(grepl("Monocyte", cell_type, ignore.case = TRUE))
  stats <- make_split_violin(
    monocytes, "Chambers_Inflammaging1", "analysis_group", c("O UT", "O uPAR"),
    "Chambers inflammaging signature in monocytes", "Seurat module score",
    file.path(outdir, "Figure_3I_Chambers_monocytes.pdf"), c("#08306B", "#FF4902")
  )
  progenitors <- obj@meta.data %>% mutate(cell = rownames(.)) %>%
    filter(grepl("HSC|MPP|LMPP", cell_type, ignore.case = TRUE))
  stats <- bind_rows(
    stats,
    make_split_violin(progenitors, "Chambers_Inflammaging1", "analysis_group", c("O UT", "O uPAR"),
                      "Chambers inflammaging signature in HSC/MPP", "Seurat module score",
                      file.path(outdir, "Figure_3U_Chambers_HSC_MPP.pdf"), c("#08306B", "#FF4902")),
    make_split_violin(progenitors, "Kovtonyuk_Myeloid_Bias1", "analysis_group", c("O UT", "O uPAR"),
                      "Kovtonyuk myeloid-bias signature in HSC/MPP", "Seurat module score",
                      file.path(outdir, "Figure_3V_Kovtonyuk_HSC_MPP.pdf"), c("#08306B", "#FF4902"))
  )
  gmps <- obj@meta.data %>% mutate(cell = rownames(.)) %>%
    filter(grepl("GMP|Gran/Mono prog", cell_type, ignore.case = TRUE))
  stats <- bind_rows(
    stats,
    make_split_violin(gmps, "Chambers_Inflammaging1", "analysis_group", c("O UT", "O uPAR"),
                      "Chambers inflammaging signature in GMPs", "Seurat module score",
                      file.path(outdir, "Figure_3W_Chambers_GMP.pdf"), c("#08306B", "#FF4902")),
    make_split_violin(gmps, "Kovtonyuk_Myeloid_Bias1", "analysis_group", c("O UT", "O uPAR"),
                      "Kovtonyuk myeloid-bias signature in GMPs", "Seurat module score",
                      file.path(outdir, "Figure_3X_Kovtonyuk_GMP.pdf"), c("#08306B", "#FF4902"))
  )
  write_csv(stats, file.path(outdir, "signature_statistics.csv"))
  deg <- run_deg(obj, "analysis_group", "O uPAR", "O UT", file.path(outdir, "DEG"), "Old_uPAR_vs_Old_UT")
  run_enrichr(deg, file.path(outdir, "DEG"), "Old_uPAR_vs_Old_UT")
  saveRDS(obj, file.path(outdir, "bone_marrow_processed.rds"), compress = FALSE)
}

read_h5ad_as_seurat <- function(h5ad, cache) {
  if (file.exists(cache)) return(readRDS(cache))
  if (!requireNamespace("anndataR", quietly = TRUE)) {
    stop("Install anndataR to import ", h5ad, ". The analysis does not alter published annotations or embeddings.")
  }
  obj <- anndataR::read_h5ad(h5ad, as = "Seurat")
  saveRDS(obj, cache, compress = FALSE)
  obj
}

run_alra_if_needed <- function(obj, assay_name = "alra") {
  if (assay_name %in% Assays(obj)) return(obj)
  if (!requireNamespace("SeuratWrappers", quietly = TRUE)) stop("SeuratWrappers is required for ALRA")
  DefaultAssay(obj) <- "RNA"
  obj <- SeuratWrappers::RunALRA(obj, k = NULL, q = 10, quantile.prob = 0.001,
                                 use.mkl = FALSE, mkl.seed = -1)
  obj
}

score_human <- function(obj) {
  human <- list(
    Myeloid_Inflammatory = read_signature(signature_files$myeloid_human, "human"),
    Hallmark_Inflammatory_Response = msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
      filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% pull(gene_symbol) %>% unique()
  )
  score_object(join_and_normalize(obj), human)
}

plaur_status <- function(obj) {
  x <- GetAssayData(obj, assay = "RNA", layer = "data")
  if (!"PLAUR" %in% rownames(x)) stop("PLAUR is absent")
  factor(ifelse(x["PLAUR", ] > 0, "PLAUR+", "PLAUR-"), levels = c("PLAUR-", "PLAUR+"))
}

run_human_brain_atlas <- function() {
  outdir <- file.path(output_root, "04_human_brain_atlas")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  h5ad <- file.path(alex_root, "human/allen/allen_non_neuronal.h5ad")
  obj <- read_h5ad_as_seurat(h5ad, file.path(outdir, "allen_non_neuronal_seurat.rds"))
  obj <- score_human(obj)
  obj$PLAUR_status <- plaur_status(obj)
  class_col <- c("subclass", "class", "cell_type")[c("subclass", "class", "cell_type") %in% names(obj@meta.data)][1]
  if (is.na(class_col)) stop("Published Allen cell class was not retained during H5AD import")
  cns <- obj@meta.data %>% mutate(cell = rownames(.), published_class = .data[[class_col]]) %>%
    filter(grepl("CNS macrophage", published_class, ignore.case = TRUE))
  stats <- make_split_violin(cns, "Myeloid_Inflammatory1", "PLAUR_status", c("PLAUR-", "PLAUR+"),
                             "Myeloid inflammatory signature in CNS macrophages", "Seurat module score",
                             file.path(outdir, "Figure_2P_CNS_macrophages.pdf"))
  write_csv(stats, file.path(outdir, "Figure_2P_statistics.csv"))
  deg <- run_deg(obj, "PLAUR_status", "PLAUR+", "PLAUR-", file.path(outdir, "DEG"), "PLAURplus_vs_PLAURminus")
  run_enrichr(deg, file.path(outdir, "DEG"), "PLAURplus_vs_PLAURminus")
}

run_freshmg <- function() {
  outdir <- file.path(output_root, "05_FreshMG")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  validated <- file.path(alex_root, "human/results/plaur_agegroup_subclass_alra_rwb/freshmg_alra_scored_object.rds")
  obj <- if (file.exists(validated)) readRDS(validated) else
    read_h5ad_as_seurat(file.path(alex_root, "human/data/DataS1_FreshMG_final.h5ad"),
                        file.path(outdir, "FreshMG_seurat.rds"))
  obj <- score_human(obj)
  obj$PLAUR_status <- plaur_status(obj)
  if (!"PLAUR_alra" %in% names(obj@meta.data) && !"PLAUR_ALRA" %in% names(obj@meta.data)) {
    obj <- run_alra_if_needed(obj)
  }
  meta <- obj@meta.data
  age_col <- c("age", "age_num", "age_group")[c("age", "age_num", "age_group") %in% names(meta)][1]
  class_col <- c("subclass", "class", "cell_type")[c("subclass", "class", "cell_type") %in% names(meta)][1]
  if (anyNA(c(age_col, class_col))) stop("FreshMG published age/cell annotations were not retained")
  meta$age_plot <- as.character(meta[[age_col]])
  meta$published_class <- as.character(meta[[class_col]])
  meta$PLAUR_ALRA <- if ("PLAUR_alra" %in% names(meta)) meta$PLAUR_alra else if ("PLAUR_ALRA" %in% names(meta))
    meta$PLAUR_ALRA else GetAssayData(obj, assay = "alra", layer = "data")["PLAUR", ]
  pvm <- meta %>% filter(grepl("PVM", published_class, ignore.case = TRUE))
  young <- unique(pvm$age_plot[grepl("20|29|young", pvm$age_plot, ignore.case = TRUE)])
  old <- unique(pvm$age_plot[grepl("90|old", pvm$age_plot, ignore.case = TRUE)])
  if (!length(young) || !length(old)) stop("FreshMG young/old labels were not recognized")
  pvm$age_group <- ifelse(pvm$age_plot %in% young, "Young", ifelse(pvm$age_plot %in% old, "Old", NA))
  stats <- make_split_violin(pvm, "PLAUR_ALRA", "age_group", c("Young", "Old"),
                             "ALRA-imputed PLAUR in FreshMG PVMs", "ALRA-imputed PLAUR expression",
                             file.path(outdir, "Figure_2Q_PLAUR_PVM_young_old.pdf"))
  old_pvm <- pvm %>% filter(age_group == "Old")
  stats <- bind_rows(stats, make_split_violin(
    old_pvm, "Myeloid_Inflammatory1", "PLAUR_status", c("PLAUR-", "PLAUR+"),
    "Myeloid inflammatory signature in old FreshMG PVMs", "Seurat module score",
    file.path(outdir, "Figure_2R_old_PVM_myeloid_inflammatory.pdf")))
  write_csv(stats, file.path(outdir, "Figure_2Q_R_statistics.csv"))
}

run_psychad <- function() {
  outdir <- file.path(output_root, "06_PsychAD")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(alex_root, "human/results/psychad_final_requested_results_age51plus/PsychAD_July_exact_methods_control_age51plus_cache.rds")
  if (!file.exists(cache)) stop("Validated PsychAD cell cache not found: ", cache)
  cells <- readRDS(cache)
  if (inherits(cells, "Seurat")) cells <- cells@meta.data
  cells <- as.data.frame(cells)
  condition_col <- c("condition", "dx")[c("condition", "dx") %in% names(cells)][1]
  class_col <- c("class", "subclass", "cell_type")[c("class", "subclass", "cell_type") %in% names(cells)][1]
  score_col <- c("Myeloid_Inflammatory", "Myeloid_Inflammatory1")[c("Myeloid_Inflammatory", "Myeloid_Inflammatory1") %in% names(cells)][1]
  alra_col <- c("PLAUR_ALRA", "PLAUR_alra")[c("PLAUR_ALRA", "PLAUR_alra") %in% names(cells)][1]
  if (anyNA(c(condition_col, class_col, score_col, alra_col))) stop("PsychAD validated cache columns changed")
  cells$condition_plot <- ifelse(grepl("AD|Alzheimer", cells[[condition_col]], ignore.case = TRUE), "AD", "Control")
  cells$published_class <- cells[[class_col]]
  cells$score <- cells[[score_col]]
  cells$plaur_alra <- cells[[alra_col]]
  cells$PLAUR_status <- factor(ifelse(grepl("\\+", cells$PLAUR_status), "PLAUR+", "PLAUR-"),
                               levels = c("PLAUR-", "PLAUR+"))
  pvm <- cells %>% filter(grepl("PVM", published_class, ignore.case = TRUE))
  stats <- make_split_violin(pvm, "plaur_alra", "condition_plot", c("Control", "AD"),
                             "ALRA-imputed PLAUR in PsychAD PVMs", "ALRA-imputed PLAUR expression",
                             file.path(outdir, "Figure_2S_PLAUR_PVM_Control_AD.pdf"), c("#4C78A8", "#E45756"))
  ad <- pvm %>% filter(condition_plot == "AD")
  stats <- bind_rows(stats, make_split_violin(
    ad, "score", "PLAUR_status", c("PLAUR-", "PLAUR+"),
    "Myeloid inflammatory signature in PsychAD AD PVMs", "Seurat module score",
    file.path(outdir, "Figure_2T_AD_PVM_myeloid_inflammatory.pdf")))
  write_csv(stats, file.path(outdir, "Figure_2S_T_statistics.csv"))
}

run_ainciburu <- function() {
  outdir <- file.path(output_root, "07_Ainciburu_GSE180298")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  source_rds <- file.path(alex_root, "bone_marrow_human/ainciburu/ainciburu_GSE180298_seurat.rds")
  obj <- readRDS(source_rds)
  obj <- score_human(obj)
  obj$PLAUR_status <- plaur_status(obj)
  obj <- run_alra_if_needed(obj)
  meta <- obj@meta.data
  population_col <- c("CellType", "population", "cell_type", "celltype", "label")[c("CellType", "population", "cell_type", "celltype", "label") %in% names(meta)][1]
  age_col <- c("age_group", "age", "condition")[c("age_group", "age", "condition") %in% names(meta)][1]
  if (anyNA(c(population_col, age_col))) stop("Ainciburu population/age metadata columns changed")
  meta$population_plot <- as.character(meta[[population_col]])
  meta$age_plot <- ifelse(grepl("young", meta[[age_col]], ignore.case = TRUE), "Young",
                          ifelse(grepl("elder", meta[[age_col]], ignore.case = TRUE), "Elderly", NA))
  meta$PLAUR_ALRA <- GetAssayData(obj, assay = "alra", layer = "data")["PLAUR", ]
  populations <- list(LMPP_MPP = "LMPP|MPP", Monocytes = "Monocyte", HSC = "HSC", GMP = "GMP")
  stats <- list()
  for (nm in names(populations)) {
    d <- meta %>% filter(!is.na(age_plot), grepl(populations[[nm]], population_plot, ignore.case = TRUE))
    stats[[paste0(nm, "_age")]] <- make_split_violin(
      d, "PLAUR_ALRA", "age_plot", c("Young", "Elderly"),
      paste("ALRA-imputed PLAUR in", nm), "ALRA-imputed PLAUR expression",
      file.path(outdir, paste0(nm, "_PLAUR_young_elderly.pdf")))
    old <- d %>% filter(age_plot == "Elderly")
    stats[[paste0(nm, "_status")]] <- make_split_violin(
      old, "Hallmark_Inflammatory_Response1", "PLAUR_status", c("PLAUR-", "PLAUR+"),
      paste("Hallmark inflammatory response in elderly", nm), "Seurat module score",
      file.path(outdir, paste0(nm, "_hallmark_inflammation_PLAUR_status.pdf")))
  }
  write_csv(bind_rows(stats, .id = "panel"), file.path(outdir, "Ainciburu_statistics.csv"))
}

runners <- list(
  sorted_brain = run_sorted_brain,
  post_treatment_brain = run_post_treatment_brain,
  bone_marrow = run_bone_marrow,
  human_brain_atlas = run_human_brain_atlas,
  freshmg = run_freshmg,
  psychad = run_psychad,
  ainciburu = run_ainciburu
)

validate_inputs <- function() {
  paths <- c(
    signature_mouse = signature_files$myeloid_mouse,
    signature_human = signature_files$myeloid_human,
    sorted_brain_labels = file.path(alex_root, "dataset1/seurat_objs/annotated_cell_types_d1.RDS"),
    post_treatment_brain_labels = file.path(alex_root, "dataset2/objects/d2_annotated_final.RDS"),
    bone_marrow_labels = file.path(alex_root, "dataset3/bone_marrow_updated_march/results/baccin_reference_transfer/baccin_reference_label_transfer_per_cell.csv"),
    bm_signatures = file.path(alex_root, "BM signatures.xlsx"),
    sorted_brain_matrix = file.path(alex_root, "dataset1/Amor_CAV17_Female_UPARneg/filtered_feature_bc_matrix"),
    post_treatment_brain_matrix = file.path(alex_root, "dataset2/Y_uPAR_F/filtered_feature_bc_matrix.h5"),
    bone_marrow_matrix = file.path(alex_root, "dataset3/counts/bone_marrow/Amor_AH05_Y_UT_1_M/filtered_feature_bc_matrix"),
    allen = file.path(alex_root, "human/allen/allen_non_neuronal.h5ad"),
    freshmg = file.path(alex_root, "human/data/DataS1_FreshMG_final.h5ad"),
    psychad_cache = file.path(alex_root, "human/results/psychad_final_requested_results_age51plus/PsychAD_July_exact_methods_control_age51plus_cache.rds"),
    ainciburu = file.path(alex_root, "bone_marrow_human/ainciburu/ainciburu_GSE180298_seurat.rds")
  )
  status <- data.frame(input = names(paths), path = unname(paths), exists = file.exists(paths))
  write_csv(status, file.path(project_root, "input_validation.csv"))
  print(status, row.names = FALSE)
  if (!all(status$exists)) stop("One or more required inputs are missing")
  message("All manuscript inputs and support files were found.")
}

run_manuscript_dataset <- function(dataset = "all") {
  if (!dataset %in% allowed) stop("Unknown dataset: ", dataset)
  if (dataset == "validate") return(validate_inputs())
  selected <- if (dataset == "all") names(runners) else dataset
  for (name in selected) {
    message("Running ", name)
    runners[[name]]()
    gc()
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
  message("Complete: ", output_root)
}

if (sys.nframe() == 0) run_manuscript_dataset(dataset)
