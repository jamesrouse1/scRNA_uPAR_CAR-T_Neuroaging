# FACS-sorted aged mouse brain
# Single-cell RNA-seq analysis used in the manuscript

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(data.table)
library(readxl)
library(openxlsx)

set.seed(1)
options(future.globals.maxSize = 16 * 1024^3)
future::plan("sequential")
if (exists("mem.maxVSize")) mem.maxVSize(130000)

# File locations
script_arg <- commandArgs()[grepl("^--file=", commandArgs())]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else normalizePath(".")
project_root <- dirname(normalizePath(script_file, mustWork = FALSE))
alex_root <- normalizePath(Sys.getenv("ALEX_ROOT", "/Users/rouse/CSH/alex"), mustWork = FALSE)
output_root <- file.path(project_root, "outputs")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

save_pdf <- function(plot, path, width = 8, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = width, height = height, family = "Helvetica", useDingbats = FALSE, onefile = FALSE, compress = TRUE,
    bg = "white", version = "1.4")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}

read_10x_any <- function(path) {
  x <- if (dir.exists(path))
    Read10X(path, gene.column = 2, unique.features = TRUE)
  else Read10X_h5(path)
  if (is.list(x)) {
    if (!"Gene Expression" %in% names(x))
      stop("Gene Expression matrix not found in ", path)
    x <- x[["Gene Expression"]]
  }
  x
}

read_signature <- function(file, species = c("mouse", "human")) {
  species <- match.arg(species)
  genes <- trimws(readLines(file, warn = FALSE))
  unique(genes[nzchar(genes)])
}

signature_files <- list(myeloid_mouse = file.path(alex_root, "manuscript_signature_validation_20260820", "tables", "new_signature_71_mouse.txt"),
  myeloid_human = file.path(alex_root, "manuscript_signature_validation_20260820", "tables", "new_signature_human_orthologs.txt"))

read_bm_signature <- function(column) {
  x <- read_excel(file.path(alex_root, "BM signatures.xlsx"), col_names = FALSE)
  unique(na.omit(trimws(as.character(x[[column]][-(1:2)]))))
}

signatures_mouse <- list(Myeloid_Inflammatory = read_signature(signature_files$myeloid_mouse, "mouse"), Hallmark_Inflammatory_Response = msigdbr::msigdbr(species = "Mus musculus",
  collection = "H") %>% filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% pull(gene_symbol) %>% unique(), Chambers_Inflammaging = read_bm_signature(4),
Kovtonyuk_Myeloid_Bias = read_bm_signature(7))

split_violin_geom <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  group <- data[1, "group"]
  data2 <- transform(data, x = if (group %% 2 == 1)
    xminv
  else xmaxv)
  data2 <- data2[order(if (group %% 2 == 1)
    data2$y
  else -data2$y), ]
  data2 <- rbind(data2[1, ], data2, data2[nrow(data2), ], data2[1, ])
  data2[c(1, nrow(data2) - 1, nrow(data2)), "x"] <- round(data2[1, "x"])
  ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(data2, ...))
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., trim = TRUE, scale = "area",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = split_violin_geom, position = position, show.legend = show.legend,
    inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, na.rm = na.rm, ...))
}

format_p <- function(p) {
  if (!is.finite(p))
    "NA"
  else format.pval(p, digits = 3, eps = 9.99999999999999e-301)
}

make_split_violin <- function(df, value, group, levels, title, ylab, path, colors = c("#7CA1CC", "#FF4902")) {
  d <- df %>% filter(.data[[group]] %in% levels, is.finite(.data[[value]]))
  d$plot_group <- factor(d[[group]], levels = levels)
  d$comparison <- factor(paste(levels, collapse = " vs "))
  d$side <- ifelse(d$plot_group == levels[[1]], 1L, 2L)
  group_counts <- table(d$plot_group)
  legend_labels <- paste0(levels, " (n=", as.integer(group_counts[levels]), ")")
  pvalue <- wilcox.test(d[[value]] ~ d$plot_group, exact = FALSE)$p.value
  ymax <- max(d[[value]], na.rm = TRUE)
  yrange <- diff(range(d[[value]], na.rm = TRUE))
  p <- ggplot(d, aes(comparison, .data[[value]], fill = plot_group, group = side)) + geom_split_violin(trim = FALSE, scale = "width",
    color = "black", linewidth = 0.25) + geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35, color = "black") +
    annotate("text", x = 1, y = ymax + max(0.08 * yrange, 0.03), label = paste0("p = ", format_p(pvalue)), size = 3.4) +
    scale_fill_manual(values = setNames(colors, levels), labels = legend_labels) + scale_y_continuous(expand = expansion(mult = c(0.03,
      0.17))) + theme_classic(base_size = 11) + labs(title = title, x = NULL, y = ylab, fill = NULL)
  save_pdf(p, path, 7.5, 5.5)
  data.frame(group_1 = levels[[1]], n_1 = sum(d$plot_group == levels[[1]]), group_2 = levels[[2]], n_2 = sum(d$plot_group ==
    levels[[2]]), wilcox_p = pvalue)
}

# Analysis

outdir <- file.path(output_root, "01_sorted_brain")

# Sample information and filtered 10x matrices
samples <- data.frame(sample = c("F_neg", "F_plus", "M_neg", "M_plus"), path = file.path(alex_root, "dataset1", c("Amor_CAV17_Female_UPARneg/filtered_feature_bc_matrix",
  "Amor_CAV17_Female_UPARpos/filtered_feature_bc_matrix", "Amor_CAV16_Male_UPARneg/filtered_feature_bc_matrix", "Amor_CAV16_Male_UPAR_pos/filtered_feature_bc_matrix")),
treatment = c("uPAR-", "uPAR+", "uPAR-", "uPAR+"), sex = c("Female", "Female", "Male", "Male"))

objects <- lapply(seq_len(nrow(samples)), function(i) {
  sample <- samples$sample[[i]]
  x <- CreateSeuratObject(read_10x_any(samples$path[[i]]), project = sample)
  x$sample <- sample
  for (field in setdiff(names(samples), c("sample", "path"))) {
    x[[field]] <- samples[[field]][[i]]
  }
  x
})

names(objects) <- samples$sample

# Read the data and perform the manuscript processing workflow
obj <- merge(objects[[1]], y = objects[-1], add.cell.ids = samples$sample)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

before <- as.data.frame(table(obj$sample), stringsAsFactors = FALSE)

names(before) <- c("sample", "cells_before")

p_before <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 1, pt.size = 0, layer = "counts")

save_pdf(p_before, file.path(outdir, "01_qc_before_filtering.pdf"), 13, 12)

keep <- obj$percent.mt < 15 & obj$nCount_RNA > 500 & obj$nCount_RNA < 20000 & obj$nFeature_RNA > 250 & obj$nFeature_RNA <
  5000

obj <- subset(obj, cells = colnames(obj)[keep])

after <- as.data.frame(table(obj$sample), stringsAsFactors = FALSE)

names(after) <- c("sample", "cells_after")

summary <- full_join(before, after, by = "sample") %>% mutate(across(starts_with("cells_"), ~ replace_na(.x, 0L)), cells_removed = cells_before -
  cells_after, percent_retained = 100 * cells_after / cells_before)

write_csv(summary, file.path(outdir, "qc_cell_counts.csv"))

p_after <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 1, pt.size = 0, layer = "counts")

save_pdf(p_after, file.path(outdir, "02_qc_after_filtering.pdf"), 13, 12)

obj <- obj

pieces <- SplitObject(obj, split.by = "sample")

pieces <- lapply(pieces, SCTransform, verbose = FALSE)

features <- SelectIntegrationFeatures(pieces, nfeatures = 2000)

pieces <- PrepSCTIntegration(pieces, anchor.features = features, verbose = FALSE)

anchors <- FindIntegrationAnchors(pieces, normalization.method = "SCT", anchor.features = features, verbose = FALSE)

obj <- IntegrateData(anchors, normalization.method = "SCT", verbose = FALSE)

obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)

obj <- FindClusters(obj, resolution = 1, verbose = FALSE)

obj <- RunUMAP(obj, dims = 1:5, n.neighbors = 30, min.dist = 0.05, seed.use = 1, verbose = FALSE)

annotated <- readRDS(file.path(alex_root, "dataset1", "seurat_objs", "annotated_cell_types_d1.RDS"))

idx <- match(colnames(obj), colnames(annotated))

if (anyNA(idx)) stop(sum(is.na(idx)), " cells are absent from ", file.path(alex_root, "dataset1", "seurat_objs", "annotated_cell_types_d1.RDS"))

for (field in c("cell_type", "cell_type_grouped")) obj[[field]] <- annotated@meta.data[[field]][idx]

rm(annotated)

gc()

obj <- obj

DefaultAssay(obj) <- "RNA"

if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

for (nm in names(signatures_mouse)) {
  genes <- intersect(signatures_mouse[[nm]], rownames(obj))
  if (length(genes) < 3)
    stop("Too few genes found for ", nm)
  set.seed(1)
  obj <- AddModuleScore(obj, features = list(genes), name = paste0("", nm), assay = "RNA", nbin = 24, ctrl = 100, seed = 1,
    search = FALSE, slot = "data")
  avg_name <- paste0("", nm, "_Average")
  obj[[avg_name]] <- Matrix::colMeans(GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE])
}

obj <- obj

p <- DimPlot(obj, reduction = "umap", group.by = "cell_type_grouped", label = TRUE, repel = TRUE, raster = FALSE, pt.size = 0.12) +
  ggtitle("Sorted aged brain") + theme_classic(base_size = 11)

save_pdf(p, file.path(outdir, "Figure_S2C_cell_type_umap.pdf"), 11, 8)

um <- as.data.frame(Embeddings(obj, "umap"))

names(um)[1:2] <- c("UMAP_1", "UMAP_2")

um$group <- factor(obj[["treatment", drop = TRUE]], levels = c("uPAR-", "uPAR+"))

um$cell_type <- as.character(obj[["cell_type_grouped", drop = TRUE]])

um <- um[!is.na(um$group), , drop = FALSE]

bins <- 100L

bx <- seq(min(um$UMAP_1), max(um$UMAP_1), length.out = bins + 1L)

by <- seq(min(um$UMAP_2), max(um$UMAP_2), length.out = bins + 1L)

um$ix <- findInterval(um$UMAP_1, bx, all.inside = TRUE)

um$iy <- findInterval(um$UMAP_2, by, all.inside = TRUE)

density <- um %>% count(ix, iy, group, name = "n") %>% complete(ix, iy, group, fill = list(n = 0)) %>% group_by(group) %>%
  mutate(value = n / sum(n)) %>% ungroup() %>% select(ix, iy, group, value) %>% pivot_wider(names_from = group, values_from = value,
    values_fill = 0) %>% mutate(x = (bx[ix] + bx[ix + 1L]) / 2, y = (by[iy] + by[iy + 1L]) / 2, difference = .data[[c("uPAR-",
    "uPAR+")[[2]]]] - .data[[c("uPAR-", "uPAR+")[[1]]]])

if (FALSE) {
  z <- matrix(0, bins, bins)
  z[cbind(density$ix, density$iy)] <- density$difference
  kernel <- outer(exp(-(-2:2)^2 / 2), exp(-(-2:2)^2 / 2))
  kernel <- kernel / sum(kernel)
  zs <- matrix(0, bins, bins)
  for (i in seq_len(bins)) for (j in seq_len(bins)) {
    ii <- pmax(1, pmin(bins, i + (-2:2)))
    jj <- pmax(1, pmin(bins, j + (-2:2)))
    zs[i, j] <- sum(z[ii, jj, drop = FALSE] * kernel)
  }
  density$difference <- zs[cbind(density$ix, density$iy)]
}

labels <- um %>% group_by(cell_type) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

lim <- max(abs(density$difference), na.rm = TRUE)

p <- ggplot(density, aes(x, y, fill = difference)) + geom_tile(linewidth = 0) + geom_text(data = labels, aes(UMAP_1, UMAP_2,
  label = cell_type), inherit.aes = FALSE, size = 2.7, fontface = "bold", check_overlap = TRUE) + scale_fill_gradient2(low = "#2C7BB6",
  mid = "white", high = "#D7191C", midpoint = 0, limits = c(-lim, lim), oob = scales::squish) + coord_equal(expand = FALSE) +
  theme_classic(base_size = 11) + labs(x = "UMAP 1", y = "UMAP 2", fill = paste0(c("uPAR-", "uPAR+")[[2]], " -\n", c("uPAR-",
    "uPAR+")[[1]]))

save_pdf(p, file.path(outdir, "Figure_2B_binned_density_difference.pdf"), 9, 7)

write_csv(density, sub("\\.pdf$", "_data.csv", file.path(outdir, "Figure_2B_binned_density_difference.pdf")))

d <- obj@meta.data %>% transmute(cell_type = as.character(.data[["cell_type_grouped"]]), group = factor(.data[["treatment"]],
  levels = c("uPAR-", "uPAR+"))) %>% filter(!is.na(group), !is.na(cell_type), nzchar(cell_type)) %>% count(cell_type, group,
  name = "n_cells") %>% complete(cell_type, group, fill = list(n_cells = 0)) %>% group_by(group) %>% mutate(fraction_within_group = n_cells / sum(n_cells)) %>%
  ungroup() %>% group_by(cell_type) %>% mutate(fraction_within_cell_type = fraction_within_group / sum(fraction_within_group)) %>%
  ungroup() %>% mutate(cell_type = factor(cell_type, levels = sort(unique(cell_type))))

p1 <- ggplot(d, aes(cell_type, fraction_within_group, fill = group)) + geom_col(position = position_dodge(width = 0.82),
  width = 0.72, color = "black", linewidth = 0.2) + scale_fill_manual(values = c(`uPAR-` = "#7CA1CC", `uPAR+` = "#FF4902"),
  drop = FALSE) + scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") + labs(x = NULL,
    y = "Fraction of all cells in group", fill = NULL)

p2 <- ggplot(d, aes(cell_type, fraction_within_cell_type, fill = group)) + geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = c(`uPAR-` = "#7CA1CC", `uPAR+` = "#FF4902"), drop = FALSE) + scale_y_continuous(labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) + theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 45,
    hjust = 1), legend.position = "top") + labs(x = NULL, y = "Fraction within cell type", fill = NULL)

save_pdf(p1, file.path(outdir, paste0("Figure_S2D_uPAR_status", "_fraction_within_group.pdf")), 12, 6)

save_pdf(p2, file.path(outdir, paste0("Figure_S2D_uPAR_status", "_fraction_within_cell_type.pdf")), 12, 6)

write_csv(d, file.path(outdir, paste0("Figure_S2D_uPAR_status", "_cell_counts_and_fractions.csv")))

bam <- obj@meta.data %>% mutate(cell = rownames(.), treatment = factor(treatment, c("uPAR-", "uPAR+"))) %>% filter(grepl("BAM",
  cell_type, ignore.case = TRUE))

# Manuscript figures and cell-level Wilcoxon tests
stats <- make_split_violin(bam, "Myeloid_Inflammatory_Average", "treatment", c("uPAR-", "uPAR+"), "Myeloid inflammatory signature in BAMs",
  "Mean LogNormalized expression", file.path(outdir, "Figure_2D_myeloid_inflammatory_BAM.pdf"))

write_csv(stats, file.path(outdir, "Figure_2D_statistics.csv"))

Idents(obj) <- factor(obj[["treatment", drop = TRUE]])

# Differential expression and Enrichr analysis
deg <- FindMarkers(obj, ident.1 = "uPAR+", ident.2 = "uPAR-", assay = "RNA", test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0,
  only.pos = FALSE)

deg$gene <- rownames(deg)

fc <- intersect(c("avg_log2FC", "avg_logFC"), names(deg))[[1]]

keep <- deg$p_val < 0.05 & abs(deg[[fc]]) > 0.5

filtered <- deg[keep, , drop = FALSE]

filtered$direction <- ifelse(filtered[[fc]] > 0, "UP", "DOWN")

dir.create(file.path(outdir, "DEG"), recursive = TRUE, showWarnings = FALSE)

write_csv(deg, file.path(file.path(outdir, "DEG"), paste0("Figure_2C_uPARplus_vs_uPARminus", "_all_genes.csv")))

write_csv(filtered, file.path(file.path(outdir, "DEG"), paste0("Figure_2C_uPARplus_vs_uPARminus", "_filtered_for_enrichr.csv")))

# Differential expression and Enrichr analysis
deg <- filtered

if (!requireNamespace("enrichR", quietly = TRUE)) stop("The enrichR package is required")

databases <- c("GO_Biological_Process_2025", "MSigDB_Hallmark_2020", "Reactome_2022")

results <- list()

for (direction in c("UP", "DOWN")) {
  genes <- unique(deg$gene[deg$direction == direction])
  if (!length(genes))
    next
  er <- enrichR::enrichr(genes, databases)
  for (db in names(er)) results[[paste(direction, db, sep = "_")]] <- er[[db]]
}

if (length(results)) {
  write.xlsx(results, file.path(file.path(outdir, "DEG"), paste0("Figure_2C_uPARplus_vs_uPARminus", "_enrichr.xlsx")),
    overwrite = TRUE)
  plot_data <- bind_rows(lapply(names(results), function(sheet) {
    x <- results[[sheet]]
    if (!nrow(x))
      return(data.frame())
    pieces <- strsplit(sheet, "_", fixed = TRUE)[[1]]
    direction <- pieces[[1]]
    database <- sub(paste0("^", direction, "_"), "", sheet)
    x %>% transmute(direction = direction, database = database, term = as.character(Term), p_value = as.numeric(P.value),
      odds_ratio = as.numeric(Odds.Ratio), combined_score = as.numeric(Combined.Score))
  })) %>% filter(is.finite(p_value), p_value > 0, is.finite(odds_ratio), odds_ratio > 0, is.finite(combined_score), combined_score >
    0) %>% group_by(direction) %>% arrange(p_value, desc(odds_ratio), .by_group = TRUE) %>% distinct(term, .keep_all = TRUE) %>%
    slice_head(n = 10) %>% ungroup() %>% mutate(term_plot = stringr::str_wrap(term, 42)) %>% arrange(direction, odds_ratio) %>%
    mutate(term_plot = factor(term_plot, levels = unique(term_plot)))
  if (nrow(plot_data)) {
    p <- ggplot(plot_data, aes(odds_ratio, term_plot)) + geom_segment(aes(x = 1, xend = odds_ratio, yend = term_plot),
      color = "#AAB4C3", linewidth = 0.5) + geom_point(aes(color = p_value, size = combined_score)) + scale_x_log10() +
      scale_color_gradient(low = "#B2182B", high = "#2166AC", trans = "log10") + scale_size_continuous(trans = "log10",
        range = c(3, 9)) + facet_grid(direction ~ ., scales = "free_y", space = "free_y") + theme_classic(base_size = 11) +
      theme(axis.text.y = element_text(size = 9), strip.text = element_text(face = "bold")) + labs(x = "Odds ratio (log scale)",
        y = NULL, color = "Nominal P value", size = "Combined score")
    save_pdf(p, file.path(file.path(outdir, "DEG"), paste0("Figure_2C_uPARplus_vs_uPARminus", "_enrichr.pdf")), 10, 9)
    write_csv(plot_data %>% mutate(term_plot = as.character(term_plot)), file.path(file.path(outdir, "DEG"), paste0("Figure_2C_uPARplus_vs_uPARminus",
      "_enrichr_plot_data.csv")))
  }
}

invisible(results)

saveRDS(obj, file.path(outdir, "sorted_brain_processed.rds"), compress = FALSE)

writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message("Complete: ", output_root)
