# Siletti adult human brain atlas
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

read_signature <- function(file, species = c("mouse", "human")) {
  species <- match.arg(species)
  genes <- trimws(readLines(file, warn = FALSE))
  unique(genes[nzchar(genes)])
}

signature_files <- list(myeloid_mouse = file.path(alex_root, "manuscript_signature_validation_20260820", "tables", "new_signature_71_mouse.txt"),
  myeloid_human = file.path(alex_root, "manuscript_signature_validation_20260820", "tables", "new_signature_human_orthologs.txt"))

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
  pvalue <- wilcox.test(d[[value]] ~ d$plot_group, exact = FALSE)$p.value
  ymax <- max(d[[value]], na.rm = TRUE)
  yrange <- diff(range(d[[value]], na.rm = TRUE))
  p <- ggplot(d, aes(comparison, .data[[value]], fill = plot_group, group = side)) + geom_split_violin(trim = FALSE, scale = "width",
    color = "black", linewidth = 0.25) + geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35, color = "black") +
    annotate("text", x = 1, y = ymax + max(0.08 * yrange, 0.03), label = paste0("p = ", format_p(pvalue)), size = 3.4) +
    scale_fill_manual(values = setNames(colors, levels)) + scale_y_continuous(expand = expansion(mult = c(0.03, 0.17))) +
    theme_classic(base_size = 11) + labs(title = title, x = NULL, y = ylab, fill = NULL)
  save_pdf(p, path, 7.5, 5.5)
  data.frame(group_1 = levels[[1]], n_1 = sum(d$plot_group == levels[[1]]), group_2 = levels[[2]], n_2 = sum(d$plot_group ==
    levels[[2]]), wilcox_p = pvalue)
}

plaur_status <- function(obj) {
  x <- GetAssayData(obj, assay = "RNA", layer = "data")
  if (!"PLAUR" %in% rownames(x))
    stop("PLAUR is absent")
  factor(ifelse(x["PLAUR", ] > 0, "PLAUR+", "PLAUR-"), levels = c("PLAUR-", "PLAUR+"))
}

# Analysis

outdir <- file.path(output_root, "04_human_brain_atlas")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

h5ad <- file.path(alex_root, "human/allen/allen_non_neuronal.h5ad")

if (!requireNamespace("anndataR", quietly = TRUE)) {
  stop("Install anndataR to import ", h5ad, ". The analysis does not alter published annotations or embeddings.")
}

# Read the data and perform the manuscript processing workflow
obj <- anndataR::read_h5ad(h5ad, as = "Seurat")

saveRDS(obj, file.path(outdir, "allen_non_neuronal_seurat.rds"), compress = FALSE)

obj <- obj

DefaultAssay(obj) <- "RNA"

if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

human <- list(Myeloid_Inflammatory = read_signature(signature_files$myeloid_human, "human"), Hallmark_Inflammatory_Response = msigdbr::msigdbr(species = "Homo sapiens",
  collection = "H") %>% filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% pull(gene_symbol) %>% unique())

for (nm in names(human)) {
  genes <- intersect(human[[nm]], rownames(obj))
  if (length(genes) < 3)
    stop("Too few genes found for ", nm)
  set.seed(1)
  obj <- AddModuleScore(obj, features = list(genes), name = nm, assay = "RNA", nbin = 24, ctrl = 100, seed = 1, search = FALSE,
    slot = "data")
  obj[[paste0(nm, "_Average")]] <- Matrix::colMeans(GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE])
}

obj <- obj

obj$PLAUR_status <- plaur_status(obj)

class_col <- c("subclass", "class", "cell_type")[c("subclass", "class", "cell_type") %in% names(obj@meta.data)][1]

if (is.na(class_col)) stop("Published Allen cell class was not retained during H5AD import")

cns <- obj@meta.data %>% mutate(cell = rownames(.), published_class = .data[[class_col]]) %>% filter(grepl("CNS macrophage",
  published_class, ignore.case = TRUE))

# Manuscript figures and cell-level Wilcoxon tests
stats <- make_split_violin(cns, "Myeloid_Inflammatory1", "PLAUR_status", c("PLAUR-", "PLAUR+"), "Myeloid inflammatory signature in CNS macrophages",
  "Seurat module score", file.path(outdir, "Figure_2P_CNS_macrophages.pdf"))

write_csv(stats, file.path(outdir, "Figure_2P_statistics.csv"))

Idents(obj) <- factor(obj[["PLAUR_status", drop = TRUE]])

# Differential expression and Enrichr analysis
deg <- FindMarkers(obj, ident.1 = "PLAUR+", ident.2 = "PLAUR-", assay = "RNA", test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0,
  only.pos = FALSE)

deg$gene <- rownames(deg)

fc <- intersect(c("avg_log2FC", "avg_logFC"), names(deg))[[1]]

keep <- deg$p_val < 0.05 & abs(deg[[fc]]) > 0.5

filtered <- deg[keep, , drop = FALSE]

filtered$direction <- ifelse(filtered[[fc]] > 0, "UP", "DOWN")

dir.create(file.path(outdir, "DEG"), recursive = TRUE, showWarnings = FALSE)

write_csv(deg, file.path(file.path(outdir, "DEG"), paste0("PLAURplus_vs_PLAURminus", "_all_genes.csv")))

write_csv(filtered, file.path(file.path(outdir, "DEG"), paste0("PLAURplus_vs_PLAURminus", "_filtered_for_enrichr.csv")))

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

if (length(results)) write.xlsx(results, file.path(file.path(outdir, "DEG"), paste0("PLAURplus_vs_PLAURminus", "_enrichr.xlsx")),
  overwrite = TRUE)

invisible(results)

writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message("Complete: ", output_root)
