# PsychAD human brain myeloid cells
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

plaur_status <- function(obj) {
  x <- GetAssayData(obj, assay = "RNA", layer = "data")
  if (!"PLAUR" %in% rownames(x))
    stop("PLAUR is absent")
  factor(ifelse(x["PLAUR", ] > 0, "PLAUR+", "PLAUR-"), levels = c("PLAUR-", "PLAUR+"))
}

# Analysis

outdir <- file.path(output_root, "06_PsychAD")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(alex_root, "human/data/DataS2_PsychAD_final.h5ad")

if (!requireNamespace("anndataR", quietly = TRUE)) {
  stop("Install anndataR to import ", input_file, ". The analysis does not alter published annotations or embeddings.")
}

# Read the data and perform the manuscript processing workflow
obj <- anndataR::read_h5ad(input_file, as = "Seurat")

saveRDS(obj, file.path(outdir, "PsychAD_seurat.rds"), compress = FALSE)

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

if (!"alra" %in% Assays(obj)) {
  if (!requireNamespace("SeuratWrappers", quietly = TRUE))
    stop("SeuratWrappers is required for ALRA")
  DefaultAssay(obj) <- "RNA"
  obj <- SeuratWrappers::RunALRA(obj, k = NULL, q = 10, quantile.prob = 0.001, use.mkl = FALSE, mkl.seed = -1)
}

obj <- obj

cells <- obj@meta.data

condition_col <- c("condition", "dx")[c("condition", "dx") %in% names(cells)][1]

class_col <- c("class", "subclass", "cell_type")[c("class", "subclass", "cell_type") %in% names(cells)][1]

age_col <- c("age_num", "age")[c("age_num", "age") %in% names(cells)][1]

if (anyNA(c(condition_col, class_col, age_col))) stop("PsychAD published metadata were not retained")

cells$condition_plot <- ifelse(grepl("AD|Alzheimer", cells[[condition_col]], ignore.case = TRUE), "AD", "Control")

cells$published_class <- cells[[class_col]]

cells$age_numeric <- suppressWarnings(as.numeric(sub("\\+", "", as.character(cells[[age_col]]))))

cells$score <- cells$Myeloid_Inflammatory1

cells$plaur_alra <- GetAssayData(obj, assay = "alra", layer = "data")["PLAUR", ]

age_matched <- cells %>% filter(condition_plot == "AD" | age_numeric >= 51)

pvm <- age_matched %>% filter(grepl("PVM", published_class, ignore.case = TRUE))

# Manuscript figures and cell-level Wilcoxon tests
stats <- make_split_violin(age_matched, "plaur_alra", "condition_plot", c("Control", "AD"), "ALRA-imputed PLAUR in age-matched PsychAD",
  "ALRA-imputed PLAUR expression", file.path(outdir, "Figure_S2R_PLAUR_all_cells_Control_AD.pdf"), c("#4C78A8", "#E45756"))

# Manuscript figures and cell-level Wilcoxon tests
stats <- bind_rows(stats, make_split_violin(pvm, "plaur_alra", "condition_plot", c("Control", "AD"), "ALRA-imputed PLAUR in PsychAD PVMs",
  "ALRA-imputed PLAUR expression", file.path(outdir, "Figure_2S_PLAUR_PVM_Control_AD.pdf"), c("#4C78A8", "#E45756")))

ad <- pvm %>% filter(condition_plot == "AD")

# Manuscript figures and cell-level Wilcoxon tests
stats <- bind_rows(stats, make_split_violin(ad, "score", "PLAUR_status", c("PLAUR-", "PLAUR+"), "Myeloid inflammatory signature in PsychAD AD PVMs",
  "Seurat module score", file.path(outdir, "Figure_2T_AD_PVM_myeloid_inflammatory.pdf")))

write_csv(stats, file.path(outdir, "Figure_2S_T_statistics.csv"))

age_cells <- rownames(cells)[cells$condition_plot == "AD" | cells$age_numeric >= 51]

age_obj <- subset(obj, cells = age_cells)

Idents(age_obj) <- factor(age_obj[["PLAUR_status", drop = TRUE]])

# Differential expression and Enrichr analysis
deg <- FindMarkers(age_obj, ident.1 = "PLAUR+", ident.2 = "PLAUR-", assay = "RNA", test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0,
  only.pos = FALSE)

deg$gene <- rownames(deg)

fc <- intersect(c("avg_log2FC", "avg_logFC"), names(deg))[[1]]

keep <- deg$p_val < 0.05 & abs(deg[[fc]]) > 0.5

filtered <- deg[keep, , drop = FALSE]

filtered$direction <- ifelse(filtered[[fc]] > 0, "UP", "DOWN")

dir.create(file.path(outdir, "DEG"), recursive = TRUE, showWarnings = FALSE)

write_csv(deg, file.path(file.path(outdir, "DEG"), paste0("Figure_S2S_PLAURplus_vs_PLAURminus", "_all_genes.csv")))

write_csv(filtered, file.path(file.path(outdir, "DEG"), paste0("Figure_S2S_PLAURplus_vs_PLAURminus", "_filtered_for_enrichr.csv")))

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
  write.xlsx(results, file.path(file.path(outdir, "DEG"), paste0("Figure_S2S_PLAURplus_vs_PLAURminus", "_enrichr.xlsx")),
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
    save_pdf(p, file.path(file.path(outdir, "DEG"), paste0("Figure_S2S_PLAURplus_vs_PLAURminus", "_enrichr.pdf")), 10,
      9)
    write_csv(plot_data %>% mutate(term_plot = as.character(term_plot)), file.path(file.path(outdir, "DEG"), paste0("Figure_S2S_PLAURplus_vs_PLAURminus",
      "_enrichr_plot_data.csv")))
  }
}

invisible(results)

writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message("Complete: ", output_root)
