#!/usr/bin/env Rscript

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
project_root <- dirname(normalizePath(script_file, mustWork = FALSE))
alex_root <- normalizePath(Sys.getenv("ALEX_ROOT", dirname(project_root)), mustWork = FALSE)
output_root <- file.path(project_root, "outputs")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

save_pdf <- function(plot, path, width = 8, height = 6) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(path, width = width, height = height, family = "Helvetica", useDingbats = FALSE, onefile = FALSE, 
        compress = TRUE, bg = "white", version = "1.4")
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
    p_before <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 1, pt.size = 0, 
        layer = "counts")
    save_pdf(p_before, file.path(outdir, "01_qc_before_filtering.pdf"), 13, 12)
    keep <- obj$percent.mt < mt_max & obj$nCount_RNA > count_min & obj$nCount_RNA < count_max & obj$nFeature_RNA > 
        feature_min & obj$nFeature_RNA < feature_max
    obj <- subset(obj, cells = colnames(obj)[keep])
    after <- as.data.frame(table(obj$sample), stringsAsFactors = FALSE)
    names(after) <- c("sample", "cells_after")
    summary <- full_join(before, after, by = "sample") %>% mutate(across(starts_with("cells_"), ~replace_na(.x, 
        0L)), cells_removed = cells_before - cells_after, percent_retained = 100 * cells_after/cells_before)
    write_csv(summary, file.path(outdir, "qc_cell_counts.csv"))
    p_after <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 1, pt.size = 0, 
        layer = "counts")
    save_pdf(p_after, file.path(outdir, "02_qc_after_filtering.pdf"), 13, 12)
    obj
}

integrate_sorted_brain <- function(obj) {
    pieces <- SplitObject(obj, split.by = "sample")
    pieces <- lapply(pieces, SCTransform, verbose = FALSE)
    features <- SelectIntegrationFeatures(pieces, nfeatures = 2000)
    pieces <- PrepSCTIntegration(pieces, anchor.features = features, verbose = FALSE)
    anchors <- FindIntegrationAnchors(pieces, normalization.method = "SCT", anchor.features = features, verbose = FALSE)
    obj <- IntegrateData(anchors, normalization.method = "SCT", verbose = FALSE)
    obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 1, verbose = FALSE)
    RunUMAP(obj, dims = 1:5, n.neighbors = 30, min.dist = 0.05, seed.use = 1, verbose = FALSE)
}

integrate_rpca <- function(obj, dims, neighbors, min_dist, resolution = 0.5) {
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
    obj <- SCTransform(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
    obj <- IntegrateLayers(obj, method = RPCAIntegration, normalization.method = "SCT", verbose = FALSE)
    obj <- FindNeighbors(obj, dims = dims, reduction = "integrated.dr", verbose = FALSE)
    obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
    RunUMAP(obj, dims = dims, reduction = "integrated.dr", n.neighbors = neighbors, min.dist = min_dist, seed.use = 1, 
        verbose = FALSE)
}

attach_labels <- function(obj, file, label_columns, target_columns = label_columns) {
    labels <- fread(file, data.table = FALSE)
    idx <- match(colnames(obj), labels$cell_id)
    if (anyNA(idx)) 
        stop(sum(is.na(idx)), " cells are absent from ", file)
    for (i in seq_along(label_columns)) {
        obj[[target_columns[[i]]]] <- labels[[label_columns[[i]]]][idx]
    }
    obj
}

attach_labels_from_rds <- function(obj, file, source_columns, target_columns = source_columns) {
    annotated <- readRDS(file)
    idx <- match(colnames(obj), colnames(annotated))
    if (anyNA(idx)) 
        stop(sum(is.na(idx)), " cells are absent from ", file)
    for (i in seq_along(source_columns)) {
        obj[[target_columns[[i]]]] <- annotated[[source_columns[[i]], drop = TRUE]][idx]
    }
    obj
}

join_and_normalize <- function(obj) {
    DefaultAssay(obj) <- "RNA"
    if (inherits(obj[["RNA"]], "Assay5")) 
        obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
    NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
}

read_signature <- function(file, species = c("mouse", "human")) {
    species <- match.arg(species)
    genes <- trimws(readLines(file, warn = FALSE))
    unique(genes[nzchar(genes)])
}

signature_files <- list(myeloid_mouse = file.path(alex_root, "manuscript_signature_validation_20260820/tables/new_signature_71_mouse.txt"), 
    myeloid_human = file.path(alex_root, "manuscript_signature_validation_20260820/tables/new_signature_human_orthologs.txt"))

read_bm_signature <- function(column) {
    x <- read_excel(file.path(alex_root, "BM signatures.xlsx"), col_names = FALSE)
    unique(na.omit(trimws(as.character(x[[column]][-(1:2)]))))
}

signatures_mouse <- list(Myeloid_Inflammatory = read_signature(signature_files$myeloid_mouse, "mouse"), Hallmark_Inflammatory_Response = msigdbr::msigdbr(species = "Mus musculus", 
    collection = "H") %>% filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% pull(gene_symbol) %>% unique(), 
    Chambers_Inflammaging = read_bm_signature(4), Kovtonyuk_Myeloid_Bias = read_bm_signature(7))

score_object <- function(obj, signatures, prefix = "") {
    for (nm in names(signatures)) {
        genes <- intersect(signatures[[nm]], rownames(obj))
        if (length(genes) < 3) 
            stop("Too few genes found for ", nm)
        set.seed(1)
        obj <- AddModuleScore(obj, features = list(genes), name = paste0(prefix, nm), assay = "RNA", nbin = 24, 
            ctrl = 100, seed = 1, search = FALSE, slot = "data")
        avg_name <- paste0(prefix, nm, "_Average")
        obj[[avg_name]] <- Matrix::colMeans(GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE])
    }
    obj
}

split_violin_geom <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    group <- data[1, "group"]
    data2 <- transform(data, x = if (group%%2 == 1) 
        xminv
    else xmaxv)
    data2 <- data2[order(if (group%%2 == 1) 
        data2$y
    else -data2$y), ]
    data2 <- rbind(data2[1, ], data2, data2[nrow(data2), ], data2[1, ])
    data2[c(1, nrow(data2) - 1, nrow(data2)), "x"] <- round(data2[1, "x"])
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(data2, ...))
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., trim = TRUE, 
    scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
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
    p <- ggplot(d, aes(comparison, .data[[value]], fill = plot_group, group = side)) + geom_split_violin(trim = FALSE, 
        scale = "width", color = "black", linewidth = 0.25) + geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35, 
        color = "black") + annotate("text", x = 1, y = ymax + max(0.08 * yrange, 0.03), label = paste0("p = ", 
        format_p(pvalue)), size = 3.4) + scale_fill_manual(values = setNames(colors, levels)) + scale_y_continuous(expand = expansion(mult = c(0.03, 
        0.17))) + theme_classic(base_size = 11) + labs(title = title, x = NULL, y = ylab, fill = NULL)
    save_pdf(p, path, 7.5, 5.5)
    data.frame(group_1 = levels[[1]], n_1 = sum(d$plot_group == levels[[1]]), group_2 = levels[[2]], n_2 = sum(d$plot_group == 
        levels[[2]]), wilcox_p = pvalue)
}

run_deg <- function(obj, group_col, ident1, ident2, outdir, prefix) {
    Idents(obj) <- factor(obj[[group_col, drop = TRUE]])
    deg <- FindMarkers(obj, ident.1 = ident1, ident.2 = ident2, assay = "RNA", test.use = "wilcox", min.pct = 0.1, 
        logfc.threshold = 0, only.pos = FALSE)
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
    if (!requireNamespace("enrichR", quietly = TRUE)) 
        stop("The enrichR package is required")
    databases <- c("GO_Biological_Process_2025", "MSigDB_Hallmark_2020", "Reactome_2022")
    results <- list()
    for (direction in c("UP", "DOWN")) {
        genes <- unique(deg$gene[deg$direction == direction])
        if (!length(genes)) 
            next
        er <- enrichR::enrichr(genes, databases)
        for (db in names(er)) results[[paste(direction, db, sep = "_")]] <- er[[db]]
    }
    if (length(results)) 
        write.xlsx(results, file.path(outdir, paste0(prefix, "_enrichr.xlsx")), overwrite = TRUE)
    invisible(results)
}

plot_cell_umap <- function(obj, field, path, title) {
    p <- DimPlot(obj, reduction = "umap", group.by = field, label = TRUE, repel = TRUE, raster = FALSE, pt.size = 0.12) + 
        ggtitle(title) + theme_classic(base_size = 11)
    save_pdf(p, path, 11, 8)
}

run_post_treatment_brain <- function() {
    outdir <- file.path(output_root, "02_post_treatment_brain")
    ids <- c("Y_uPAR_F", "Y_UT_F", "O_UT_F", "Y_uPAR_M", "Y_UT_M", "O_uPAR_F", "O_uPAR_M", "O_UT_M")
    samples <- data.frame(sample = ids, path = file.path(alex_root, "dataset2", ids, "filtered_feature_bc_matrix.h5"), 
        age = ifelse(grepl("^Y_", ids), "Young", "Old"), treatment = ifelse(grepl("uPAR", ids), "uPAR", "UT"), 
        sex = ifelse(grepl("_F$", ids), "Female", "Male"))
    obj <- qc_filter(make_objects(samples), 5, 500, 20000, 250, 5000, outdir)
    obj <- integrate_rpca(obj, 1:6, 30, 0.3)
    obj <- attach_labels_from_rds(obj, file.path(alex_root, "dataset2/objects/d2_annotated_final.RDS"), c("cell_type_grouped", 
        "cell_type_grouped"), c("cell_type", "cell_type_grouped"))
    obj$analysis_group <- factor(paste(substr(obj$age, 1, 1), obj$treatment), levels = c("Y UT", "Y uPAR", "O UT", 
        "O uPAR"))
    obj <- score_object(join_and_normalize(obj), signatures_mouse)
    plot_cell_umap(obj, "cell_type_grouped", file.path(outdir, "cell_type_umap.pdf"), "Post-treatment brain")
    bam <- obj@meta.data %>% mutate(cell = rownames(.)) %>% filter(grepl("BAM", cell_type, ignore.case = TRUE))
    stats <- make_split_violin(bam, "Myeloid_Inflammatory_Average", "analysis_group", c("O UT", "O uPAR"), "Myeloid inflammatory signature in BAMs", 
        "Mean LogNormalized expression", file.path(outdir, "Figure_2G_myeloid_inflammatory_BAM.pdf"), c("#08306B", 
            "#FF4902"))
    write_csv(stats, file.path(outdir, "Figure_2G_statistics.csv"))
    deg <- run_deg(obj, "analysis_group", "O uPAR", "O UT", file.path(outdir, "DEG"), "Old_uPAR_vs_Old_UT")
    run_enrichr(deg, file.path(outdir, "DEG"), "Old_uPAR_vs_Old_UT")
    saveRDS(obj, file.path(outdir, "post_treatment_brain_processed.rds"), compress = FALSE)
}

run_post_treatment_brain()
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message("Complete: ", output_root)

