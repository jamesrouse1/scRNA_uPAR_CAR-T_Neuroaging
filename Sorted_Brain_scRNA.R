#!/usr/bin/env Rscript

script_arg <- commandArgs()[grepl("^--file=", commandArgs())]
project_root <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
Sys.setenv(ALEX_MANUSCRIPT_PROJECT_ROOT = project_root)
source(file.path(project_root, "Analysis_Functions.R"))
run_manuscript_dataset("sorted_brain")
