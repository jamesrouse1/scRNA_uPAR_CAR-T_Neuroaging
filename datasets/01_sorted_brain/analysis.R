#!/usr/bin/env Rscript

script_arg <- commandArgs()[grepl("^--file=", commandArgs())]
script_file <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
Sys.setenv(ALEX_MANUSCRIPT_PROJECT_ROOT = project_root)
source(file.path(project_root, "R", "manuscript_functions.R"))
run_manuscript_dataset("sorted_brain")
