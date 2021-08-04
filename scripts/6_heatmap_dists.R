#! /usr/bin/env Rscript

msg <- file("logs/logfile_heatmaps.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "fossil", "tidyr", "plyr", "dplyr", "readr", 
          "testit", "tibble", "reshape2", "RColorBrewer", "gplots", "data.table", "R6")
y <- suppressWarnings(
  suppressPackageStartupMessages(
    lapply(libs, require, character.only = TRUE)))

source("scripts/ECC/classes_ecc.R")
source("scripts/ECC/dist_functions.R")
source("scripts/ECC/ecc_functions.R")
source("report_specific/epi-helper-no-source.R")

getCoeff <- function(x, y) {
  x %>% strsplit("-") %>% unlist() %>% extract2(y) %>% as.double()
}

dir.create("report_specific/heatmaps/", showWarnings = FALSE)

fnames <- list.files("intermediate_data/TP2/dists/", full.names = TRUE)
distfiles <- lapply(fnames, function(f) readRDS(f))
extremes <- readRDS("intermediate_data/dist_extremes.Rds")

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-a", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  make_option(c("-c", "--CGMs"), metavar = "file", 
              default = "results/CGM_strain_results.tsv", help = "CGM result file"),
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))

arg <- parse_args(OptionParser(option_list=option_list))

# Extract threshold of interest and coefficients from form inputs ----------------------------------------------
params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ")

test_params <- c("Region of interest", "Country of interest", "Has defined lineage information", 
                 "Has defined date information (day, month, and year)", "Has province-level data", 
                 "Province of interest", "Threshold of interest", "Is in a non-singleton cluster (at TP1)", 
                 "Is in a non-singleton cluster (at TP2)", "Filtering by date", "Column names", 
                 "Source-temporal-geographic coefficents", "Generate heatmaps for top __ largest clusters")

assert("Input parameters are correctly labelled", identical(sapply(params, '[[', 1), test_params))

params %<>% set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
                        "th","nsTP1","nsTP2", "temp_win","cnames","coeffs", "numcl"))

number_of_clusters <- as.double(params$numcl[2])

if (number_of_clusters > 0) {
  source("scripts/Misc/generate_heatmaps.R")
}
