#! /usr/bin/env Rscript

msg <- file("logs/logfile_epiquant1.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source))

assert("Distances were collected and saved", file.exists("intermediate_data/dist_extremes.Rds"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ",
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))

option_list <- list(
  make_option(c("-m", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-t", "--trio"), metavar = "character", default = "0-1-0",
              help = "source, temporal, geographic coefficients"))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$trio %>% strsplit(., ",") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp2$height_list %>% set_names("2")

m <- read_tsv(params$strains) %>% processedStrains()

extremes <- readRDS("intermediate_data/dist_extremes.Rds")

dfx <- expand.grid(x = combos, k = 2) %>% as.data.frame()

for (k in unique(dfx$k)) {
  dir.create(paste0("intermediate_data/TP", k, "/ecc_groups/"), showWarnings = FALSE)
}

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)
for (i in 1:nrow(dfx)) {
  cat(paste0("\n\nStep ", i, " / ", nrow(dfx) + 2, ":"))
  c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "-") %>% 
    unlist() %>% as.numeric()
  
  paste0("intermediate_data/TP", dfx$k[i], "/ecc_groups/", as.character(dfx$x[i])) %>% 
    dir.create(., showWarnings = FALSE)
  
  outputMessages(paste0("\nCollecting and saving ECCs for groups of clusters at TP", 
                        dfx$k[i], ", for ECC coefficient triple ", as.character(dfx$x[i])))
  
  k <- dfx$k[i]
  parts <- sectionClusters(k, typing_data, m)
  read_from <- paste0("intermediate_data/TP", k, "/dists/")
  save_to <- paste0("intermediate_data/TP", k, "/ecc_groups/", as.character(dfx$x[i]), "/")
  
  collectECCs(k, m, parts, extremes, c1, read_from, save_to)
  
  c1 %>% paste0(., collapse = ", ") %>% 
    paste0("\nMerging ECC files at TP", dfx$k[i], ", for ECC parameters ", .) %>% 
    outputMessages()
}

cat(paste0("\n\nStep ", nrow(dfx) + 1, " / ", nrow(dfx) + 2, ":"))
for (i in 1:nrow(dfx)) {
  fnames <- list.files(paste0("intermediate_data/TP", dfx$k[i], 
                              "/ecc_groups/", dfx$x[i], "/"), full.names = TRUE)
  
  tpk <- lapply(fnames, function(f) {readRDS(f)}) %>% bind_rows()
  
  saveRDS(tpk, paste0("intermediate_data/TP", dfx$k[i], "/", dfx$x[i], "-eccs.Rds"))
}

# Generating ECC results file ----------------------------------------------------
# source("scripts/5_merging_eccs.R")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ",
           paste0(rep("-", 31), collapse = ""), "||\n"))
