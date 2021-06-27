#! /usr/bin/env Rscript

msg <- file("logs/logfile_epiquant.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- paste0("scripts/ECC") %>% 
  list.files(., full.names = TRUE) %>% 
  grep("dist_matrices.R", ., invert = TRUE, value = TRUE)
invisible(sapply(files, source))

# Title: "EpiQuant - Salmonella Enteritidis Project (2019-2020)"
# Authors of original work and initial modifications: Ben Hetman, Elissa Giang, Dillon Barker
# Responsible for changes made during and after merge with CGM process: Vasena Jayamanna

option_list <- list(
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"),
  make_option(c("-t", "--trio"), metavar = "character", default = "010",
              help = "temporal, geographic coefficients"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ", 
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$trio %>% strsplit(., "-") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

m <- read_tsv(params$strains) %>% processedStrains()

# Note: dr stands for data representative
# in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
# cat(paste0("\n\nStep 1:"))
# cat(paste0("\n   Note that the source coefficent is always 0 in this version"))
# outputMessages("   Removing redundancy (comparing date, lat-long, etc. - not every pair of strains)")
# outputMessages("   Identifying which strains match which non-redundant 'data representatives'")

# TIMEPOINT 2 ANALYSIS -----------------------------------------------------------------------
c1 <- unlist(strsplit(combos[1], split = "")) %>% as.numeric()

for (k in 1:2) {
  outputMessages(paste0("\nCollecting and saving ECCs for groups of clusters at TP", k))
  extremes <- readRDS(paste0("results/tmp/TP", k, "-dists/extremes.Rds"))
  
  sectionClusters(k, typing_data, m) %>% collectECCs(k, m, ., extremes, c1)
  
  outputMessages(paste0("\nMerging ECC files at TP", k))
  dir_k <- paste0("results/tmp/TP", k, "-dists/")
  fnames <- grep("eccs", list.files(dir_k), value = TRUE)
  
  tpk <- lapply(fnames, function(f) {
    ecc_j <- readRDS(paste0(dir_k, f))
    file.remove(paste0(dir_k, f))
    return(ecc_j)
  }) %>% bind_rows()
  
  saveRDS(tpk, paste0("results/tmp/TP", k, "-dists/all_eccs.Rds"))
}

# ---------------------------------------------------------------------------------------------
a1 <- readRDS("results/tmp/TP1-dists/all_eccs.Rds")
b1 <- readRDS("results/tmp/TP2-dists/all_eccs.Rds")

tp1_eccs <- tp1$proc[,c("Strain",colnames(a1)[1])] %>% left_join(., a1)
tp2_eccs <- tp2$proc[,c("Strain",colnames(b1)[1])] %>% left_join(., b1)

all_eccs <- right_join(tp1_eccs, tp2_eccs, by = "Strain")

cat(paste0("\n\nStep ", length(combos) + 2, ":"))
outputMessages("   Merging collected ECCs ...\n")
write.table(all_eccs, file = "results/ECCs.tsv", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ",
           paste0(rep("-", 31), collapse = ""), "||\n"))