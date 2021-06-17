#! /usr/bin/env Rscript

# msg <- file("logs/logfile_epiquant.txt", open="wt")
# sink(msg, type="message")

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
extremes <- readRDS("results/tmp/TP2-dists/extremes.Rds")

c1 <- unlist(strsplit(combos[1], split = "")) %>% as.numeric()

for (k in 1:2) {
  outputMessages(paste0("\nCollecting and saving ECCs for groups of clusters at TP", k))
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

# temporal_dists <- strain_data %>% 
#   select(dr, Date) %>%
#   distMatrix(., "temp", "Date") %>% 
#   transformData(., "temp") %>%
#   formatData(., c("dr1", "dr2", "Temp.Dist"))

# # Distance matrix, then transformation step, then formatting data
# tr_temp <- m$assignments %>% 
#   select(dr, Date) %>%
#   distMatrix(., "temp", "Date") %>% 
#   transformData(., "temp") %>%
#   formatData(., c("dr1", "dr2", "Temp.Dist"))
# 
# tr_geo <- m$assignments %>% select(dr, Latitude, Longitude) %>%
#   distMatrix(., "geo", c("Latitude", "Longitude")) %>%
#   pairwiseDists(., "geo", c("dr1", "dr2", "Geog.Dist"))
# 
# tr_dists <- merge.data.table(tr_temp, tr_geo)
# original_eccs <- tr_dists %>%
#   epiCollection(m$strain_data, tau, gamma, typing_data, ., dm_temp, dm_geo, m$dr_matches)
# # saveRDS(original_eccs, "results/tmp/original_eccs.Rds")
# 
# tp1_eccs <- readRDS("results/tmp/TP1-dists/all_eccs.Rds")
# colnames(tp1_eccs)[grepl("ECC", colnames(tp1_eccs))] %<>% paste0("New_", .)
# 
# inner_join(original_eccs[[1]], tp1_eccs) %>% 
#   mutate(dif = abs(New_TP1_T0_ECC.0.1.0 - TP1_T0_ECC.0.1.0)) %>% 
#   filter(!is.na(dif)) %>% 
#   filter(dif > 0)
# 
# # original_eccs <- readRDS("results/tmp/original_eccs.Rds")[[k]]
# colnames(new_eccs)[grepl("ECC", colnames(new_eccs))] %<>% paste0("New_", .)
# eccs <- inner_join(new_eccs, original_eccs)
# eccs$dif <- abs(pull(eccs, 3) - pull(eccs, 4))
# assert("Original ECCs match the new ones (cluster sectioning)", nrow(eccs[!is.na(dif)][dif > 1e-10]) == 0)
#   
# eccs[!is.na(dif)][dif > 1e-10] %>% as.data.table()



# # a1 <- readRDS("results/collected_eccs.Rds")
# collected_eccs <- lapply(1:length(combos), function(j) {
#   c1 <- unlist(strsplit(combos[j], split = "")) %>% as.numeric()
#   cat(paste0("\n\nStep ", j + 1, ":"))
#   tau <- c1[2]
#   gamma <- c1[3]
#   epiCollectionByCluster(m$strain_data, tau, gamma, transformed_dists, 2, cluster_x)
#   # epiCollection(m$strain_data, tau, gamma, typing_data, 
#   #               transformed_dists, dm_temp, dm_geo, m$dr_matches, j)
# })
# all_equal(a1, original_eccs)
# 
# all_equal(a1[[1]][[2]], collected_eccs[[1]])
# all_equal(a1[[2]][[2]], collected_eccs[[2]])

# cat(paste0("\n\nStep ", length(combos) + 2, ":"))
# outputMessages("   Merging collected ECCs ...\n")
# full_set <- mergeECCs(collected_eccs, 1, tp1$proc) %>%
#   merge.data.table(., mergeECCs(collected_eccs, 2, tp2$proc), by = "Strain", all.y = TRUE) %>%
#   mutate(TP1 = ifelse(is.na(TP1), 0, TP1))
# 
# write.table(full_set, file = "results/ECCs.tsv", col.names = TRUE, 
#             row.names = FALSE, quote = FALSE, sep = "\t")
# 
# stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
# cat(timeTaken(pt = "ECC data collection", stopwatch))
# 
# cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ", 
#            paste0(rep("-", 31), collapse = ""), "||\n"))
