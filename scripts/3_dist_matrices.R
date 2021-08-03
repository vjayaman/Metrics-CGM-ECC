#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_generatedists.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source))

option_list <- list(
  make_option(c("-m", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-a", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

m <- read_tsv(params$strains) %>% processedStrains()

cat(paste0("\n||", paste0(rep("-", 23), collapse = ""), 
           " Generating non-redundant pairwise distances ", 
           paste0(rep("-", 23), collapse = ""), "||\nStarted process at: ", Sys.time()))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

# Note: dr stands for data representative
# in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
# cat(paste0("\n\nStep 1:"))
# cat(paste0("\n   Note that the source coefficent is always 0 in this version"))
# outputMessages("   Removing redundancy (comparing date, lat-long, etc. - not every pair of strains)")
# outputMessages("   Identifying which strains match which non-redundant 'data representatives'")

# COLLECT dist matrices using TP2 clusters ---------------------------------------------------------

for (k in c(2,1)) {
  dir.create(paste0("intermediate_data/TP", k), showWarnings = FALSE)
  dir.create(paste0("intermediate_data/TP", k, "/dists"), showWarnings = FALSE)
  
  outputMessages(paste0("\nCollecting and saving distances for groups at TP", k))
  dists <- paste0("intermediate_data/TP", k, "/dists/")
  parts <- sectionClusters(k, typing_data, m)
  
  if (k == 2) {
    geo_dists <- m$assignments %>% select(Longitude, Latitude) %>% unique() %>% 
      rownames_to_column("id") %>% distMatrix(., "geo", c("Longitude", "Latitude"))
    temp_dists <- m$assignments %>% select(Date) %>% unique() %>% 
      rownames_to_column("id") %>% distMatrix(., "temp", "Date")
    
    extremes <- list(maxt = max(temp_dists), mint = min(temp_dists), 
                     maxg = max(geo_dists), ming = min(geo_dists))
    saveRDS(extremes, "intermediate_data/dist_extremes.Rds")
  }
  
  outputMessages("\nGenerating intra-cluster distances:")
  # assignments <- m$assignments; fpaths <- dists
  collectDistances(m$assignments, parts, fpaths = dists)
  # collectDistances(TRUE, m$assignments, parts, fpaths = dists, extremes)
}

outputMessages("\nFinished saving non-redundant pairwise distances, in groups.")

# Average distances --------------------------------------------------------------------------
assert("Distances were collected and saved", file.exists("intermediate_data/dist_extremes.Rds"))

outputMessages("\nCollecting average distances ...")
avg_dists <- lapply(1:2, function(tpx) {
  tpx_dists <- paste0("intermediate_data/TP", tpx, "/dists/") %>% 
    list.files(., pattern = "group", full.names = TRUE)
  
  groups <- sectionClusters(tpx, typing_data, m)
  groups$drs %<>% set_colnames(c("Strain", "Th", "dr"))
  
  dr_td1 <- typing_data[[tpx]] %>% 
    rownames_to_column("Strain") %>% as_tibble() %>%
    left_join(., m$dr_matches, by = "Strain") %>%
    mutate(across(dr, as.character)) %>% select(-Strain)  
  
  g_cuts <- countDataReps(dr_td1)
  
  temp_dists <- avgsFromDM(tpx_dists, groups, g_cuts, "temp", "Temp.Dist", tpx)
  geo_dists <- avgsFromDM(tpx_dists, groups, g_cuts, "geo", "Geog.Dist", tpx)
  
  inner_join(temp_dists, geo_dists, 
             by = intersect(colnames(temp_dists), colnames(geo_dists))) %>% return()
}) %>% set_names(c("TP1", "TP2"))

tp1_avg_dists <- tp1$proc %>% select(-TP1) %>% 
  left_join(., avg_dists[["TP1"]], by = intersect(colnames(.), colnames(avg_dists[["TP1"]])))
tp2_avg_dists <- tp2$proc %>% select(-TP2) %>% 
  left_join(., avg_dists[["TP2"]], by = intersect(colnames(.), colnames(avg_dists[["TP2"]])))

all_avg_dists <- right_join(tp1_avg_dists, tp2_avg_dists, by = "Strain")
saveRDS(all_avg_dists, "intermediate_data/average_dists.Rds")

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), 
           " End of distances collection ", paste0(rep("-", 31), collapse = ""), "||\n"))
