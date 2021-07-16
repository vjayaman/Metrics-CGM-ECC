#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_distmatrices.txt", open="wt")
sink(msg, type="message")

source("scripts/ECC/ecc_opener.R")
source("scripts/ECC/dist_functions.R")

# Note: dr stands for data representative
# in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
# cat(paste0("\n\nStep 1:"))
# cat(paste0("\n   Note that the source coefficent is always 0 in this version"))
# outputMessages("   Removing redundancy (comparing date, lat-long, etc. - not every pair of strains)")
# outputMessages("   Identifying which strains match which non-redundant 'data representatives'")

# COLLECT dist matrices using TP2 clusters ---------------------------------------------------------

for (k in 1:2) {
  dir.create(paste0("results/TP", k), showWarnings = FALSE)
  dir.create(paste0("results/TP", k, "/dists"), showWarnings = FALSE)
  
  outputMessages(paste0("\nCollecting and saving distances for groups at TP", k))
  dists <- paste0("results/TP", k, "/dists/")
  extremes <- "results/"
  parts <- sectionClusters(k, typing_data, m)
  collectDistances(TRUE, m$assignments, parts, fpaths = list(dists, extremes))
}

outputMessages("\nFinished saving distance matrices.")

# Average distances --------------------------------------------------------------------------
assert("Distances were collected and saved", file.exists("results/TP1/dists/group1.Rds"))

avg_dists <- lapply(1:2, function(tpx) {
  fpath1 <- paste0("results/TP", tpx, "/dists/")
  tpx_dists <- list.files(fpath1) %>% paste0(fpath1, .)
  
  groups <- sectionClusters(tpx, typing_data, m)
  groups$drs %<>% set_colnames(c("Strain", "Th", "dr"))
  
  dr_td1 <- typing_data[[tpx]] %>% 
    rownames_to_column("Strain") %>% as_tibble() %>%
    left_join(., m$dr_matches, by = "Strain") %>%
    mutate(across(dr, as.character)) %>% select(-Strain)  
  
  g_cuts <- countDataReps(dr_td1)
  
  temp_dists <- avgsFromDM(tpx_dists, groups, g_cuts, "temp", "Temp.Dist", tpx)
  geo_dists <- avgsFromDM(tpx_dists, groups, g_cuts, "geo", "Geog.Dist", tpx)
  
  inner_join(temp_dists, geo_dists) %>% return()
}) %>% set_names(c("TP1", "TP2"))

tp1_avg_dists <- tp1$proc %>% select(-TP1) %>% left_join(., avg_dists[["TP1"]])
tp2_avg_dists <- tp2$proc %>% select(-TP2) %>% left_join(., avg_dists[["TP2"]])

all_avg_dists <- right_join(tp1_avg_dists, tp2_avg_dists, by = "Strain")

saveRDS(all_avg_dists, "results/average_dists.Rds")
