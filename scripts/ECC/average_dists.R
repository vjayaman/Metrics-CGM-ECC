#! /usr/bin/env Rscript

# msg <- file("logs/logfile_epiquant.txt", open="wt")
# sink(msg, type="message")
# 
# source("scripts/ecc_opener.R")
assert("Distances were collected and saved", file.exists("results/TP1/dists/group1.Rds"))

source("scripts/ECC/functions/dist_functions.R")
# Average distances --------------------------------------------------------------------------

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
