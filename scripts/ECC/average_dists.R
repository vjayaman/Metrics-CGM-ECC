#! /usr/bin/env Rscript

# msg <- file("logs/logfile_epiquant.txt", open="wt")
# sink(msg, type="message")
# 
# source("scripts/ecc_opener.R")
assert("Distances were collected and saved", file.exists("results/TP1/dists/group1.Rds"))

countDataReps <- function(typing_data, tpx, dr_matches) {
  dr_td1 <- typing_data[[tpx]] %>% rownames_to_column("Strain") %>% as_tibble() %>%
    left_join(., dr_matches, by = "Strain") %>%
    mutate(across(dr, as.character)) %>% select(-Strain)  
  
  # Counting data representatives (so we know how much to 
  # multiply each ECC value by to represent all strains)
  cx <- colnames(dr_td1)[1]
  tallied_reps <- dr_td1 %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
  g_cuts <- left_join(dr_td1, tallied_reps, 
                      by = intersect(colnames(tallied_reps), colnames(dr_td1))) %>%
    unique() %>% mutate(across(dr, as.character))
  return(g_cuts)
}

avgsFromDM <- function(tpx_dists, groups, g_cuts, type, cname, tpx) {
  lapply(1:length(tpx_dists), function(j) {
    group_j <- pull(groups$results[[j]], 1)
    dm_type <- readRDS(tpx_dists[j])[[type]]
    cluster_i <- groups$drs[Th %in% group_j]
    unidrs <- unique(cluster_i$dr)
    dists_i <- dm_type[unidrs, unidrs]
    
    # outputMessages(paste0("      Calculating average (not transformed) distances for timepoint ", i))
    g_cuts %>% set_colnames(c("Th", "dr", "n")) %>% 
      filter(Th %in% group_j) %>% 
      set_colnames(colnames(g_cuts)) %>% 
      avgDists(., dists_i, cname, paste0("TP", tpx, "_", colnames(g_cuts)[1])) %>% return()
  }) %>% bind_rows() %>% return()
}
# Average distances --------------------------------------------------------------------------

avg_dists <- lapply(1:2, function(tpx) {
  fpath1 <- paste0("results/TP", tpx, "/dists/")
  tpx_dists <- list.files(fpath1) %>% paste0(fpath1, .)
  
  groups <- sectionClusters(tpx, typing_data, m)
  groups$drs %<>% set_colnames(c("Strain", "Th", "dr"))
  g_cuts <- countDataReps(typing_data, tpx, m$dr_matches)
  
  temp_dists <- avgsFromDM(tpx_dists, groups, g_cuts, "temp", "Temp.Dist", tpx)
  geo_dists <- avgsFromDM(tpx_dists, groups, g_cuts, "geo", "Geog.Dist", tpx)
  
  inner_join(temp_dists, geo_dists) %>% return()
}) %>% set_names(c("TP1", "TP2"))

tp1_avg_dists <- tp1$proc %>% select(-TP1) %>% left_join(., avg_dists[["TP1"]])
tp2_avg_dists <- tp2$proc %>% select(-TP2) %>% left_join(., avg_dists[["TP2"]])

all_avg_dists <- right_join(tp1_avg_dists, tp2_avg_dists, by = "Strain")
