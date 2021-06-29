#! /usr/bin/env Rscript

msg <- file("logs/logfile_epiquant.txt", open="wt")
sink(msg, type="message")

source("scripts/ecc_opener.R")
assert("Distances were collected and saved", file.exists("results/TP1/dists/group1.Rds"))

# Average distances --------------------------------------------------------------------------
tpx <- 2

fpath1 <- paste0("results/TP", tpx, "/dists/")
tpx_dists <- list.files(fpath1) %>% paste0(fpath1, .)

groups <- sectionClusters(tpx, typing_data, m)
groups$drs %<>% set_colnames(c("Strain", "Th", "dr"))

dr_td1 <- typing_data[[tpx]] %>% rownames_to_column("Strain") %>% as_tibble() %>%
  left_join(., m$dr_matches, by = "Strain") %>%
  mutate(across(dr, as.character)) %>% select(-Strain)  

# Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
cx <- colnames(dr_td1)[1]
tallied_reps <- dr_td1 %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
g_cuts <- left_join(dr_td1, tallied_reps, by = intersect(colnames(tallied_reps), colnames(dr_td1))) %>%
  unique() %>% mutate(across(dr, as.character))

temp_avg_dists <- lapply(1:length(tpx_dists), function(j) {
  group_j <- pull(groups$results[[j]], 1)
  dm_temp <- readRDS(tpx_dists[j])[["temp"]]
  cluster_i <- groups$drs[Th %in% group_j]
  unidrs <- unique(cluster_i$dr)
  dists_i <- dm_temp[unidrs, unidrs]

  # outputMessages(paste0("      Calculating average (not transformed) distances for timepoint ", i))
  g_cuts %>% set_colnames(c("Th", "dr", "n")) %>% 
    filter(Th %in% group_j) %>% 
    set_colnames(colnames(g_cuts)) %>% 
    avgDists(., dists_i, "Temp.Dist", paste0("TP", tpx, "_", colnames(g_cuts)[1])) %>% return()
}) %>% bind_rows()


geo_avg_dists <- lapply(1:length(tpx_dists), function(j) {
  group_j <- pull(groups$results[[j]], 1)
  dm_geo <- readRDS(tpx_dists[j])[["geo"]]
  cluster_i <- groups$drs[Th %in% group_j]
  unidrs <- unique(cluster_i$dr)
  dists_i <- dm_geo[unidrs, unidrs]
  
  # outputMessages(paste0("      Calculating average (not transformed) distances for timepoint ", i))
  g_cuts %>% set_colnames(c("Th", "dr", "n")) %>% 
    filter(Th %in% group_j) %>% 
    set_colnames(colnames(g_cuts)) %>% 
    avgDists(., dists_i, "Geog.Dist", paste0("TP", tpx, "_", colnames(g_cuts)[1])) %>% return()
}) %>% bind_rows()

inner_join(geo_avg_dists, temp_avg_dists) %>% return()
