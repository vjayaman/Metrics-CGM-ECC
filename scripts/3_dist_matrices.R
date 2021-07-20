#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_distmatrices.txt", open="wt")
sink(msg, type="message")

source("scripts/ECC/ecc_opener.R")
# source("scripts/ECC/dist_functions.R")
source("scripts/all_distances.R")

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

for (k in 1:2) {
  dir.create(paste0("intermediate_data/TP", k), showWarnings = FALSE)
  dir.create(paste0("intermediate_data/TP", k, "/dists"), showWarnings = FALSE)
  
  outputMessages(paste0("\nCollecting and saving distances for groups at TP", k))
  dists <- paste0("intermediate_data/TP", k, "/dists/")
  parts <- sectionClusters(k, typing_data, m)
  
  # inter-cluster detour
  dr_clusters <- typing_data[[k]] %>% rownames_to_column("Strain") %>%
    as.data.table() %>% left_join(., m$dr_matches, by = "Strain") %>% 
    select(-Strain) %>% unique() %>% set_colnames(c("hx", "dr"))
  b1 <- dr_clusters$dr %>% unique()
  b2 <- t(combn(as.factor(b1), 2)) %>% as.data.table()
  b3 <- b2 %>% set_colnames(c("dr1", "dr2"))
  
  b4 <- left_join(dr_clusters, b3, by = c("dr" = "dr1")) %>% rename("hx1" = "hx", "dr1" = "dr")# %>% 
  b4 <- b4[!is.na(dr1) & !is.na(dr2)]
  b5 <- b4 %>% left_join(., dr_clusters, by = c("dr2" = "dr")) %>% rename("hx2" = "hx")
  # drs that are NOT in the same cluster:
  b6 <- b5[hx1 != hx2] %>% select(-hx1, -hx2) %>% unique()

  inter_dists <- c(b6$dr1, b6$dr2) %>% unique() %>% sort() %>% as_tibble() %>% 
    set_colnames("dr") %>% left_join(., m$assignments, by = "dr")
  
  temp_dists <- inter_dists %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
  geo_dists <- inter_dists %>% select(dr, Longitude, Latitude) %>% 
    distMatrix(., "geo", c("Longitude", "Latitude"))
  extremes <- list(maxt = max(temp_dists), mint = min(temp_dists), 
                   maxg = max(geo_dists), ming = min(geo_dists))
  
  collectDistances(TRUE, m$assignments, parts, fpaths = dists, extremes)
}

outputMessages("\nFinished saving pairwise distances, in groups.")

if (save_dist_matrix) {
  outputMessages("\nSaving full distance matrix (strains as rows and columns)...")
  
  read_dists <- "intermediate_data/TP2/dists/"
  grouped_dists <- list.files(read_dists, full.names = TRUE)
  
  temp_path <- paste0("results/temporal_distance/")
  dir.create(temp_path, showWarnings = FALSE)
  
  for (i in 1:length(grouped_dists)) {
    type <- c("temp", "Temp.Dist")
    dm <- readRDS(grouped_dists[i])[[type[1]]] %>% 
      as.data.table() %>% 
      rownames_to_column("dr") %>% melt.data.table(id.vars = "dr") %>% 
      set_colnames(c("dr1", "dr2", "Dist"))
    
    a1 <- parts$drs[,c("Strain","dr")] %>% right_join(., dm, by = c("dr" = "dr1")) %>% 
      rename(Strain1 = Strain) %>% select(-dr)
    rm(dm)
    
    a2 <- right_join(parts$drs[,c("Strain","dr")], a1, by = c("dr" = "dr2")) %>% 
      rename(Strain2 = Strain) %>% select(Strain1, Strain2, Dist) %>% 
      set_colnames(c("Strain1", "Strain2", type[2]))
    rm(a1)
    
    grouped_dists[i] %>% strsplit(.,"/") %>% unlist() %>% extract2(4) %>% 
      paste0(temp_path, .) %>% 
      saveRDS(a2, .)
    rm(a2)
  }
  
  geo_path <- paste0("results/geographical_distances/")
  dir.create(geo_path, showWarnings = FALSE)
  
  for (i in 1:length(grouped_dists)) {
    type <- c("geo", "Geog.Dist")
    # print(i)
    # print("a")
    dm <- readRDS(grouped_dists[i])[[type[1]]] %>% 
      as.data.table() %>% 
      rownames_to_column("dr") %>% melt.data.table(id.vars = "dr") %>% 
      set_colnames(c("dr1", "dr2", "Dist"))
    # print("b")
    a1 <- parts$drs[,c("Strain","dr")] %>% right_join(., dm, by = c("dr" = "dr1")) %>% 
      rename(Strain1 = Strain) %>% select(-dr)
    rm(dm)
    # print("c")
    a2 <- right_join(parts$drs[,c("Strain","dr")], a1, by = c("dr" = "dr2")) %>% 
      rename(Strain2 = Strain) %>% select(Strain1, Strain2, Dist) %>% 
      set_colnames(c("Strain1", "Strain2", type[2]))
    rm(a1)
    # print("d")
    grouped_dists[i] %>% strsplit(.,"/") %>% unlist() %>% extract2(4) %>% 
      paste0(geo_path, .) %>% saveRDS(a2, .)
    rm(a2)
  }
  outputMessages(paste0("\nSee results/geographical_distances/ or ", 
                        "results/temporal_distances/ \nfor the distances. Note they are ", 
                        "saved in groups of clusters, to help with in-memory constraints."))
}else {
  outputMessages(paste0("\nNot saving full distance matrix (strains as rows and columns); ", 
                        "saving reduced redundancy \npairwise distances. To save the full ", 
                        "distance matrix, add --distmat TRUE to input arguments. "))
}

# Average distances --------------------------------------------------------------------------
assert("Distances were collected and saved", file.exists("intermediate_data/dist_extremes.Rds"))

avg_dists <- lapply(1:2, function(tpx) {
  fpath1 <- paste0("intermediate_data/TP", tpx, "/dists/")
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

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), 
           " End of distances collection ", paste0(rep("-", 31), collapse = ""), "||\n")

# # Check the dist results for a few sampled at random -----------------------------------------
# a1 <- readRDS("results/geographical_distances/group1.Rds")
# strains <- a1$Strain1[1:100] %>% sort()
# a2 <- a1[Strain1 %in% strains & Strain2 %in% strains] %>% arrange(Strain1, Strain2)
# group1 <- m$strain_data %>% filter(Strain %in% strains) %>% 
#   select(Strain, Longitude, Latitude) %>% arrange(Strain)
# 
# dm <- group1 %>% column_to_rownames("Strain") %>% as.data.frame() %>% 
#   earth.dist(dist = TRUE) %>% as.matrix() %>% 
#   set_rownames(group1$Strain) %>% set_colnames(group1$Strain)
# 
# dm2 <- dm %>% 
#   as.data.frame() %>% rownames_to_column("Strain1") %>% as.data.table() %>% 
#   melt.data.table(., id.vars = "Strain1", variable.name = "Strain2", value.name = "Geog.Dist2") %>% 
#   arrange(Strain1, Strain2)
# 
# testers <- c("Fuzhou/JX2012/2020", "Foshan/20SF207/2020")
# tmp1 <- m$strain_data %>% filter(Strain %in% testers) %>% 
#   select(Strain, Longitude, Latitude) %>% as.data.frame() %>% 
#   column_to_rownames("Strain") %>% earth.dist(dist = TRUE) %>% as.matrix()
# 
# b1 <- m$dr_matches %>% filter(Strain %in% testers) %>% pull(dr)
# m$assignments[dr %in% b1]
# 
# # Note: with this method, we don't have per-strain inter-cluster distances, 
# # we infer the intracluster distances using data representatives
# # So we can check our process but checking the end results, not the intermediate
# # results along the way
# # Note: should randomly sample and check every file saved to the "results/" directory