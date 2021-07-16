#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_distmatrices.txt", open="wt")
sink(msg, type="message")

source("scripts/ECC/ecc_opener.R")
source("scripts/ECC/dist_functions.R")

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
  extremes <- "intermediate_data/"
  parts <- sectionClusters(k, typing_data, m)
  collectDistances(TRUE, m$assignments, parts, fpaths = list(dists, extremes))
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
      paste0(geo_path, .) %>% 
      saveRDS(a2, .)
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
           " End of distances collection ", 
           paste0(rep("-", 31), collapse = ""), "||\nStarted process at: ", Sys.time()))
