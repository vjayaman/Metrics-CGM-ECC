#! /usr/bin/env Rscript

# The number of clusters to generate pairwise distances for
# Will collect the largest ones, looking at the full data set (last time point, 
# when all strains are present)
# If this method seems alright, I'll update the form.html input file to include this parameter
num_cl <- 20

# General setup - calling libraries, etc. --------------------------------------------------------------
msg <- file("logs/pairwise_strain_distances.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "fossil", "tidyr", "plyr", "dplyr", "readr", 
          "testit", "tibble", "reshape2", "RColorBrewer", "gplots", "data.table", "R6")
y <- suppressWarnings(
  suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE)))

cat(paste0("\n||", paste0(rep("-", 21), collapse = ""), " Starting pairwise (strain) distance collection ", 
           paste0(rep("-", 22), collapse = ""), "||\n"))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

# Adjust input parameters here: ------------------------------------------------------------------------
# This section collects user-specified inputs
option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", 
              help = "TP2 cluster assignments"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))

# Reading in required functions, checking for pre-processed results (e.g. distances) -------------------
source("scripts/ECC/classes_ecc.R")
source("scripts/ECC/dist_functions.R")
source("scripts/ECC/ecc_functions.R")
source("report_specific/epi-helper-no-source.R")
cluster_mapping <- readRDS("inputs/processed/allTP2.Rds")

dir.create("report_specific/epitables", showWarnings = FALSE, recursive = TRUE)

fnames <- list.files("intermediate_data/TPN/dists/", full.names = TRUE)
assert("Distances were collected, script '3_dist_matrices.R' was run correctly", length(fnames) >= 1)
distfiles <- lapply(fnames, function(f) readRDS(f))

assert("Max and min of distances collected", file.exists("intermediate_data/TPN/extreme_dists.Rds"))
extremes <- readRDS("intermediate_data/TPN/extreme_dists.Rds")

arg <- parse_args(OptionParser(option_list=option_list))

params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

# Reading in the threshold column (to get appropriate clusters)
hx <- strsplit(as.character(params$th[2]), split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

# Reading the CGM results (weekly/monthly/multiset)
results_files <- list.files("results/", full.names = TRUE) %>% grep(params$int_type[2], ., value = TRUE)
cgms <- grep("CGM", results_files, value = TRUE) %>% readRDS()

# NON-REDUNDANT METHOD
last_ivl <- unique(cgms$interval) %>% last()

# Collecting the top x clusters by size (after all strains are in the data, last time point)
# x is given by the user in the form_input.txt file (the "Generate heatmaps" line)
size_details <- cgms[interval == last_ivl] %>% 
  arrange(-tp2_cl_size) %>% select(first_tp2_flag, tp2_cl_size) %>% unique()
colnames(size_details)[which(colnames(size_details) == "tp2_cl_size")] <- "TP2_cluster_size"

top_clusters <- size_details %>% arrange(-TP2_cluster_size) %>% 
  slice(1:num_cl) %>% pull(first_tp2_flag)

# Matching these clusters to their original cluster names, for saving purposes
clusters <- data.table(
  new_id = top_clusters, 
  new_h = top_clusters %>% strsplit(., split = "_") %>% sapply(., '[[', 2) %>% 
    gsub("h", "", .) %>% as.double(), 
  new_cl = top_clusters %>% strsplit(., split = "_") %>% sapply(., '[[', 3) %>% 
    gsub("c", "", .) %>% as.integer()
) %>% inner_join(cluster_mapping$lookup_table, ., by = c("new_h", "new_cl")) %>% 
  select(new_id, new_cl, old_h, old_cl) %>% 
  set_colnames(c("new_id", "new_cl", "original_h", "original_cl")) %>% 
  inner_join(., size_details, by = c("new_id" = "first_tp2_flag")) %>% 
  arrange(-TP2_cluster_size)

# Reading in the last timepoint data, as well as the strain metadata and the 
# data representative assignments - used for processing distance files from intermediate_data/
tpn <- tp2 <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- list(tp2$height_list)

m <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()

assignments <- typing_data[[1]] %>% as.data.frame() %>% set_colnames("tp_cl") %>% 
  rownames_to_column("Strain") %>% as.data.table() %>% 
  filter(tp_cl %in% clusters$new_cl)

cl_drs <- lapply(1:nrow(clusters), function(i) {
  x1 <- assignments[tp_cl %in% clusters$new_cl[i]] %>% pull(Strain)
  m$dr_matches %>% filter(Strain %in% x1) %>% pull(dr) %>% unique()
}) %>% set_names(clusters$new_cl)

# Collectiing epitable data for each of the clusters, with column names: 
#       ||---------------------------------------------||
#       || Strain.1 | Strain.2 | Temp.Dist | Geog.Dist ||
#       ||---------------------------------------------||

outputDetails(paste0("Collecting pairwise distances (from non-redundant set), for each of ", 
                     nrow(clusters), " clusters"), newcat = TRUE)
for (j in 1:nrow(clusters)) {
  clx <- clusters$new_cl[j]
  outputDetails(paste0("Cluster ", gsub("TP2", "TPN", clusters$new_id[j]), " = ", 
                       clusters$original_cl[j], " (of size ", clusters$TP2_cluster_size[j], "), "))
  
  rawdists <- lapply(1:length(distfiles), function(i) {
    drs <- cl_drs[[as.character(clx)]]
    
    tdm <- distfiles[[i]]$temp
    tdm_subset <- tdm[rownames(tdm) %in% drs, colnames(tdm) %in% drs]
    
    gdm <- distfiles[[i]]$geo  
    gdm_subset <- gdm[rownames(gdm) %in% drs, colnames(gdm) %in% drs]
    
    if (length(tdm_subset) == 0 & length(gdm_subset) == 0) {
      data.table("dr1" = NA, "dr2" = NA, "Temp.Dist" = NA, "Geog.Dist" = NA) 
    }else {
      ctdm <- tdm_subset %>% 
        # COMMENT OUT THIS PART IF YOU DO NOT WANT TRANSFORMED (SCALED) DATA
        # This is done using a function from scripts/ECC/ecc_functions.R
        transformData2(., "temp", extremes$mint, extremes$maxt) %>%
        formatData(., c("dr1","dr2","Temp.Dist"))
      
      cgdm <- gdm_subset %>% 
        # COMMENT OUT THIS PART IF YOU DO NOT WANT TRANSFORMED (SCALED) DATA
        # This is done using a function from scripts/ECC/ecc_functions.R
        transformData2(., "geo", extremes$ming, extremes$maxg) %>%
        formatData(., c("dr1","dr2","Geog.Dist"))
      
      merge.data.table(ctdm, cgdm)
    }
  }) %>% bind_rows()
  
  strains_in_cluster <- assignments[tp_cl %in% clx] %>% pull(Strain)
  strain_drs <- m$dr_matches %>% filter(Strain %in% strains_in_cluster) %>% as.data.table()
  cluster_raws <- rawdists[!is.na(dr1)] %>% unique()
  
  first_strain <- inner_join(strain_drs, cluster_raws, by = c("dr" = "dr1")) %>% select(-dr)
  
  rm(strains_in_cluster); rm(cluster_raws); rm(rawdists); gc()
  
  colnames(first_strain)[which(colnames(first_strain)=="Strain")] <- "Strain.1"
  
  pw_strains <- strain_drs %>% set_colnames(c("Strain", "dr2")) %>% 
    merge.data.table(first_strain, ., by = "dr2", allow.cartesian = TRUE)
  
  colnames(pw_strains)[which(colnames(pw_strains)=="Strain")] <- "Strain.2"
  
  pw_strains <- pw_strains %>% select(Strain.1, Strain.2, Temp.Dist, Geog.Dist)
  
  outputDetails("saving ...", TRUE)
  saveRDS(pw_strains, file = paste0("report_specific/epitables/", 
                                    clusters$new_id[j], "-", clusters$original_cl[j], ".Rds"))
}

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
timeTaken(pt = "pairwise distances collection", stopwatch) %>% outputDetails(., newcat = TRUE)

cat(paste0("\n||", paste0(rep("-", 22), collapse = ""), " End of pairwise (strain) distance collection ", 
           paste0(rep("-", 23), collapse = ""), "||\n"))
