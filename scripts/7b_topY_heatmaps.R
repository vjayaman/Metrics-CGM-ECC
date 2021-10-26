#! /usr/bin/env Rscript

# Will not generate heatmaps for cluster of size > max_cluster_size
max_cluster_size <- 2000

# General setup - calling libraries, etc. --------------------------------------------------------------
msg <- file("logs/heatmaps.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "fossil", "tidyr", "plyr", "dplyr", "readr", 
          "testit", "tibble", "reshape2", "RColorBrewer", "gplots", "data.table", "R6")
y <- suppressWarnings(
  suppressPackageStartupMessages(
    lapply(libs, require, character.only = TRUE)))

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), " Starting heatmap generation ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))

# Adjust input parameters here: ------------------------------------------------------------------------
# This section collects user-specified inputs
option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))

# Reading in required functions, checking for pre-processed results (e.g. distances) -------------------
source("scripts/ECC/classes_ecc.R")
source("scripts/ECC/dist_functions.R")
source("scripts/ECC/ecc_functions.R")
source("report_specific/epi-helper-no-source.R")
cluster_mapping <- readRDS("inputs/processed/allTP2.Rds")

dir.create("report_specific/heatmaps/clustering/", showWarnings = FALSE, recursive = TRUE)

fnames <- list.files("intermediate_data/TPN/dists/", full.names = TRUE)
assert("Distances were collected, script '3_dist_matrices.R' was run correctly", length(fnames) >= 1)
distfiles <- lapply(fnames, function(f) readRDS(f))

assert("Max and min of distances collected", file.exists("intermediate_data/TPN/extreme_dists.Rds"))
extremes <- readRDS("intermediate_data/TPN/extreme_dists.Rds")

arg <- parse_args(OptionParser(option_list=option_list))

params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

# The number of clusters to generate heatmaps for
# Will collect the largest ones, looking at the full data set (last time point, 
# when all strains are present)
# This is set in the "Generate heatmaps" line in the form_inputs.txt file
num_cl <- as.integer(params$numcl[2])

if (num_cl > 0) {
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
  
  # Filtering so we only look at clusters of size < max_cluster_size
  # This is because the heatmap function hangs when the data sets are too large
  # May need to find a different function for this
  # Note that we are using: 
  #   heatmap.2() called in EpiHeatmap_pdf() in report_specific/epi-helper-no-source.R
  top_clusters <- size_details %>% arrange(-TP2_cluster_size) %>% 
    filter(TP2_cluster_size < max_cluster_size) %>% 
    slice(1:num_cl) %>% pull(first_tp2_flag)
  
  # Matching these clusters to their original cluster names, for saving purposes
  clusters <- data.table(
    chr = top_clusters, 
    new_h = top_clusters %>% strsplit(., split = "_") %>% sapply(., '[[', 2) %>% 
      gsub("h", "", .) %>% as.double(), 
    new_cl = top_clusters %>% strsplit(., split = "_") %>% sapply(., '[[', 3) %>% 
      gsub("c", "", .) %>% as.integer()
  ) %>% inner_join(cluster_mapping$lookup_table, ., by = c("new_h", "new_cl")) %>% 
    select(chr, new_cl, old_h, old_cl) %>% 
    set_colnames(c("chr", "int", "original_h", "original_cl")) %>% 
    inner_join(., size_details, by = c("chr" = "first_tp2_flag")) %>% 
    arrange(-TP2_cluster_size)
  
  # # Reading in the last timepoint data, as well as the strain metadata and the 
  # # data representative assignments - used for processing distance files from intermediate_data/
  # tpn <- tp2 <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)
  # typing_data <- list(tp2$height_list)
  # 
  # m <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()
  # 
  # assignments <- typing_data[[1]] %>% as.data.frame() %>% set_colnames("tp_cl") %>% 
  #   rownames_to_column("Strain") %>% as.data.table() %>% 
  #   filter(tp_cl %in% clusters$int)
  # 
  # cl_drs <- lapply(1:nrow(clusters), function(i) {
  #   x1 <- assignments[tp_cl %in% clusters$int[i]] %>% pull(Strain)
  #   m$dr_matches %>% filter(Strain %in% x1) %>% pull(dr) %>% unique()
  # }) %>% set_names(clusters$int)
  # 
  for (i in 1:nrow(clusters)) {
    epitable <- readRDS(paste0("report_specific/epitables/", 
                               clusters$chr[i], "-", clusters$original_cl[i], ".Rds"))
    
    # The following two steps use:
    #   - distEpiMatrix() and 
    #   - EpiHeatmap_pdf(), 
    #   which can be found in report_specific/epi-helper-no-source.R
    
    # Processing epitables into temporal distance matrices, one for each cluster
    epi.matrix <- epitable %>% select(Strain.1, Strain.2, Temp.Dist)
    temp_mat <- distEpiMatrix(epi.matrix, "Temp.Dist")  
    
    # Processing epitables into geographical distance matrices, one for each cluster
    epi.matrix <- epitable %>% select(Strain.1, Strain.2, Geog.Dist)
    geo_mat <- distEpiMatrix(epi.matrix, "Geog.Dist")
    
    # Saving the temporal and geographical distance matrices as heatmaps
    outputDetails(paste0("Saving temporal and geographical heatmaps for cluster ", 
                         clusters$chr[i], " (", clusters$original_cl[i], ")", "..."), newcat = TRUE)
    
    if (clusters$TP2_cluster_size[i] > 1) {
      cl_id <- paste0(clusters$chr[i], "-", clusters$original_cl[i], "-after_last_TP")
      
      png(paste0("report_specific/heatmaps/", cl_id, "-temp.png"))
      tplot <- EpiHeatmap_pdf(temp_mat)
      tplot
      dev.off()
      
      png(paste0("report_specific/heatmaps/", cl_id, "-geo.png"))
      gplot <- EpiHeatmap_pdf(geo_mat)
      gplot
      dev.off()
      
      # Used this method (https://stackoverflow.com/questions/18354501/how-to-get-member-of-clusters-from-rs-hclust-heatmap-2)
      # to extract the clustering from the heatmaps
      temp_clustering <- cutree(as.hclust(tplot$rowDendrogram), 1:dim(temp_mat)[1])
      geo_clustering <- cutree(as.hclust(gplot$rowDendrogram), 1:dim(geo_mat)[1])
      list("temp" = temp_clustering, "geo" = geo_clustering) %>% 
        saveRDS(., paste0("report_specific/heatmaps/clustering/", 
                          clusters$chr[i], "-", clusters$original_cl[i], ".Rds"))
    }else {
      outputDetails(paste0("TP2 cluster ", clusters$chr[i], " has only one strain, no heatmap generated"))
    }  
  }

  outputDetails("\nHeatmaps saved to report_specific/heatmaps/", TRUE)
  outputDetails("Epitables saved to report_specific/epitables/", TRUE)
  outputDetails("Clustering done in the heatmaps saved to report_specific/heatmaps/clustering/", TRUE)
}

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), " Completed heatmap generation ", 
           paste0(rep("-", 30), collapse = ""), "||\n"))
