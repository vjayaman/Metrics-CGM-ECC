#! /usr/bin/env Rscript

# Will not generate heatmaps for cluster of size > max_cluster_size
# max_cluster_size <- 2000

# # General setup - calling libraries, etc. --------------------------------------------------------------
msg <- file("logs/heatmaps.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "fossil", "tidyr", "plyr", "dplyr", "readr", "heatmap3",
          "testit", "tibble", "reshape2", "RColorBrewer", "gplots", "data.table", "R6")
y <- suppressWarnings(
  suppressPackageStartupMessages(
    lapply(libs, require, character.only = TRUE)))

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), " Starting heatmap generation ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))

# This section collects user-specified inputs ----------------------------------------------------------
option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"), 
  make_option(c("-r", "--tp2RData"), metavar = "file", default = "inputs/processed/allTP2.Rds", 
              help = "TP2 processed data"))

# Reading in required functions, checking for pre-processed results (e.g. distances) -------------------
source("scripts/Misc/plotting_functions.R")
basedir <- "report_specific/heatmaps/"
dir.create(paste0(basedir, "clustering/"), showWarnings = FALSE, recursive = TRUE)

arg <- parse_args(OptionParser(option_list=option_list))
params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov", 
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl", 
              "clustby", "lowcol", "midcol", "highcol"))

# The number of clusters to generate heatmaps for:
# Will collect the largest ones, looking at the full data set (last time point,
# when all strains are present). Set in "Generate heatmaps" line in form_inputs.txt file
num_cl <- as.integer(params$numcl[2])

# the column to use for clustering information
cluster_by <- params$clustby[2]
# color scheme - should be user inputs
heatcolor <- colorRampPalette(c(params$lowcol[2],params$midcol[2],params$highcol[2]))(512)
col_dir <- c(params$lowcol[2], params$highcol[2]) #c("Low similarity", "High similarity")

if (num_cl > 0) {
  clustering_fname <- "intermediate_data/heatmap_cluster_labs.Rds"
  assert("7a_topX_distances.R was run beforehand", file.exists(clustering_fname))
  clusters <- readRDS(clustering_fname) %>% arrange(TP2_cluster_size)
  
  metadata <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()
  
  for (i in 1:nrow(clusters)) {
    epitable <- readRDS(paste0("report_specific/epitables/", 
                               clusters$chr[i], "-", clusters$original_cl[i], ".Rds")) %>% 
      mutate(Temp_Geo.Dist = (Temp.Dist + Geog.Dist)/2)
    
    # euclidean distances, later using single linkage clustering (these should be user inputs)
    # Processing epitables into temporal distance matrices, one for each cluster
    temp_mat <- prepEpiMatrix(epitable, "Temp.Dist")
    ordered_by_date <- metadata$strain_data %>% select(Strain, Date) %>% arrange(Date) %>%
      filter(Strain %in% rownames(temp_mat)) %>% pull(Strain)
    temp_mat <- temp_mat[ordered_by_date,ordered_by_date]
    
    # Processing epitables into geographical distance matrices, one for each cluster
    geo_mat <- prepEpiMatrix(epitable, "Geog.Dist")
    geo_mat <- geo_mat[ordered_by_date,ordered_by_date]
    
    # Processing epitables into geographical distance matrices, one for each cluster
    temp_geo_mat <- prepEpiMatrix(epitable, "Temp_Geo.Dist")
    temp_geo_mat <- temp_geo_mat[ordered_by_date,ordered_by_date]
    
    # Saving the temporal and geographical distance matrices as heatmaps
    outputDetails(paste0("Saving (temporal, geographical and average of both) heatmaps for cluster ", 
                         clusters$chr[i], " (", clusters$original_cl[i], ")", "..."), newcat = TRUE)
    
    if (clusters$TP2_cluster_size[i] > 1) {
      cl_id <- paste0(clusters$chr[i], "-", clusters$original_cl[i], "-after_last_TP")
      
      clustering_mat <- switch(cluster_by, "temp" = temp_mat, "geo" = geo_mat, "temp_geo" = temp_geo_mat)
      hc_choice <- hclust(as.dist(clustering_mat), method="single")
      
      # TEMPORAL ---------------------------------------------------------------------------------------
      outputDetails(paste0("   Generating heatmap for temporal data ..."), newcat = TRUE)
      heatmapPlot(paste0(basedir, cl_id, "-temp.png"), temp_mat, heatcolor, hc_choice, col_dir)
      
      outputDetails(paste0("   Generating density plot and frequency histogram ..."), newcat = TRUE)
      densityPlot(paste0(basedir, cl_id, "-temp-density.png"), "Temporal", temp_mat)
      frequencyPlot(paste0(basedir, cl_id, "-temp-frequencies.png"), "Temporal", temp_mat)

      
      # GEOGRAPHICAL -----------------------------------------------------------------------------------
      outputDetails(paste0("   Generating heatmap for geographical data ..."), newcat = TRUE)
      heatmapPlot(paste0(basedir, cl_id, "-geo.png"), geo_mat, heatcolor, hc_choice, col_dir)
      
      outputDetails(paste0("   Generating density plot and frequency histogram ..."), newcat = TRUE)
      densityPlot(paste0(basedir, cl_id, "-geo-density.png"), "Geographical", geo_mat)
      frequencyPlot(paste0(basedir, cl_id, "-geo-frequencies.png"), "Geographical", geo_mat)
      
      
      # TEMP AND GEO -----------------------------------------------------------------------------------
      outputDetails(paste0("   Generating heatmap for (temp+geo)/2 data ..."), newcat = TRUE)
      heatmapPlot(paste0(basedir, cl_id, "-tempgeo.png"), temp_geo_mat, heatcolor, hc_choice, col_dir)
      
      outputDetails(paste0("   Generating density plot and frequency histogram ..."), newcat = TRUE)
      densityPlot(paste0(basedir, cl_id, "-tempgeo-density.png"), "Temp and geo", temp_geo_mat)
      frequencyPlot(paste0(basedir, cl_id, "-tempgeo-frequencies.png"), "Temp and geo", temp_geo_mat)
      
      # Save clustering separately ---------------------------------------------------------------------
      # Used this method (https://stackoverflow.com/questions/18354501/how-to-get-member-of-clusters-from-rs-hclust-heatmap-2)
      # to extract the clustering from the heatmaps
      
      user_clustering <- cutree(hc_choice, 1:dim(clustering_mat)[1])
      saveRDS(user_clustering, paste0(basedir, "clustering/", clusters$chr[i], "-", clusters$original_cl[i], ".Rds"))
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
