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

dir.create("report_specific/heatmaps/clustering/", showWarnings = FALSE, recursive = TRUE)

arg <- parse_args(OptionParser(option_list=option_list))
params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

# The number of clusters to generate heatmaps for:
num_cl <- as.integer(params$numcl[2])
# Will collect the largest ones, looking at the full data set (last time point,
# when all strains are present). Set in "Generate heatmaps" line in form_inputs.txt file

if (num_cl > 0) {
  clustering_fname <- "intermediate_data/heatmap_cluster_labs.Rds"
  assert("7a_topX_distances.R was run beforehand", file.exists(clustering_fname))
  clusters <- readRDS(clustering_fname) %>% arrange(TP2_cluster_size)
  
  metadata <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()
  
  # color scheme - should be user inputs
  heatcolor <- colorRampPalette(c("white","yellowgreen","darkgreen"))(512)
  
  for (i in 1:nrow(clusters)) {
    epitable <- readRDS(paste0("report_specific/epitables/", 
                               clusters$chr[i], "-", clusters$original_cl[i], ".Rds"))
    epitable <- epitable %>% mutate(Temp_Geo.Dist = (Temp.Dist + Geog.Dist)/2)
    
    # euclidean distances, later using single linkage clustering (these should be user inputs)
    # Processing epitables into temporal distance matrices, one for each cluster
    temp_epi.matrix <- epitable %>% select(Strain.1, Strain.2, Temp.Dist)
    temp_mat <- distEpiMatrix(temp_epi.matrix, "Temp.Dist")  
    rm(temp_epi.matrix); gc()
    
    ordered_by_date <- metadata$strain_data %>% select(Strain, Date) %>% arrange(Date) %>%
      filter(Strain %in% rownames(temp_mat)) %>% pull(Strain)
    # if you want to print the matrices to check: 
    temp_mat <- temp_mat[ordered_by_date,ordered_by_date]
    
    
    # Processing epitables into geographical distance matrices, one for each cluster
    geo_epi.matrix <- epitable %>% select(Strain.1, Strain.2, Geog.Dist)
    geo_mat <- distEpiMatrix(geo_epi.matrix, "Geog.Dist")
    rm(geo_epi.matrix); gc()
    geo_mat <- geo_mat[ordered_by_date,ordered_by_date]
    
    
    # Processing epitables into geographical distance matrices, one for each cluster
    tg_epi.matrix <- epitable %>% select(Strain.1, Strain.2, Temp_Geo.Dist)
    temp_geo_mat <- distEpiMatrix(tg_epi.matrix, "Temp_Geo.Dist")
    rm(tg_epi.matrix); gc()
    temp_geo_mat <- temp_geo_mat[ordered_by_date,ordered_by_date]
    
    
    # Saving the temporal and geographical distance matrices as heatmaps
    outputDetails(paste0("Saving (temporal, geographical and average of both) heatmaps for cluster ", 
                         clusters$chr[i], " (", clusters$original_cl[i], ")", "..."), newcat = TRUE)
    
    if (clusters$TP2_cluster_size[i] > 1) {
      cl_id <- paste0(clusters$chr[i], "-", clusters$original_cl[i], "-after_last_TP")
      
      heatcolor1 <- colorRampPalette(c("white","yellowgreen","darkgreen"))(512)
      cluster_by <- "temp"
      clustering_mat <- switch(cluster_by, "temp" = temp_mat, "geo" = geo_mat, "temp_geo" = temp_geo_mat)
      thc <- hclust(as.dist(clustering_mat), method="single")
      # thc <- hclust(as.dist(temp_mat), method="single")
      
      
      # TEMPORAL ---------------------------------------------------------------------------------------
      outputDetails(paste0("   Generating heatmap for temporal data ..."), newcat = TRUE)
      
      png(paste0("report_specific/heatmaps/", cl_id, "-temp.png"))
      tplot <- heatmap3(temp_mat, col=rev(heatcolor1), labRow = NA, labCol = NA, 
                        Rowv=as.dendrogram(thc), Colv=as.dendrogram(thc), revC = T, scale = 'none', 
                        margins = c(10,10), xlab=NULL, ylab=NULL, method = "complete", useRaster = TRUE, 
                        # main = paste0("Pairwise distances for temporal data"), 
                        legendfun = function() showLegend(legend=c("Low similarity", "High similarity"), 
                                                          col = c("white", "darkgreen")))
      tplot
      dev.off()
      
      outputDetails(paste0("   Generating density plot and frequency histogram ..."), newcat = TRUE)
      png(paste0("report_specific/heatmaps/", cl_id, "-temp-density.png"))
      d <- density(temp_mat, adjust = 0.25)
      plot(d, main = "Temporal distances density", xlab = "Temporal pairwise distances")
      dev.off()
      
      png(paste0("report_specific/heatmaps/", cl_id, "-temp-frequencies.png"))
      half_size <- temp_mat
      half_size[upper.tri(half_size)] <- 0
      hist(as.vector(half_size), freq = TRUE, main = "Temporal distances counts", 
           xlab = "Temporal pairwise distances")
      dev.off()
      
      # GEOGRAPHICAL -----------------------------------------------------------------------------------
      outputDetails(paste0("   Generating heatmap for geographical data ..."), newcat = TRUE)
      
      png(paste0("report_specific/heatmaps/", cl_id, "-geo.png"))
      gplot <- heatmap3(geo_mat, col=rev(heatcolor1), labRow = NA, labCol = NA, 
                        Rowv=as.dendrogram(thc), Colv=as.dendrogram(thc), revC = T, scale = 'none', 
                        margins = c(10,10), xlab=NULL, ylab=NULL, method = "complete", useRaster = TRUE, 
                        # main = paste0("Pairwise distances for geographical data"), 
                        legendfun = function() showLegend(legend=c("Low similarity", "High similarity"), 
                                                          col = c("white", "darkgreen")))
      gplot
      dev.off()  
      
      outputDetails(paste0("   Generating density plot and frequency histogram ..."), newcat = TRUE)
      png(paste0("report_specific/heatmaps/", cl_id, "-geo-density.png"))
      d <- density(geo_mat, adjust = 0.25)
      plot(d, main = "Geographical distances density", xlab = "Geographical pairwise distances")
      dev.off()
      
      png(paste0("report_specific/heatmaps/", cl_id, "-geo-frequencies.png"))
      half_size <- geo_mat
      half_size[upper.tri(half_size)] <- 0
      hist(as.vector(half_size), freq = TRUE, main = "Geographical distances counts", 
           xlab = "Geographical pairwise distances")
      dev.off()
      
      
      # TEMP AND GEO -----------------------------------------------------------------------------------
      outputDetails(paste0("   Generating heatmap for (temp+geo)/2 data ..."), newcat = TRUE)
      
      png(paste0("report_specific/heatmaps/", cl_id, "-tempgeo.png"))
      tgplot <- heatmap3(temp_geo_mat, col=rev(heatcolor1), labRow = NA, labCol = NA, 
                        Rowv=as.dendrogram(thc), Colv=as.dendrogram(thc), revC = T, scale = 'none', 
                        margins = c(10,10), xlab=NULL, ylab=NULL, method = "complete", useRaster = TRUE, 
                        # main = paste0("Pairwise distances for temp and geo data"), 
                        legendfun = function() showLegend(legend=c("Low similarity", "High similarity"), 
                                                          col = c("white", "darkgreen")))
      tgplot
      dev.off()  
      
      outputDetails(paste0("   Generating density plot and frequency histogram ..."), newcat = TRUE)
      png(paste0("report_specific/heatmaps/", cl_id, "-tempgeo-density.png"))
      d <- density(temp_geo_mat, adjust = 0.25)
      plot(d, main = "Temp and geo distances density", xlab = "Temp and geo pairwise distances")
      dev.off()
      
      png(paste0("report_specific/heatmaps/", cl_id, "-tempgeo-frequencies.png"))
      half_size <- temp_geo_mat
      half_size[upper.tri(half_size)] <- 0
      hist(as.vector(half_size), freq = TRUE, main = "Temp and geo distances counts", 
           xlab = "Temp and geo pairwise distances")
      dev.off()
      
      
      # Save clustering separately ---------------------------------------------------------------------
      # Used this method (https://stackoverflow.com/questions/18354501/how-to-get-member-of-clusters-from-rs-hclust-heatmap-2)
      # to extract the clustering from the heatmaps
      
      temp_clustering <- cutree(thc, 1:dim(temp_mat)[1])
      # geo_clustering <- cutree(ghc, 1:dim(geo_mat)[1])
      
      list("temp" = temp_clustering) %>% #, "geo" = geo_clustering) %>%
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
