#! /usr/bin/env Rscript

# General setup - calling libraries, etc. --------------------------------------------------------------
msg <- file("logs/manual_heatmap_check.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "fossil", "tidyr", "plyr", "dplyr", "readr", "heatmap3",
          "testit", "tibble", "reshape2", "RColorBrewer", "gplots", "data.table", "R6")
y <- suppressWarnings(
  suppressPackageStartupMessages(
    lapply(libs, require, character.only = TRUE)))

cat(paste0("\n||", paste0(rep("-", 24), collapse = ""), " Starting (manual check) heatmap generation ", 
           paste0(rep("-", 23), collapse = ""), "||\n"))

# Adjust input parameters here: ------------------------------------------------------------------------
# This section collects user-specified inputs
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

fnames <- list.files("intermediate_data/TPN/dists/", full.names = TRUE)
assert("Distances were collected, script '3_dist_matrices.R' was run correctly", length(fnames) >= 1)
distfiles <- lapply(fnames, function(f) readRDS(f))

assert("Max and min of distances collected", file.exists("intermediate_data/TPN/extreme_dists.Rds"))
extremes <- readRDS("intermediate_data/TPN/extreme_dists.Rds")

arg <- parse_args(OptionParser(option_list=option_list))

cluster_mapping <- readRDS(arg$tp2RData)

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
  
# Filtering so we only look at small clusters: this is for manually checking the heatmaps
manual_check <- size_details %>% 
  filter(TP2_cluster_size < 25 & TP2_cluster_size >= 10) %>% slice(1) %>% 
  pull(first_tp2_flag)
  
# Matching these clusters to their original cluster names, for saving purposes
clusters <- data.table(
  chr = manual_check, 
  new_h = manual_check %>% strsplit(., split = "_") %>% sapply(., '[[', 2) %>% 
    gsub("h", "", .) %>% as.double(), 
  new_cl = manual_check %>% strsplit(., split = "_") %>% sapply(., '[[', 3) %>% 
    gsub("c", "", .) %>% as.integer()) %>% 
  
  inner_join(cluster_mapping$lookup_table, ., by = c("new_h", "new_cl")) %>% 
  
  select(chr, new_cl, old_h, old_cl) %>% 
  set_colnames(c("chr", "int", "original_h", "original_cl")) %>% 
  
  inner_join(., size_details, by = c("chr" = "first_tp2_flag")) %>% 
  arrange(-TP2_cluster_size)


  
metadata <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()
  
# color scheme - should be user inputs
heatcolor <- colorRampPalette(c("white","yellowgreen","darkgreen"))(512)
  
i <- 1
epitable <- readRDS(paste0("report_specific/epitables/", 
                           clusters$chr[i], "-", clusters$original_cl[i], ".Rds"))
    
# euclidean distances, later using single linkage clustering (these should be user inputs)
    
# Processing epitables into temporal distance matrices, one for each cluster
epi.matrix <- epitable %>% select(Strain.1, Strain.2, Temp.Dist)
temp_mat <- distEpiMatrix(epi.matrix, "Temp.Dist")  
    
ordered_by_date <- metadata$strain_data %>% select(Strain, Date) %>% arrange(Date) %>%
  filter(Strain %in% rownames(temp_mat)) %>% pull(Strain)
temp_mat <- temp_mat[ordered_by_date,ordered_by_date]
    
# Processing epitables into geographical distance matrices, one for each cluster
epi.matrix <- epitable %>% select(Strain.1, Strain.2, Geog.Dist)
geo_mat <- distEpiMatrix(epi.matrix, "Geog.Dist")
geo_mat <- geo_mat[ordered_by_date,ordered_by_date]
    
# Saving the temporal and geographical distance matrices as heatmaps
outputDetails(paste0("Saving temporal and geographical heatmaps for cluster ", 
                     clusters$chr[i], " (", clusters$original_cl[i], ")", "..."), newcat = TRUE)

cl_id <- paste0(clusters$chr[i], "-", clusters$original_cl[i], "-after_last_TP")
heatcolor1<- colorRampPalette(c("white","yellowgreen","darkgreen"))(512)
thc <- hclust(as.dist(temp_mat), method="single")  

png(paste0("report_specific/heatmaps/", cl_id, "-temp-manual-check.png"))
heatmap3(temp_mat, col=rev(heatcolor1), 
                  Rowv=as.dendrogram(thc), Colv=as.dendrogram(thc), revC = T, scale = 'none', 
                  margins = c(10,10), xlab=NULL, ylab=NULL, method = "complete", useRaster = TRUE, 
                  main = paste0("Pairwise distances for temporal data"), 
                  legendfun = function() showLegend(legend=c("Low similarity", "High similarity"), 
                                                    col = c("white", "darkgreen")))
dev.off()
  
png(paste0("report_specific/heatmaps/", cl_id, "-geo-manual-check.png"))
heatmap3(geo_mat, col=rev(heatcolor1), 
                  Rowv=as.dendrogram(thc), Colv=as.dendrogram(thc), revC = T, scale = 'none', 
                  margins = c(10,10), xlab=NULL, ylab=NULL, method = "complete", useRaster = TRUE, 
                  main = paste0("Pairwise distances for geographical data"), 
                  legendfun = function() showLegend(legend=c("Low similarity", "High similarity"), 
                                                    col = c("white", "darkgreen")))
dev.off()

outputDetails("\nHeatmaps saved to report_specific/heatmaps/", TRUE)
cat(paste0("\n||", paste0(rep("-", 23), collapse = ""), " Completed (manual check) heatmap generation ", 
           paste0(rep("-", 23), collapse = ""), "||\n"))


  