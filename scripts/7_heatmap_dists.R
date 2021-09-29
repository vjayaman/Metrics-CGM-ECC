#! /usr/bin/env Rscript

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
# THIS IS THE SECTION THAT USERS CAN/SHOULD CHANGE
option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  # make_option(c("-t", "--tau"), metavar = "numeric", help = "temporal coefficient", default = 0.0), 
  # make_option(c("-g", "--gamma"), metavar = "numeric", help = "geographic coefficient", default = 1.0), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"), 
  make_option(c("-s", "--maxclsize"), metavar = "numeric", 
              help = "Max cluster size to make heatmaps for", default = 1000))

# Reading in required functions, checking for pre-processed results (e.g. distances) -------------------
source("scripts/ECC/classes_ecc.R")
source("scripts/ECC/dist_functions.R")
source("scripts/ECC/ecc_functions.R")
source("report_specific/epi-helper-no-source.R")

dir.create("report_specific/heatmaps/", showWarnings = FALSE, recursive = TRUE)

fnames <- list.files("intermediate_data/TPN/dists/", full.names = TRUE)
assert("Distances were collected, script '3_dist_matrices.R' was run correctly", length(fnames) >= 1)
distfiles <- lapply(fnames, function(f) readRDS(f))

assert("Max and min of distances collected", file.exists("intermediate_data/TPN/extreme_dists.Rds"))
extremes <- readRDS("intermediate_data/TPN/extreme_dists.Rds")

arg <- parse_args(OptionParser(option_list=option_list))

params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

num_cl <- as.integer(params$numcl[2])

if (num_cl > 0) {
  # Reading in the threshold column (to get appropriate clusters)
  hx <- strsplit(as.character(params$th[2]), split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))
  
  # Reading the CGM results (weekly/monthly/multiset)
  results_files <- list.files("results/", full.names = TRUE) %>% grep(params$int_type[2], ., value = TRUE)
  cgms <- grep("CGM", results_files, value = TRUE) %>% readRDS()
  
  # NON-REDUNDANT METHOD
  last_ivl <- unique(cgms$interval) %>% last()
  
  size_details <- cgms[interval == last_ivl] %>% 
    arrange(-tp2_cl_size) %>% select(first_tp2_flag, tp2_cl_size) %>% unique()
  colnames(size_details)[which(colnames(size_details) == "tp2_cl_size")] <- "TP2_cluster_size"
  
  top_clusters <- size_details %>% arrange(-TP2_cluster_size) %>% 
    filter(TP2_cluster_size < arg$maxclsize) %>%
    slice(1:num_cl) %>% pull(first_tp2_flag)
  
  clusters <- strsplit(top_clusters, "_") %>% sapply(., `[[`, 3) %>% gsub("c", "", .) %>% 
    as.integer() %>% tibble(chr = top_clusters, int = .)
  
  tpn <- tp2 <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)
  typing_data <- list(tp2$height_list)
  
  m <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()
  
  assignments <- typing_data[[1]] %>% as.data.frame() %>% set_colnames("tp_cl") %>% 
    rownames_to_column("Strain") %>% as.data.table() %>% 
    filter(tp_cl %in% clusters$int)
  
  cl_drs <- lapply(1:nrow(clusters), function(i) {
    x1 <- assignments[tp_cl %in% clusters$int[i]] %>% pull(Strain)
    m$dr_matches %>% filter(Strain %in% x1) %>% pull(dr) %>% unique()
  }) %>% set_names(clusters$int)
  
  epi.tables <- lapply(1:nrow(clusters), function(j) {
    clx <- clusters$int[j]
    
    rawdists <- lapply(1:length(distfiles), function(i) {
      drs <- cl_drs[[as.character(clx)]]
      
      tdm <- distfiles[[i]]$temp
      ctdm <- tdm[rownames(tdm) %in% drs, colnames(tdm) %in% drs] %>% 
        # FEEL FREE to comment out this part if you don't want transformed data
        # This is done using a function from scripts/ECC/ecc_functions.R
        transformData2(., "temp", extremes$mint, extremes$maxt) %>%
        formatData(., c("dr1","dr2","Temp.Dist"))  
      
      gdm <- distfiles[[i]]$geo  
      cgdm <- gdm[rownames(gdm) %in% drs, colnames(gdm) %in% drs] %>% 
        # FEEL FREE to comment out this part if you don't want transformed data
        # This is done using a function from scripts/ECC/ecc_functions.R
        transformData2(., "geo", extremes$ming, extremes$maxg) %>%
        formatData(., c("dr1","dr2","Geog.Dist"))
      
      merge.data.table(ctdm, cgdm)
    }) %>% bind_rows()
    
    x1 <- assignments[tp_cl %in% clx] %>% pull(Strain)
    x2 <- m$dr_matches %>% filter(Strain %in% x1) %>% as.data.table()
    x3 <- inner_join(x2, rawdists, by = c("dr" = "dr1")) %>% select(-dr)
    colnames(x3)[which(colnames(x3)=="Strain")] <- "Strain.1"
      
    dists <- inner_join(x3, x2, by = c("dr2" = "dr"))
    colnames(dists)[which(colnames(dists)=="Strain")] <- "Strain.2"
    
    dists %>% select(Strain.1, Strain.2, Temp.Dist, Geog.Dist)
  }) %>% set_names(clusters$chr)
  
  saveRDS(epi.tables, "report_specific/heatmaps/epitables_for_heatmaps.Rds")
  
  # distEpiMatrix() and EpiHeatmap_pdf() can be found in report_specific/epi-helper-no-source.R
  
  temp_matrices <- lapply(clusters$chr, function(cl_j) {
    epi.matrix <- epi.tables[[cl_j]] %>% select(Strain.1, Strain.2, Temp.Dist)
    distEpiMatrix(epi.matrix, "Temp.Dist")
  }) %>% set_names(clusters$chr)
  
  geo_matrices <- lapply(clusters$chr, function(cl_j) {
    epi.matrix <- epi.tables[[cl_j]] %>% select(Strain.1, Strain.2, Geog.Dist)
    distEpiMatrix(epi.matrix, "Geog.Dist")
  }) %>% set_names(clusters$chr)
  
  for (i in 1:nrow(clusters)) {
    outputDetails(paste0("Saving temporal and geographical distances for cluster ", 
                         clusters$chr[i], "..."), newcat = TRUE)
    cl_strains <- size_details %>% filter(first_tp2_flag %in% clusters$chr[i])
    if (cl_strains$TP2_cluster_size > 1) {
      cl_id <- cgms %>% filter(first_tp2_flag == clusters$chr[i]) %>% slice(1) %>% 
        pull(first_tp2_flag)
      
      cl_epi_temp <- temp_matrices[[clusters$chr[i]]]
      cl_epi_geo <- geo_matrices[[clusters$chr[i]]]
      
      png(paste0("report_specific/heatmaps/", cl_id, "-temp.png"))
      EpiHeatmap_pdf(cl_epi_temp)
      dev.off()
      
      png(paste0("report_specific/heatmaps/", cl_id, "-geo.png"))
      EpiHeatmap_pdf(cl_epi_geo)
      dev.off()
    }else {
      outputDetails(paste0("TP2 cluster ", top_clusters[i], " has only one strain, no heatmap generated"))
    }
  }
  
  num_results <- list.files("report_specific/heatmaps/") %>% length()
  if (num_results == num_cl+1) {
    outputDetails("All heatmaps (and epitables file) saved to report_specific/heatmaps/", newcat = TRUE)
  }else {
    outputDetails(paste0("Not all expected files saved to report_specific/heatmaps/, ", 
                         "something went wrong (see logs)"), newcat = TRUE)
  }
  
}else {
  outputDetails(paste0("From the form.html results, we have 'Generate heatmaps for top 0 ", 
                       "largest clusters'.\nTo change this, change the corresponding value in ", 
                       "inputs/form_inputs.txt"), newcat = TRUE)
}

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), " Completed heatmap generation ", 
           paste0(rep("-", 30), collapse = ""), "||\n"))
