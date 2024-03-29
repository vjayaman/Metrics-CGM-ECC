#! /usr/bin/env Rscript

# General setup - calling libraries, etc. --------------------------------------------------------------
# msg <- file("logs/pairwise_strain_distances.txt", open="wt")
# sink(msg, type="message")

libs <- c("optparse", "magrittr", "fossil", "tidyr", "plyr", "dplyr", "readr", 
          "testit", "tibble", "reshape2", "RColorBrewer", "gplots", "data.table", "R6")
y <- suppressWarnings(
  suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE)))

cat(paste0("\n||", paste0(rep("-", 21), collapse = ""), 
           " Starting pairwise (strain) distance collection ", 
           paste0(rep("-", 22), collapse = ""), "||\n"))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

# Adjust input parameters here: ------------------------------------------------------------------------
# This section collects user-specified inputs
source("scripts/arguments.R")

# Reading in required functions, checking for pre-processed results (e.g. distances) -------------------
source("scripts/ECC/classes_ecc.R")
source("scripts/ECC/dist_functions.R")
source("scripts/ECC/ecc_functions.R") # transformData2() is here
source("report_specific/epi-helper-no-source.R")

# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------
# COLLECTING DISTANCES 
# save_to <- file.path("intermediate_data", params$int_type[2], "TPN", "dists", "/")
# dir.create(save_to,  recursive = TRUE, showWarnings = FALSE)

hx <- strsplit(as.character(params$th[2]), split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

fdata <- readRDS(arg$tpn)$new_cols %>% column_to_rownames("Strain")
tp2 <- Timepoint$new(arg$tpn, "tp2", fdata)$Process(hx)$listHeights(hx)

m <- read_tsv(arg$metadata) %>% processedStrains()
# basedir <- file.path("intermediate_data", params$int_type[2], "TPN")


# Extremes, for scaling ----------------------------------------------------------------------
ext_geo_dists <- m$assignments %>% select(Longitude, Latitude) %>% unique() %>% 
  rownames_to_column("id") %>% distMatrix(., "geo", c("Longitude", "Latitude"))
ext_temp_dists <- m$assignments %>% select(Date) %>% unique() %>% 
  rownames_to_column("id") %>% distMatrix(., "temp", "Date")

extremes <- list(maxt = max(ext_temp_dists), mint = min(ext_temp_dists), 
                 maxg = max(ext_geo_dists), ming = min(ext_geo_dists))
# saveRDS(extremes, file.path(basedir, "extreme_dists.Rds"))
# rm(ext_geo_dists); rm(ext_temp_dists); rm(extremes)

# Pairwise distances within clusters at TPN --------------------------------------------------
outputMessages("Collecting non-redundant pairwise distances (at last timepoint, when dataset is full)")
metadata <- m$strain_data %>% as.data.table()

clustersets <- file.path("intermediate_data", params$int_type[2], "clustersets.Rds") %>% readRDS(.)
interval_list <- names(clustersets)
rm(clustersets)

k <- last(interval_list)

if (params$int_type[2] == "multiset") {
  interval <- "Multiset"
}else if (params$int_type[2] == "monthly") {
  interval <- "YearMonth"
}else if (params$int_type[2] == "weekly") {
  interval <- "YearWeek"
}

typing_data <- lapply(1:length(interval_list), function(i) {
  n1 <- as.character(interval_list[i])
  tpkstrains <- metadata[get(interval) <= n1]$Strain
  dfz <- tp2$filedata %>% rownames_to_column("isolate") %>%
    select(isolate, all_of(hx$h)) %>%
    filter(isolate %in% tpkstrains) %>% column_to_rownames("isolate")
  dfz[,hx$h[1],drop=FALSE] %>% set_colnames(hx$th[1])
}) %>% set_names(as.character(interval_list))

td <- typing_data[[length(typing_data)]] %>% rownames_to_column("Strain") %>% as.data.table()
rm(typing_data)

parts <- m$dr_matches %>% filter(Strain %in% td$Strain) %>% 
  left_join(td, ., by = "Strain") %>% sectionClusters(.)
saveRDS(parts, file.path("intermediate_data", params$int_type[2], "TPN", "parts.Rds"))

outputMessages(paste0("  Collecting and saving distances for cluster groups at TP", k, ":\n"))

tpkstrains <- metadata[get(interval) <= k]$Strain
collectDistances(parts$drs, parts$results, m$dr_matches, m$assignments, tpkstrains, save_to)

# --------------------------------------------------------------------------------------------
# fnames <- file.path("intermediate_data", params$int_type[2], "TPN/dists/") %>%
#   list.files(., full.names = TRUE)
# assert("Distances were collected, script '3_dist_matrices.R' was run correctly", length(fnames) >= 1)
# distfiles <- lapply(fnames, function(f) readRDS(f))

cluster_mapping <- readRDS(arg$tpn)

ext_file <- file.path("intermediate_data", params$int_type[2], "TPN/extreme_dists.Rds")
assert("Max and min of distances collected", file.exists(ext_file))
extremes <- readRDS(ext_file)

dir.create(file.path("results", params$int_type[2], "epitables"), 
           showWarnings = FALSE, recursive = TRUE)

# Reading in the threshold column (to get appropriate clusters)
hx <- strsplit(as.character(params$th[2]), split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

# Reading the CGM results (weekly/monthly/multiset)
cgms <- list.files(file.path("results", params$int_type[2]), full.names = TRUE) %>% 
  grep("CGM", ., value = TRUE) %>% readRDS()

num_cl <- as.numeric(params$numcl[2])

if (num_cl == 0) {
  cat(paste0("No data generated for heatmaps; user input gave ", 
             "'Generate heatmaps for top 0 largest clusters'"))
}else {
  # NON-REDUNDANT METHOD
  last_ivl <- unique(cgms$interval) %>% last()
  
  # Collecting the top x clusters by size (after all strains are in the data, last time point)
  # or by growth rate. x is given by the user in the form_input.txt file (the "Generate heatmaps" line)
  
  # "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)" = actual_growth_rate
  # "Novel growth = (TP2 size) / (TP2 size - number of novels)" = new_growth
  
  cols_to_consider <- cgms[interval == last_ivl] %>% 
    select(first_tp2_flag, tp2_cl_size, actual_growth_rate, new_growth) %>% unique() %>% 
    filter(tp2_cl_size >= (as.integer(params$mincl[2]) + 1)) 
  # we add 1 now because we used cluster size + 1 for the cgm columns
  
  top_x_by_size <- cols_to_consider %>% arrange(-tp2_cl_size) %>% slice(1:num_cl) %>% 
    mutate(by_size = 1) %>% select(first_tp2_flag, tp2_cl_size, by_size)
  
  top_x_by_gr <- cols_to_consider %>% arrange(-actual_growth_rate) %>% slice(1:num_cl) %>% 
    mutate(by_gr = 1) %>% select(first_tp2_flag, tp2_cl_size, by_gr)
  
  top_x_by_nr <- cols_to_consider %>% arrange(-new_growth) %>% slice(1:num_cl) %>% 
    mutate(by_nr = 1) %>% select(first_tp2_flag, tp2_cl_size, by_nr)
  
  manual_check <- cgms[interval == last_ivl] %>% 
    arrange(-tp2_cl_size) %>% select(first_tp2_flag, tp2_cl_size) %>% unique() %>% 
    filter(tp2_cl_size < 25 & tp2_cl_size >= 10) %>% slice(1) %>% 
    add_column(by_size = NA, by_gr = NA, by_nr = NA)
  
  top_x <- top_x_by_size %>% 
    merge.data.table(., top_x_by_gr, all.x = TRUE, all.y = TRUE) %>% 
    merge.data.table(., top_x_by_nr, all.x = TRUE, all.y = TRUE) %>% 
    bind_rows(manual_check) %>% select(-tp2_cl_size)
  
  x_clusters <- top_x %>% pull(first_tp2_flag)
  
  raw_sizes <- cluster_mapping$pango_clusters %>% select(Strain, old_h, old_cl, new_h, new_cl) %>% 
    group_by(old_h, old_cl, new_h, new_cl) %>% count() %>% ungroup() %>% as.data.table()
  colnames(raw_sizes)[colnames(raw_sizes) == "n"] <- "TP2_cluster_size"
  
  # Matching these clusters to their old cluster names, for saving purposes
  clusters <- data.table(
    chr = x_clusters,
    new_h = x_clusters %>% strsplit(., split = "_") %>% sapply(., '[[', 2) %>% 
      gsub("h", "", .) %>% as.double(), 
    new_cl = x_clusters %>% strsplit(., split = "_") %>% sapply(., '[[', 3) %>% 
      gsub("c", "", .) %>% as.integer()
  ) %>% inner_join(cluster_mapping$lookup_table, ., by = c("new_h", "new_cl")) %>% 
    select(chr, new_h, new_cl, old_h, old_cl) %>%
    set_colnames(c("chr", "new_h", "new_cl", "old_h", "old_cl")) %>% 
    full_join(., top_x, by = c("chr" = "first_tp2_flag")) %>% 
    merge.data.table(., raw_sizes, by = c("old_h", "old_cl", "new_h", "new_cl")) %>% 
    arrange(-TP2_cluster_size)
  
  file.path("results", params$int_type[2], "cluster_labs.Rds") %>% saveRDS(clusters, .)
  
  outputDetails(paste0("Collecting pairwise distances (from non-redundant set), ", 
                       "for each of the \n     ", nrow(clusters)-1, " clusters ", 
                       "from user input (top ", num_cl, " by size and top ", num_cl, 
                       " by rate (potential overlap)) \n          ", "Top ", num_cl, " by size: ", 
                       paste0(clusters[!is.na(by_size)]$chr, collapse = ", "), 
                       "\n          ", "Top ", num_cl, " by rate: ", 
                       paste0(clusters[!is.na(by_gr)]$chr, collapse = ", "), 
                       "\n     + 1 manual check\n"), 
                newcat = TRUE)
  
  # Reading in the last timepoint data, as well as the strain metadata and the 
  # data representative assignments - used for processing distance files from intermediate_data/
  
  fdata <- cluster_mapping$new_cols %>% column_to_rownames("Strain")
  tpn <- tp2 <- Timepoint$new(arg$tpn, "tp2", fdata)$Process(hx)$listHeights(hx)
  
  
  typing_data <- list(tp2$height_list)
  
  m <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()
  
  assignments <- typing_data[[1]] %>% as.data.frame() %>% set_colnames("tp_cl") %>% 
    rownames_to_column("Strain") %>% as.data.table() %>% 
    filter(tp_cl %in% clusters$new_cl)
  
  cl_drs <- lapply(1:nrow(clusters), function(i) {
    x1 <- assignments[tp_cl %in% clusters$new_cl[i]] %>% pull(Strain)
    m$dr_matches %>% filter(Strain %in% x1) %>% pull(dr) %>% unique()
  }) %>% set_names(clusters$new_cl)
  
  # Collecting epitable data for each of the clusters, with column names: 
  #       ||---------------------------------------------||
  #       || Strain.1 | Strain.2 | Temp.Dist | Geog.Dist ||
  #       ||---------------------------------------------||
  
  for (j in 1:nrow(clusters)) {
    clx <- clusters$new_cl[j]
    outputDetails(paste0("(", j, "/", nrow(clusters), ") Cluster ", 
                         gsub("TP2", "TPN", clusters$new_cl[j]), " (new) = ", 
                         clusters$old_cl[j], " (old), (of size ", 
                         clusters$TP2_cluster_size[j], "), "))
    
    if (length(cl_drs[[as.character(clx)]]) == 1) {
      outputDetails("This cluster is a singleton") 
    }
    
    rawdists <- lapply(1:length(distfiles), function(i) {
      drs <- cl_drs[[as.character(clx)]]
      
      tdm <- distfiles[[i]]$temp
      tdm_subset <- tdm[rownames(tdm) %in% drs, colnames(tdm) %in% drs]
      
      gdm <- distfiles[[i]]$geo  
      gdm_subset <- gdm[rownames(gdm) %in% drs, colnames(gdm) %in% drs]
      
      if (length(tdm_subset) == 0 & length(gdm_subset) == 0) {
        data.table("dr1" = NA, "dr2" = NA, "Temp.Dist" = NA, "Geog.Dist" = NA) 
      }else {
        if (as.logical(params$trheatmaps[2])) {
          ctdm <- tdm_subset %>% 
            # This is done using a function from scripts/ECC/ecc_functions.R
            transformData2(., "temp", extremes$mint, extremes$maxt) %>%
            formatData(., c("dr1","dr2","Temp.Dist"))
          
          cgdm <- gdm_subset %>% 
            # This is done using a function from scripts/ECC/ecc_functions.R
            transformData2(., "geo", extremes$ming, extremes$maxg) %>%
            formatData(., c("dr1","dr2","Geog.Dist"))
        }else {
          # DEFAULT: NOT TRANSFORMED - ADD TO THE FORM INPUTS FILE
          ctdm <- tdm_subset %>% formatData(., c("dr1","dr2","Temp.Dist"))
          cgdm <- gdm_subset %>% formatData(., c("dr1","dr2","Geog.Dist"))
        }
        
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
    file.path("results", params$int_type[2], "epitables", 
              paste0(clusters$chr[j], "-", clusters$old_cl[j], ".Rds")) %>% 
      saveRDS(pw_strains, file = .)
  }
}
# beep(3)

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
timeTaken(pt = "pairwise distances collection", stopwatch) %>% outputDetails(., newcat = TRUE)

cat(paste0("\n||", paste0(rep("-", 22), collapse = ""), " End of pairwise (strain) distance collection ", 
           paste0(rep("-", 23), collapse = ""), "||\n"))
