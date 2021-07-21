#! /usr/bin/env Rscript

msg <- file("logs/logfile_alldists.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", 
           "scripts/ECC/dist_functions.R")
invisible(sapply(files, source))

option_list <- list(
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-t", "--typedist"), metavar = "character", default = "temp", 
              help = paste0("Provide either 'temp' or 'geo', and the full set (strain) ", 
                            "pairwise distances will be saved. Note that they are saved ", 
                            "in parts for very large datasets, to save memory and processing power.")))

cat(paste0("\n||", paste0(rep("-", 22), collapse = ""), 
           " Saving full set of (strain) pairwise distances ", 
           paste0(rep("-", 21), collapse = ""), "||\nStarted process at: ", Sys.time()))

params <- parse_args(OptionParser(option_list=option_list))

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

m <- read_tsv(params$strains) %>% processedStrains()

# "Save distances as a full matrix, e.g. for heatmaps. ", 
# "Default is FALSE, uses reduced redundancy methods for distance metrics"))

outputMessages(paste0("\nNote that this may take some time and memory, as we are merging ", 
                      "the smaller, \nnon-redundant set with all possible pairs ", 
                      "of strains and then saving the results."))
k <- 2

intra_cl <- paste0("intermediate_data/TP", k, "/dists/") %>% 
  list.files(., full.names = TRUE, pattern = "group") %>% 
  lapply(., function(fname) {readRDS(fname)})

parts <- sectionClusters(k, typing_data, m)
rm(typing_data)
rm(tp1)
rm(tp2)

dr_matches <- parts$drs[,c("Strain","dr")] %>% as_tibble()

if (params$typedist == "temp") {
  type <- c("temp", "Temp.Dist", "Date", "temporal_distances")
}else if (params$typedist == "geo") {
  type <- c("geo", "Geog.Dist", "Longitude,Latitude", "geographical_distances")
}else {
  outputMessages("Incorrect 'type' argument provided. Either 'temp' or 'geo' work at the moment.")
}

dir.create(paste0("results/", type[4]), showWarnings = FALSE)

outputMessages(paste0("\nFormatting and merging ", gsub("_", " ", type[4]), " ...\n\n"))
  
intra_temps <- lapply(1:length(intra_cl), function(i) {
  intra_cl[[i]][[type[1]]] %>% as.data.frame() %>% rownames_to_column("dr") %>% 
    as.data.table() %>% melt.data.table(id.vars = "dr") %>% 
    set_colnames(c("dr1", "dr2", type[2]))
}) %>% bind_rows() %>% unique()
  
inter_temps <- paste0("intermediate_data/TP", k, "/dists/inter_dists_", type[1], ".Rds") %>% 
  readRDS() %>% as.data.frame() %>% rownames_to_column("dr") %>% as.data.table() %>% 
  melt.data.table(id.vars = "dr") %>% 
  set_colnames(c("dr1", "dr2", type[2]))
  
compiled_temps <- bind_rows(inter_temps, intra_temps) %>% unique() %>% 
  mutate(across(dr2, as.character)) %>% as.data.table()
rm(inter_temps)
rm(intra_temps)
gc()
  
p1 <- left_join(compiled_temps, m$dr_matches, by = c("dr1" = "dr")) %>% rename(Strain1 = Strain)
rm(compiled_temps)
gc()
gc()
  
if (nrow(p1) < 5000000) {
  p2 <- right_join(dr_matches, p1, by = c("dr" = "dr2")) %>% rename(Strain2 = Strain)
  p2 <- right_join(dr_matches, p1, by = c("dr" = "dr2")) %>% rename(Strain2 = Strain)
  p3 <- p2 %>% as.data.table() %>% select(Strain1, Strain2, all_of(type[2])) %>% 
    set_colnames(c("Strain1", "Strain2", type[2]))
  saveRDS(p3, paste0("results/", type[4], "/alldists.Rds"))
}else {
  indices <- split(1:nrow(p1), ceiling(seq_along(1:nrow(p1)) / 5000000))
  gc(full = TRUE, verbose = FALSE)
  
  for (j in 1:length(indices)) {
    p2 <- right_join(dr_matches, p1[indices[[j]],], by = c("dr" = "dr2")) %>% 
      rename(Strain2 = Strain)
    gc(full = TRUE, verbose = FALSE)
    
    p3 <- p2 %>% as.data.table() %>% select(Strain1, Strain2, all_of(type[2])) %>% 
      set_colnames(c("Strain1", "Strain2", type[2]))
    gc(full = TRUE, verbose = FALSE)
    
    saveRDS(p3, paste0("results/", type[4], "/part", j, "of", length(indices), ".Rds"))
    
    rm(p2); rm(p3); gc(full = TRUE, verbose = FALSE)
    }
}

outputMessages(paste0("See ", paste0("results/", type[4], "/"), 
                      " for the saved (strain) pairwise distances."))

cat(paste0("\n||", paste0(rep("-", 27), collapse = ""), 
           " End of processing full distances set ", 
           paste0(rep("-", 26), collapse = ""), "||\n"))
