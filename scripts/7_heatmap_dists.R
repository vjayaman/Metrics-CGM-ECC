#! /usr/bin/env Rscript

msg <- file("logs/heatmaps.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "fossil", "tidyr", "plyr", "dplyr", "readr", 
          "testit", "tibble", "reshape2", "RColorBrewer", "gplots", "data.table", "R6")
y <- suppressWarnings(
  suppressPackageStartupMessages(
    lapply(libs, require, character.only = TRUE)))

source("scripts/ECC/classes_ecc.R")
source("scripts/ECC/dist_functions.R")
source("scripts/ECC/ecc_functions.R")
source("report_specific/epi-helper-no-source.R")

dir.create("report_specific/heatmaps/", showWarnings = FALSE)

fnames <- list.files("intermediate_data/TPN/dists/", full.names = TRUE)
distfiles <- lapply(fnames, function(f) readRDS(f))
extremes <- readRDS("intermediate_data/TPN/extreme_dists.Rds")

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  make_option(c("-c", "--CGMs"), metavar = "file", 
              default = "results/CGM-monthly-intervals.Rds", help = "CGM result file"),
  make_option(c("-t", "--tau"), metavar = "numeric", help = "temporal coefficient", default = 0.0), 
  make_option(c("-g", "--gamma"), metavar = "numeric", help = "geographic coefficient", default = 1.0), 
  make_option(c("-l", "--clusters"), metavar = "numeric", 
              help = "Number of clusters to get heatmaps for", default = 5),
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"), 
  make_option(c("-s", "--maxclsize"), metavar = "numeric", 
              help = "Max cluster size to make heatmaps for", default = 1000))

arg <- parse_args(OptionParser(option_list=option_list))

params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

hx <- strsplit(as.character(params$th[2]), split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

cgms <- readRDS(arg$CGMs)

strain_data <- suppressMessages(read_tsv(arg$metadata)) %>% 
  mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
         Location = paste(Country, Province, City, sep = "_"))

# NON-REDUNDANT METHOD
last_ivl <- unique(cgms$interval) %>% last()

size_details <- cgms[interval == last_ivl] %>% 
  arrange(-tp2_cl_size) %>% select(first_tp2_flag, tp2_cl_size) %>%
  rename(TP2_cluster_size = tp2_cl_size) %>% unique()

top_clusters <- size_details %>% arrange(-TP2_cluster_size) %>% 
  filter(TP2_cluster_size < arg$maxclsize) %>%
  slice(1:arg$clusters) %>% pull(first_tp2_flag)

clusters <- substr(top_clusters, 11, 13) %>% as.integer()

tpn <- tp2 <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- list(tp2$height_list)

m <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()

assignments <- typing_data[[1]] %>% as.data.frame() %>% set_colnames("tp_cl") %>% 
  rownames_to_column("Strain") %>% as.data.table() %>% 
  filter(tp_cl %in% clusters)

cl_drs <- lapply(1:length(clusters), function(i) {
  x1 <- assignments[tp_cl %in% clusters[i]] %>% pull(Strain)
  m$dr_matches %>% filter(Strain %in% x1) %>% pull(dr) %>% unique()
}) %>% set_names(clusters)

epi.tables <- lapply(1:length(clusters), function(j) {
  clx <- clusters[j]
  
  rawdists <- lapply(1:length(distfiles), function(i) {
    drs <- cl_drs[[as.character(clx)]]
    
    tdm <- distfiles[[i]]$temp
    ctdm <- tdm[rownames(tdm) %in% drs, colnames(tdm) %in% drs] %>% 
      transformData2(., "temp", extremes$mint, extremes$maxt) %>% 
      formatData(., c("dr1","dr2","Temp.Dist"))
    
    gdm <- distfiles[[i]]$geo
    cgdm <- gdm[rownames(gdm) %in% drs, colnames(gdm) %in% drs] %>% 
      transformData2(., "geo", extremes$ming, extremes$maxg) %>%
      formatData(., c("dr1","dr2","Geog.Dist"))
    
    merge.data.table(ctdm, cgdm)
  }) %>% bind_rows()
  
  x1 <- assignments[tp_cl %in% clx] %>% pull(Strain)
  x2 <- m$dr_matches %>% filter(Strain %in% x1) %>% as.data.table()
  x3 <- inner_join(x2, rawdists, by = c("dr" = "dr1")) %>% select(-dr) %>% rename(Strain.1 = Strain)
  
  dists <- inner_join(x3, x2, by = c("dr2" = "dr")) %>% rename(Strain.2 = Strain) %>% 
    select(Strain.1, Strain.2, Temp.Dist, Geog.Dist) %>% 
    mutate(Total.Dist = sqrt( (((Temp.Dist^2)*arg$tau) + ((Geog.Dist^2)*arg$gamma)) ),
           Epi.Sym = 1 - Total.Dist)
  
  return(dists)
})

saveRDS(epi.tables, "report_specific/heatmaps/epitables_for_heatmaps.Rds")

epi.matrix <- lapply(1:length(clusters), function(j) EpiMatrix(epi.tables[[j]]))

for (i in 1:length(top_clusters)) {
  cl_strains <- size_details %>% filter(first_tp2_flag %in% top_clusters[i])
  if (cl_strains$TP2_cluster_size > 1) {
    cl_epi <- epi.matrix[[i]]
    cl_id <- cgms %>% filter(first_tp2_flag == top_clusters[i]) %>% slice(1) %>% 
      pull(first_tp2_flag)
    
    png(paste0("report_specific/heatmaps/", cl_id, ".png"))
    EpiHeatmap_pdf(cl_epi)
    dev.off()
  }else {
    print(paste0("TP2 cluster ", top_clusters[i], " has only one strain, no heatmap generated"))
  }
}
