#! /usr/bin/env Rscript

msg <- file("logs/logfile_heatmaps.txt", open="wt")
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

fnames <- list.files("intermediate_data/TP2/dists/", full.names = TRUE)
distfiles <- lapply(fnames, function(f) readRDS(f))
extremes <- readRDS("intermediate_data/dist_extremes.Rds")

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-a", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  make_option(c("-c", "--CGMs"), metavar = "file", 
              default = "results/CGM_strain_results.tsv", help = "CGM result file"),
  make_option(c("-t", "--tau"), metavar = "numeric", help = "temporal coefficient", 
              default = 0.0), 
  make_option(c("-g", "--gamma"), metavar = "numeric", help = "geographic coefficient", 
              default = 1.0), 
  make_option(c("-l", "--clusters"), metavar = "numeric", 
              help = "Number of clusters to get heatmaps for", default = 5),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for")
  )

arg <- parse_args(OptionParser(option_list=option_list))

cgms <- suppressMessages(read_tsv(arg$CGMs))

strain_data <- suppressMessages(read_tsv(arg$metadata)) %>% 
  mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
         Location = paste(Country, Province, City, sep = "_"))

# NON-REDUNDANT METHOD
size_details <- cgms %>% arrange(-tp2_cl_size) %>% select(Strain, tp2_cl, tp2_cl_size) %>%
  rename(TP2_cluster_size = tp2_cl_size)

top_clusters <- size_details %>% select(-Strain) %>% unique() %>%
  arrange(-TP2_cluster_size) %>% filter(TP2_cluster_size < 10000) %>%
  slice(1:arg$clusters) %>% pull(tp2_cl)

clusters <- top_clusters %>% gsub("c", "", .) %>% as.integer()

hx <- arg$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(arg$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

m <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()

assignments <- typing_data %>% extract2(2) %>% set_colnames("tp_cl") %>% 
  rownames_to_column("Strain") %>% filter(tp_cl %in% clusters) %>% as.data.table()

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
  cl_strains <- size_details %>% filter(tp2_cl %in% top_clusters[i])
  if (length(cl_strains$Strain) > 1) {
    cl_epi <- epi.matrix[[i]]
    cl_id <- cgms %>% filter(tp2_cl == top_clusters[i]) %>% slice(1) %>% pull(first_tp2_flag)
    
    png(paste0("report_specific/heatmaps/", cl_id, ".png"))
    EpiHeatmap_pdf(cl_epi)
    dev.off()
  }else {
    print(paste0("TP2 cluster ", top_clusters[i], " has only one strain, no heatmap generated"))
  }
}
