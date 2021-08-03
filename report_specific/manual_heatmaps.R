#! /usr/bin/env Rscript

libs <- c("optparse", "magrittr", "fossil", "tidyr", "plyr", "dplyr", "readr", 
          "testit", "tibble", "reshape2", "RColorBrewer", "gplots")
y <- suppressWarnings(
  suppressPackageStartupMessages(
    lapply(libs, require, character.only = TRUE)
  )
)
# assert("All packages loaded correctly", all(unlist(y)))

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-c", "--CGMs"), metavar = "file", 
              default = "results/CGM_strain_results.tsv", help = "CGM result file"),
  make_option(c("-t", "--tau"), metavar = "numeric", help = "temporal coefficient", 
              default = 0.0), 
  make_option(c("-g", "--gamma"), metavar = "numeric", help = "geographic coefficient", 
              default = 1.0), 
  make_option(c("-l", "--clusters"), metavar = "numeric", 
              help = "Number of clusters to get heatmaps for", default = 5)
)

arg <- parse_args(OptionParser(option_list=option_list))

strain_data <- suppressMessages(read_tsv(arg$metadata)) %>% 
  mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
         Location = paste(Country, Province, City, sep = "_"))

# MANUAL METHOD: 
source("report_specific/epi-helper-no-source.R")
dir.create("report_specific/manual_heatmaps/", showWarnings = FALSE)
epi.table <- EpiTableNoSource(strain_data, temp_coeff = arg$tau, geog_coeff = arg$gamma)
saveRDS(epi.table, "report_specific/manual_heatmaps/epitable.Rds")

epi.matrix <- EpiMatrix(epi.table)

### Section 3: Heatmap for selected clusters
cgms <- suppressMessages(read_tsv("results/CGM_strain_results.tsv"))
size_details <- cgms %>%
  select(Strain, tp2_cl, tp2_cl_size) %>%
  rename(TP2_cluster_size = tp2_cl_size)# %>% as.data.table()

top_five <- size_details %>% select(-Strain) %>% unique() %>%
  arrange(-TP2_cluster_size) %>% filter(TP2_cluster_size < 10000) %>%
  slice(1:arg$clusters) %>% pull(tp2_cl)

saveRDS(epi.matrix, "report_specific/manual_heatmaps/epimatrix.Rds")

for (i in 1:length(top_five)) {
  cl_strains <- size_details %>% filter(tp2_cl %in% top_five[i])
  if (length(cl_strains$Strain) > 1) {
    cl_epi <- epi.matrix[rownames(epi.matrix) %in% cl_strains$Strain,
                         colnames(epi.matrix) %in% cl_strains$Strain]
    cl_id <- cgms %>% filter(tp2_cl == top_five[i]) %>% slice(1) %>% pull(first_tp2_flag)
    
    png(paste0("report_specific/manual_heatmaps/", cl_id, ".png"))
    EpiHeatmap_pdf(cl_epi)
    dev.off()
  }else {
    print(paste0("TP2 cluster ", top_five[i], " has only one strain, no heatmap generated"))
  }
}

