#! /usr/bin/env Rscript
x <- c("plyr", "dplyr", "tidyr", "readr", "stringr", "magrittr", "tibble", "purrr", "readxl", 
       "reshape2", "optparse")
lapply(x, require, character.only = TRUE)

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "Europe TP2 metadata.xlsx", help = "Metadata file"),
  
  make_option(c("-a", "--tp1"), metavar = "file", 
              default = "european-t1_clusters.csv", help = "TP1 cluster assignments"), 
  
  make_option(c("-b", "--tp2"), metavar = "file", 
              default = "european-t2_clusters.csv", help = "TP2 cluster assignments"), 
  
  make_option(c("-h", "--heights"), metavar = "character", default = "0,5", 
              help = "String of comma-delimited numbers for heights you want ECC input files for"))

arg <- parse_args(OptionParser(option_list=option_list))

checkAndSave <- function(x = TRUE, msg, filename, filedata) {
  if (x) {
    write.table(filedata, paste0(filename, ".txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  }else {
    stop(paste0("Problem with processing: ", msg))
  }
}

# Reading and processing the data in the full metadata excel file ----------------------------------------------
strain_data <- read_excel(arg$metadata, sheet = 4, col_names = TRUE, .name_repair = "minimal") %>% 
  select(c("Strain", "Source", "Country", "Province", "City", "Latitude", "Longitude", "Day", 
           "Month", "Year", "TP1", "TP2")) %>% 
  filter(Province != "England") %>% filter(Country != "United Kingdom") %>% 
  filter(TP2 == 1) %>% arrange(Strain)

# Reading and processing the cluster data for TP1 and TP2, making sure they match the metadata file ------------
# TP1 DATA -----------------------------------------------------------------------------------------------------
time1 <- read.csv(file = arg$tp1, stringsAsFactors = FALSE, 
                       numerals = "no.loss", check.names = FALSE, sep = ",") %>% as_tibble() %>% 
  set_colnames(c("Strain", 0:ncol(.))) %>% 
  mutate(Strain = gsub("hCoV-19/", "", Strain) %>% sub("\\|.*", "", .)) %>% 
  filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

checkAndSave(
  identical(strain_data$Strain[strain_data$TP1 == 1], time1$Strain), 
  "Not all TP1 strains accounted for", "tp1_clusters", time1)

# TP1 DATA -----------------------------------------------------------------------------------------------------
time2 <- read.csv(file = arg$tp2, stringsAsFactors = FALSE,
                       numerals = "no.loss", check.names = FALSE, sep = ",") %>% as_tibble() %>% 
  set_colnames(c("Strain", 0:ncol(.))) %>% 
  mutate(Strain = gsub("hCoV-19/", "", Strain) %>% sub("\\|.*", "", .)) %>% 
  filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

checkAndSave(
  identical(strain_data$Strain, time2$Strain), 
  "Not all TP2 strains accounted for", "tp2_clusters", time2)

# Extracting clusters for T0 from both TP1 and TP2, to run EpiQuant on -----------------------------------------
for (h in unlist(strsplit(arg$heights, split = ","))) {
  time1 %>% select(Strain, all_of(h)) %>% set_colnames(c("Strain", paste0("T", h))) %>% 
    checkAndSave(filename = paste0("TP1_T", h), filedata = .)
  
  time2 %>% select(Strain, all_of(h)) %>% set_colnames(c("Strain", paste0("T", h))) %>% 
    checkAndSave(filename = paste0("TP2_T", h), filedata = .)
}





