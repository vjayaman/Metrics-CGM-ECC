#! /usr/bin/env Rscript

msg <- file("logs/logfile_inputs.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "readr", "dplyr")
y <- lapply(libs, require, character.only = TRUE)

option_list <- list(
  make_option(c("-i", "--inputdir"), metavar = "dir", default = "inputs/", help = "Raw inputs directory"), 
  make_option(c("-m", "--metadata"), metavar = "file", default = "strain_data.tsv", help = "Metadata file"),
  make_option(c("-a", "--tp1"), metavar = "file", default = "european-t1_clusters.csv", help = "TP1 cluster assignments"), 
  make_option(c("-b", "--tp2"), metavar = "file", default = "european-t2_clusters.csv", help = "TP2 cluster assignments"), 
  make_option(c("-h", "--heights"), metavar = "character", default = "0,5", 
              help = "String of comma-delimited numbers for heights you want ECC input files for"))

arg <- parse_args(OptionParser(option_list=option_list))

checkAndSave <- function(x = TRUE, msg, filepath, filedata, ext) {
  if (x) {
    write.table(filedata, filepath, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  }else {
    stop(paste0("Problem with processing: ", msg))
  }
}

# Reading and processing the data in the full metadata excel file ----------------------------------------------
# read_excel("input_data/Europe TP2 metadata.xlsx", sheet = 4, col_names = TRUE, .name_repair = "minimal")
strain_data <- file.path(arg$inputdir, arg$metadata) %>% read_tsv() %>% 
  select(c("Strain", "Source", "Country", "Province", "City", "Latitude", "Longitude", "Day", 
           "Month", "Year", "TP1", "TP2")) %>% 
  filter(Province != "England") %>% filter(Country != "United Kingdom") %>% 
  filter(TP2 == 1) %>% arrange(Strain)

# Reading and processing the cluster data for TP1 and TP2, making sure they match the metadata file ------------
# TP1 DATA -----------------------------------------------------------------------------------------------------
time1 <- file.path(arg$inputdir, arg$tp1) %>% 
  read.csv(., stringsAsFactors = FALSE, check.names = FALSE, sep = ",") %>% as_tibble() %>% 
  set_colnames(c("Strain", 0:ncol(.))) %>% 
  mutate(Strain = gsub("hCoV-19/", "", Strain) %>% sub("\\|.*", "", .)) %>% 
  filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

checkAndSave(
  identical(strain_data$Strain[strain_data$TP1 == 1], time1$Strain), 
  "Not all TP1 strains accounted for", file.path(arg$inputdir, "processed", "tp1_clusters.txt"), time1)

# TP1 DATA -----------------------------------------------------------------------------------------------------
time2 <- file.path(arg$inputdir, arg$tp2) %>% 
  read.csv(., stringsAsFactors = FALSE, check.names = FALSE, sep = ",") %>% as_tibble() %>% 
  set_colnames(c("Strain", 0:ncol(.))) %>% 
  mutate(Strain = gsub("hCoV-19/", "", Strain) %>% sub("\\|.*", "", .)) %>% 
  filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

checkAndSave(
  identical(strain_data$Strain, time2$Strain), 
  "Not all TP2 strains accounted for", file.path(arg$inputdir, "processed", "tp2_clusters.txt"), time2)

# ECC-SPECIFIC INPUT FILES -------------------------------------------------------------------------------------
# placeholder source file --------------------------------------------------------------------------------------

tibble(Source.1 = "Placeholder1", Source.2 = "Placeholder2", value = 0) %>% 
  write.table(., file.path("inputs", "processed", "source_data.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)




