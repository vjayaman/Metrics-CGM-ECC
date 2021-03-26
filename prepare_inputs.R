#! /usr/bin/env Rscript

msg <- file("logs/logfile_inputs.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "readr", "dplyr")
y <- lapply(libs, require, character.only = TRUE)

option_list <- list(
  make_option(c("-i", "--inputdir"), metavar = "dir", default = "inputs/", help = "Raw inputs directory"), 
  make_option(c("-m", "--metadata"), metavar = "file", default = "strain_info.txt", help = "Metadata file"),
  make_option(c("-a", "--tp1"), metavar = "file", default = "tp1_base_clusters.txt", help = "TP1 cluster assignments"), 
  make_option(c("-b", "--tp2"), metavar = "file", default = "tp2_base_clusters.txt", help = "TP2 cluster assignments"))

arg <- parse_args(OptionParser(option_list=option_list))

intClusters <- function(df) {
  list_lookup <- lapply(2:ncol(df), function(j) {
    original_values <- df %>% pull(j) %>% unique()
    tibble(original_values, 1:length(original_values)) %>% 
      set_names(c(names(df)[j], paste0("new_", j-2)))
  })
  
  lookup_tbl <- suppressMessages(plyr::join_all(append(list(df), list_lookup), type = "left")) %>% as_tibble()
  cnames <- grep("Strain|new", colnames(lookup_tbl), value = TRUE)
  
  list(lookup_tbl, lookup_tbl[,cnames] %>% set_colnames(gsub("new_", "", colnames(.)))) %>% 
    set_names(c("lookup_table", "new_cols"))%>% return()
}

writeData <- function(df, filepath) {
  write.table(df, filepath, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

# checkAndSave <- function(x = TRUE, msg, filepath, filedata, ext) {
#   if (x) {
#     write.table(filedata, filepath, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#   }else {
#     stop(paste0("Problem with processing: ", msg))
#   }
# }

# Reading and processing the data in the full metadata excel file ----------------------------------------------
# read_excel("input_data/Europe TP2 metadata.xlsx", sheet = 4, col_names = TRUE, .name_repair = "minimal")
strain_data <- file.path(arg$inputdir, arg$metadata) %>% 
  read.table(sep = "\t", header = TRUE, fill = TRUE, quote = "") %>% as_tibble() %>% 
  select(c("Strain", "Source", "Country", "Province", "City", "Latitude", 
           "Longitude", "Day", "Month", "Year", "TP1", "TP2")) %>%
  filter(TP2 == 1) %>%
  arrange(Strain) %>% 
  mutate(TP1 = ifelse(is.na(TP1), 0, TP1))

# writeData(strain_data, file.path(arg$inputdir, "processed", "strain_info.txt"))
# strain_data <- file.path("inputs/strain_data.tsv") %>% read_tsv() %>%
#   select(c("Strain", "Source", "Country", "Province", "City", "Latitude", "Longitude", "Day",
#            "Month", "Year", "TP1", "TP2")) %>% filter(TP2 == 1) %>% arrange(Strain) %>%
#   filter(Province != "England") %>% filter(Country != "United Kingdom")
# writeData(strain_data, filepath = file.path("inputs/processed/strain_info.txt"))
# Reading and processing the cluster data for TP1 and TP2, making sure they match the metadata file ------------
# TP1 DATA -----------------------------------------------------------------------------------------------------
time1 <- file.path(arg$inputdir, arg$tp1) %>% 
  read.table(sep = "\t", header = TRUE) %>% as_tibble() %>% 
  filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

processed_tp1 <- intClusters(time1)

writeData(processed_tp1$lookup_table, file.path(arg$inputdir, "processed", "tp1_lookup.txt"))
writeData(processed_tp1$new_cols, file.path(arg$inputdir, "processed", "tp1_clusters.txt"))

# TP1 DATA -----------------------------------------------------------------------------------------------------
time2 <- file.path(arg$inputdir, arg$tp2) %>% 
  read.table(sep = "\t", header = TRUE) %>% as_tibble() %>% 
  filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

processed_tp2 <- intClusters(time2)

writeData(processed_tp2$lookup_table, file.path(arg$inputdir, "processed", "tp2_lookup.txt"))
writeData(processed_tp2$new_cols, file.path(arg$inputdir, "processed", "tp2_clusters.txt"))

# time2 <- file.path(arg$inputdir, arg$tp2) %>%
#   read.csv(., stringsAsFactors = FALSE, check.names = FALSE, sep = ",") %>% as_tibble() %>%
#   rename(Strain = isolate) %>% 
  # set_colnames(c("Strain", 0:ncol(.))) %>%
  # mutate(Strain = gsub("hCoV-19/", "", Strain) %>% sub("\\|.*", "", .))
# writeData(time2, file.path(arg$inputdir, "tp2_base_clusters.txt"))
  # filter(Strain %in% strain_data$Strain) %>% arrange(Strain)
# checkAndSave(identical(strain_data$Strain, time2$Strain), 
#   "Not all TP2 strains accounted for", file.path(arg$inputdir, "processed", "tp2_clusters.txt"), time2)

# ECC-SPECIFIC INPUT FILES -------------------------------------------------------------------------------------
# placeholder source file --------------------------------------------------------------------------------------

tibble(Source.1 = "Placeholder1", Source.2 = "Placeholder2", value = 0) %>% 
  write.table(., file.path("inputs", "processed", "source_data.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)




