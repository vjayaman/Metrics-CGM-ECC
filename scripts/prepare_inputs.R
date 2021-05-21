#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_inputs.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "readr", "dplyr", "testit")
y <- lapply(libs, require, character.only = TRUE)

option_list <- list(
  make_option(c("-i", "--inputdir"), metavar = "dir", default = "inputs", help = "Raw inputs directory"), 
  make_option(c("-m", "--metadata"), metavar = "file", default = "inputs/strain_info.txt", help = "Metadata file"),
  make_option(c("-a", "--tp1"), metavar = "file", default = "inputs/tp1_clusters_init.txt", help = "TP1 cluster assignments"), 
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/tp2_clusters_init.txt", help = "TP2 cluster assignments"))

arg <- parse_args(OptionParser(option_list=option_list))

dir.create("results", showWarnings = FALSE)
dir.create(file.path(arg$inputdir, "processed"))

# FUNCTIONS ----------------------------------------------------------------------------------------------------
checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

# Outputs the same message in two ways, one is directed to standard output and one to a log file
outputDetails <- function(msg, newcat = FALSE) {
  cat(msg)
  if (newcat) {cat("\n")}
  message(msg)
}

# Writes data to a given location, saves as tab-delimited text file
writeData <- function(df, filepath) {
  write.table(df, filepath, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

# Creates a lookup table of integer cluster (in case clusters with character names are provided)
# Returns lookup table and new integer columns, with integer thresholds as well, from 0 to # of thresholds - 1
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

outputDetails(paste0("\n||", paste0(rep("-", 20), collapse = ""), " Preparing inputs for CGM and ECC metric generation ", 
           paste0(rep("-", 20), collapse = ""), "||\nStarted process at: ", Sys.time()), TRUE)
outputDetails(paste0("\nWill save formatted inputs to new 'processed' directory in ", arg$inputdir, " directory."), TRUE)

# Reading and processing the data in the full metadata excel file ----------------------------------------------
outputDetails("Reading in strain data ...", TRUE)

filtering_params <- readLines("analysis_inputs.txt", warn = FALSE) %>% 
  strsplit(., split = ": ") %>% 
  set_names(c("geo", "lin", "date", "prov", "th", "ns", "temp_win"))

strain_data <- file.path(arg$metadata) %>% 
  read.table(sep = "\t", header = TRUE, fill = TRUE, quote = "", 
             fileEncoding = checkEncoding(.)) %>% as_tibble() %>% 
  select(c("Strain", "Source", "Country", "Province", "City", "Latitude", 
           "Longitude", "Day", "Month", "Year", "TP1", "TP2")) %>%
  filter(TP2 == 1) %>%
  arrange(Strain) %>% 
  mutate(TP1 = ifelse(is.na(TP1), 0, TP1))

# Reading and processing the cluster data for TP1 and TP2, making sure they match the metadata file ------------
# TP1 DATA -----------------------------------------------------------------------------------------------------
time1 <- file.path(arg$tp1) %>% 
  read.table(sep = "\t", header = TRUE, fileEncoding = checkEncoding(.)) %>% as_tibble() %>% 
  filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

outputDetails("Making table for matching TP1 clusters to integers (for metrics process) ...", TRUE)
processed_tp1 <- intClusters(time1)

# TP1 DATA -----------------------------------------------------------------------------------------------------
time2 <- file.path(arg$tp2) %>% 
  read.table(sep = "\t", header = TRUE, fileEncoding = checkEncoding(.)) %>% as_tibble() %>% 
  filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

outputDetails("Making table for matching TP2 clusters to integers (for metrics process) ...", TRUE)
processed_tp2 <- intClusters(time2)

# SAVE PROCESSED FILES -----------------------------------------------------------------------------------------

# geographical area of interest
# if (filtering_params$geo[2] != "All") {
#   strain_data %>% filter()
# }

# has defined lineage information
if (filtering_params$lin[2] == "true") {
  tp1_lineage_strains <- processed_tp1$new_cols$Strain
  tp2_lineage_strains <- processed_tp2$new_cols$Strain

  strain_data <- strain_data %>% filter(Strain %in% c(tp1_lineage_strains, tp2_lineage_strains))  
}

# has defined date information (day, month, and year)
if (filtering_params$date[2] == "true") {
  strain_data <- strain_data %>% 
    filter(!is.na(Day)) %>% filter(!is.na(Month)) %>% filter(!is.na(Year)) %>% 
    filter(Day != "") %>% filter(Month != "") %>% filter(Year != "")
}

# has provincial-level data
if (filtering_params$prov[2] == "true") {
  strain_data <- strain_data %>% 
    filter(!is.na(Province)) %>% 
    filter(Province != "")
}

# is in a non-singleton cluster (at TP1)
if (filtering_params$ns[2] == "true") {
  assert("Threshold provided", !is.null(filtering_params$th[2]))
  singletons <- processed_tp1$lookup_table %>% select(Strain, filtering_params$th[2]) %>% 
    set_colnames(c("Strains", "cluster")) %>% group_by(cluster) %>% 
    summarise(cluster_size = n()) %>% 
    filter(cluster_size == 1) %>% pull(cluster)  
  
  in_singletons <- processed_tp1$lookup_table %>% 
    filter(!!as.symbol(filtering_params$th[2]) %in% singletons) %>% 
    pull(Strain)
  
  processed_tp1$new_cols <- processed_tp1$new_cols %>% filter(!(Strain %in% in_singletons))
  processed_tp2$new_cols <- processed_tp2$new_cols %>% filter(!(Strain %in% in_singletons))
  strain_data <- strain_data %>% filter(!(Strain %in% in_singletons))
}

# within a specified temporal window


writeData(strain_data, file.path(arg$inputdir, "processed", "strain_info.txt"))
writeData(processed_tp1$new_cols, file.path(arg$inputdir, "processed", "tp1_clusters.txt"))
writeData(processed_tp2$new_cols, file.path(arg$inputdir, "processed", "tp2_clusters.txt"))

outputDetails(paste0("\nFinished process at: ", Sys.time(), "\n||", paste0(rep("-", 14), collapse = ""), " Saved formatted inputs to 'processed' in the ", 
                     arg$inputdir, " directory", paste0(rep("-", 15), collapse = ""), "||"), TRUE)


