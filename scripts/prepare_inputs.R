#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_inputs.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "readr", "dplyr", "testit")
y <- lapply(libs, require, character.only = TRUE)

option_list <- list(
  make_option(c("-i", "--inputdir"), metavar = "dir", 
              default = "inputs", help = "Raw inputs directory"), 
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/strain_info.txt", help = "Metadata file"),
  make_option(c("-a", "--tp1"), metavar = "file", 
              default = "inputs/tp1_clusters_init.txt", help = "TP1 cluster assignments"), 
  make_option(c("-b", "--tp2"), metavar = "file", 
              default = "inputs/tp2_clusters_init.txt", help = "TP2 cluster assignments"), 
  make_option(c("-x", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))

arg <- parse_args(OptionParser(option_list=option_list))

dir.create("results", showWarnings = FALSE)
dir.create(file.path(arg$inputdir, "processed"))

# FUNCTIONS ----------------------------------------------------------------------------------------------------
source("scripts/prepfunctions.R")

outputDetails(paste0("\n||", paste0(rep("-", 20), collapse = ""), " Prepping inputs for metric generation ", 
                     paste0(rep("-", 20), collapse = ""), "||\nStarted process at: ", Sys.time()), TRUE)
outputDetails(paste0("\nWill save formatted inputs to 'processed' directory in ", arg$inputdir, " directory."), TRUE)

# Results of "Form for analysis inputs" ------------------------------------------------------------------------
filtering_params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>% 
  set_names(c("reg","cou","has_date","has_prov","prov","th","ns","temp_win","cnames"))

strain_data <- file.path(arg$metadata) %>% 
  read.table(sep = "\t", header = TRUE, fill = TRUE, quote = "", 
             fileEncoding = checkEncoding(.)) %>% as_tibble()

time1 <- file.path(arg$tp1) %>% 
  read.table(sep = "\t", header = TRUE, fileEncoding = checkEncoding(.)) %>% as_tibble()

time2 <- file.path(arg$tp2) %>% 
  read.table(sep = "\t", header = TRUE, fileEncoding = checkEncoding(.)) %>% as_tibble()

# LINEAGE INFO -------------------------------------------------------------------------------------------------

# Strains with metadata and defined lineage info at TP1
tp1strains <- intersect(strain_data$Strain, time1$Strain)

time1 <- time1 %>% filter(Strain %in% tp1strains)

# Strains with metadata and defined lineage info at TP1
tp2strains <- intersect(strain_data$Strain, time2$Strain)

time2 <- time2 %>% filter(Strain %in% tp2strains)

# Strains that have defined lineage info
strain_data <- strain_data %>% filter(Strain %in% c(tp1strains, tp2strains))

# COLUMN NAMES -------------------------------------------------------------------------------------------------
cnames <- filtering_params$cnames[2] %>% strsplit(split = ",") %>% unlist()
if (cnames != "none") {
  strain_data <- strain_data %>% 
    select("Strain", "Latitude", "Longitude", "Day", "Month", "Year", all_of(cnames))
}

# REGION OF INTEREST -------------------------------------------------------------------------------------------
if (filtering_params$reg[2] != "All" & "Region" %in% colnames(strain_data)) {
  strain_data <- strain_data %>% filter(Region %in% filtering_params$reg[2])
}

# COUNTRY OF INTEREST ------------------------------------------------------------------------------------------
if (filtering_params$cou[2] != "All" & "Country" %in% colnames(strain_data)) {
  strain_data <- strain_data %>% filter(Country %in% filtering_params$cou[2])
}

# HAS DEFINED DATE INFO ----------------------------------------------------------------------------------------
# if (!as.logical(filtering_params$has_date[2])) {
#   
# }

# HAS PROVINCE-LEVEL DATA --------------------------------------------------------------------------------------
if (as.logical(filtering_params$has_prov[2])) {
  if (filtering_params$prov[2] != "All") {
    strain_data <- strain_data %>% filter(Province %in% filtering_params$prov[2])
  }
}

# NON-SINGLETON CLUSTERS ---------------------------------------------------------------------------------------
if (as.logical(filtering_params$ns[2])) {
  time1
}

# TEMPORAL WINDOW ----------------------------------------------------------------------------------------------


# reg "Region of interest" "All"
# cou "Country of interest" "All"
# has_date "Has defined date information (day, month, and year)" "true"
# has_prov "Has province-level data" "true"
# prov "Province of interest" "All"
th "Threshold of interest" "T0"
ns "Is in a non-singleton cluster (at TP1)" "false"
temp_win "Filtering by date" "none"
# cnames "Column names" "Source,Country,Province,City,TP2"



# REQUIRED
# Strain, Latitude, Longitude, Day (for now), Month, Year

# Optional
# Source, Region, Country, Province, City, TP1, TP2

# if ("TP2" %in% cnames) {
#   strain_data <- strain_data %>% filter(TP2 == 1)
# }else {
#   
# }
# 
# if ("TP1" %in% cnames) {
#   strain_data <- strain_data %>% mutate(TP1 = ifelse(is.na(TP1), 0, TP1))
# }else {
#   
# }

# Reading and processing the cluster data for TP1 and TP2, making sure they match the metadata file ------------
# TP1 DATA -----------------------------------------------------------------------------------------------------

# provided TP1 lineages

  
filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

# x1 <- time1 %>% select(filtering_params$th[2]) %>% pull()
# 
# if (!all(is.integer(x1))) {
outputDetails("Making table for matching TP1 clusters to integers (for metrics process) ...", TRUE)
processed_tp1 <- intClusters(time1)  

# TP1 DATA -----------------------------------------------------------------------------------------------------



  filter(Strain %in% strain_data$Strain) %>% arrange(Strain)

outputDetails("Making table for matching TP2 clusters to integers (for metrics process) ...", TRUE)
processed_tp2 <- intClusters(time2)

# SAVE PROCESSED FILES -----------------------------------------------------------------------------------------

# geographical area of interest
if (filtering_params$geo[2] != "All") {
  strain_data %>% filter()
}else {
  
}

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


