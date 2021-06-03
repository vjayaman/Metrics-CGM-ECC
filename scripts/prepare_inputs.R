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

outputDetails(paste0("\n||", paste0(rep("-", 26), collapse = ""), " Prepping inputs for metric generation ", 
                     paste0(rep("-", 26), collapse = ""), "||\nStarted process at: ", Sys.time()), TRUE)
outputDetails(paste0("\nWill save formatted inputs to 'processed' directory in ", arg$inputdir, " directory."), TRUE)

# Results of "Form for analysis inputs" ------------------------------------------------------------------------
filtering_params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>% 
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP1","nsTP2","temp_win","cnames"))

strain_data <- readData(arg$metadata, TRUE)
time1 <- readData(arg$tp1)
time2 <- readData(arg$tp2)

# LINEAGE INFO -------------------------------------------------------------------------------------------------
initial_sizes <- tibble(type="initial", a=nrow(strain_data), b=nrow(time1), d=nrow(time2))

assert("Has lineage info", as.logical(filtering_params$has_lin[2]))
x <- updateStrains("lin_info", strain_data, time1, time2, initial_sizes)
strain_data <- x$sd; time1 <- x$t1; time2 <- x$t2; initial_sizes <- x$sizes

# COLUMN NAMES -------------------------------------------------------------------------------------------------
cnames <- filtering_params$cnames[2] %>% strsplit(split = ",") %>% unlist()
if (!("none" %in% cnames)) {
  fullcnames <- c("Strain", "Latitude", "Longitude", "Day", "Month", "Year", cnames)
}else {
  fullcnames <- c("Strain", "Latitude", "Longitude", "Day", "Month", "Year")
}
strain_data <- strain_data %>% select(all_of(fullcnames))

# add column to show which strains are found in TP1
strain_data %<>% mutate(TP1 = ifelse(Strain %in% time1$Strain, 1, 0))

# add column to show which strains are found in TP2
strain_data %<>% mutate(TP2 = ifelse(Strain %in% time2$Strain, 1, 0))

# REGION OF INTEREST -------------------------------------------------------------------------------------------
if (filtering_params$reg[2] != "All" & "Region" %in% colnames(strain_data)) {
  strain_data <- strain_data %>% filter(Region %in% filtering_params$reg[2])
}

# COUNTRY OF INTEREST ------------------------------------------------------------------------------------------
if (filtering_params$cou[2] != "All" & "Country" %in% colnames(strain_data)) {
  strain_data <- strain_data %>% filter(Country %in% filtering_params$cou[2])
}

# HAS DEFINED DATE INFO ----------------------------------------------------------------------------------------
if (as.logical(filtering_params$has_date[2])) {
  strain_data <- strain_data %>% 
    filter(!is.na(Day)) %>% filter(!is.na(Month)) %>% filter(!is.na(Year)) %>% 
    filter(Day != "") %>% filter(Month != "") %>% filter(Year != "")
}

if (nchar(filtering_params$temp_win[2]) > nchar("[,]")) {
  tempwindow <- filtering_params$temp_win[2] %>% gsub("\\[|\\]", "", .) %>% 
    strsplit(., ",") %>% unlist()
  
  strain_data %<>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
    filter(Date >= tempwindow[1] & Date <= tempwindow[2])

  # Strains with metadata and defined lineage info at TP1
  y <- updateStrains("temp_win", strain_data, time1, time2, initial_sizes)
  strain_data <- y$sd; time1 <- y$t1; time2 <- y$t2; initial_sizes <- y$sizes
}

# HAS PROVINCE-LEVEL DATA --------------------------------------------------------------------------------------

if (as.logical(filtering_params$has_prov[2])) {
  strain_data <- strain_data %>% filter(!is.na(Province)) %>% filter(Province != "")
  
  if (filtering_params$prov[2] != "All") {
    strain_data <- strain_data %>% filter(Province %in% filtering_params$prov[2])
  }
  
  # Strains with metadata and defined lineage info at TP2
  z <- updateStrains("aft_prov", strain_data, time1, time2, initial_sizes)
  strain_data <- z$sd; time1 <- z$t1; time2 <- z$t2; initial_sizes <- z$sizes
}

# NON-SINGLETON CLUSTERS ---------------------------------------------------------------------------------------
outputDetails("Making table for matching TP1 clusters to integers (for metrics process) ...", TRUE)
processed_tp1 <- intClusters(time1)

outputDetails("Making table for matching TP2 clusters to integers (for metrics process) ...", TRUE)
processed_tp2 <- intClusters(time2)

th <- filtering_params$th[2]

if (as.logical(filtering_params$nsTP1[2])) {
  remove_strains <- strainsInSingletons(processed_tp1$new_cols, th)
  processed_tp1$new_cols %<>% filter(!(Strain %in% remove_strains))
  processed_tp2$new_cols %<>% filter(!(Strain %in% remove_strains))
  strain_data %<>% filter(!(Strain %in% remove_strains))
  
  initial_sizes %<>% add_row(tibble(type="tp1_ns", a=nrow(strain_data), b=nrow(time1), d=nrow(time2)))
}

if (as.logical(filtering_params$nsTP2[2])) {
  remove_strains <- strainsInSingletons(processed_tp2$new_cols, th)
  processed_tp1$new_cols %<>% filter(!(Strain %in% remove_strains))
  processed_tp2$new_cols %<>% filter(!(Strain %in% remove_strains))
  strain_data %<>% filter(!(Strain %in% remove_strains))
  
  initial_sizes %<>% add_row(tibble(type="tp2_ns", a=nrow(strain_data), b=nrow(time1), d=nrow(time2)))
}

writeData(strain_data, file.path(arg$inputdir, "processed", "strain_info.txt"))

writeData(processed_tp1$new_cols, file.path(arg$inputdir, "processed", "tp1_clusters.txt"))
saveRDS(processed_tp1, file.path(arg$inputdir, "processed", "allTP1.Rds"))

writeData(processed_tp2$new_cols, file.path(arg$inputdir, "processed", "tp2_clusters.txt"))
saveRDS(processed_tp2, file.path(arg$inputdir, "processed", "allTP2.Rds"))

outputDetails(paste0("\nFinished process at: ", Sys.time(), "\n||", paste0(rep("-", 14), collapse = ""), " Saved formatted inputs to 'processed' in the ",
                     arg$inputdir, " directory", paste0(rep("-", 15), collapse = ""), "||"), TRUE)
