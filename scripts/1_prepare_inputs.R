#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_inputs.txt", open="wt")
sink(msg, type="message")

libs <- c("optparse", "magrittr", "readr", "dplyr", "testit", "data.table", "tibble")
y <- lapply(libs, require, character.only = TRUE); rm(y); rm(libs)

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/strain_info.txt", help = "Metadata file"),
  make_option(c("-b", "--tp2"), metavar = "file", 
              default = "inputs/tp2_clusters_init.txt", help = "TP2 cluster assignments"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))

arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)

# FUNCTIONS ----------------------------------------------------------------------------------------------------
source("scripts/Misc/formatprep.R")

outputDetails(paste0("\n||", paste0(rep("-", 26), collapse = ""), " Prepping inputs for metric generation ",
                     paste0(rep("-", 26), collapse = ""), "||\nStarted process at: ", Sys.time()), TRUE)
outputDetails(paste0("\nWill save formatted inputs to 'processed' directory in inputs directory."), TRUE)

# Results of "Form for analysis inputs" ------------------------------------------------------------------------
params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ")

test_params <- c("Region of interest", "Country of interest", "Has defined lineage information", 
                 "Has defined date information (day, month, and year)", "Has province-level data", 
                 "Province of interest", "Threshold of interest", 
                 "Is in a non-singleton cluster (at TP2)", "Filtering by date", "Column names", 
                 "Interval type", "Dividers", "Source-temporal-geographic coefficents", 
                 "Generate heatmaps for top __ largest clusters")

assert("Input parameters are correctly labelled", identical(sapply(params, '[[', 1), test_params))
rm(test_params)

params %<>% set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
                        "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

a1 <- readData(arg$metadata, FALSE)
a2 <- suppressWarnings(readData(arg$metadata, check_enc = TRUE))
if (nrow(a1) > nrow(a2)) {strain_data <- a1}else {strain_data <- a2}
rm(a1); rm(a2)

time2 <- suppressWarnings(readData(arg$tp2, check_enc = TRUE))
if (!exists("time2")) {time2 <- readData(arg$tp2, FALSE)}
rm(arg)

# Check that coefficient sets each add up to 1 -----------------------------------------------------------------
coeffset <- params$coeffs[2] %>% strsplit(",") %>% unlist()

for (i in 1:length(coeffset)) {
  x1 <- coeffset[i] %>% strsplit("-") %>% unlist() %>% as.double() %>% sum()
  assert(paste0("Parameter set ", i, " sums to 1"), x1 == 1)  
}

# Check that number of clusters (to generate heatmaps for) is a positive integer (or 0) ------------------------
assert("Number of clusters (for heatmap generation) is an integer >= 0", as.integer(params$numcl[2]) >= 0)

# LINEAGE INFO -------------------------------------------------------------------------------------------------

assert("Has lineage info", as.logical(params$has_lin[2]))
x <- updateStrains("lin_info", strain_data, time2)
strain_data <- x$sd; time2 <- x$t2
rm(x)

# COLUMN NAMES -------------------------------------------------------------------------------------------------
reqnames <- c("Strain", "Latitude", "Longitude", "Day", "Month", "Year")
cnames <- params$cnames[2] %>% strsplit(split = ",") %>% unlist()
if (nchar(params$cnames[2]) > 3) { # nchar([,]) == 3
  if (!("none" %in% cnames)) {
    fullcnames <- c(reqnames, cnames)
  }else {
    fullcnames <- reqnames  
  }
}else {
  fullcnames <- reqnames
}

strain_data <- strain_data %>% select(all_of(fullcnames)) %>% 
  na.omit(Strain) %>% na.omit(Latitude) %>% na.omit(Longitude) %>% na.omit(Day) %>% 
  na.omit(Month) %>% na.omit(Year)

# add column to show which strains are found in TP2
strain_data %<>% mutate(TP2 = ifelse(Strain %in% time2$Strain, 1, 0))

# REGION OF INTEREST -------------------------------------------------------------------------------------------
if (params$reg[2] != "All" & "Region" %in% colnames(strain_data)) {
  strain_data <- strain_data %>% filter(Region %in% params$reg[2])
}

# COUNTRY OF INTEREST ------------------------------------------------------------------------------------------
if (params$cou[2] != "All" & "Country" %in% colnames(strain_data)) {
  strain_data <- strain_data %>% filter(Country %in% params$cou[2])
}

# HAS DEFINED DATE INFO ----------------------------------------------------------------------------------------
if (as.logical(params$has_date[2])) {
  strain_data <- strain_data %>% 
    filter(!is.na(Day)) %>% filter(!is.na(Month)) %>% filter(!is.na(Year)) %>% 
    filter(Day != "") %>% filter(Month != "") %>% filter(Year != "")
}

if (nchar(params$temp_win[2]) > nchar("[,]")) {
  tempwindow <- params$temp_win[2] %>% gsub("\\[|\\]", "", .) %>% 
    strsplit(., ",") %>% unlist()
  
  strain_data %<>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
    filter(Date >= tempwindow[1] & Date <= tempwindow[2])
  
  y <- updateStrains("temp_win", strain_data, time2)
  strain_data <- y$sd; time2 <- y$t2
}

# HAS PROVINCE-LEVEL DATA --------------------------------------------------------------------------------------

if (as.logical(params$has_prov[2])) {
  strain_data <- strain_data %>% filter(!is.na(Province)) %>% filter(Province != "")
  
  if (params$prov[2] != "All") {
    strain_data <- strain_data %>% filter(Province %in% params$prov[2])
  }
  
  # Strains with metadata and defined lineage info at TP2
  z <- updateStrains("aft_prov", strain_data, time2)
  strain_data <- z$sd; time2 <- z$t2
}

# NON-SINGLETON CLUSTERS ---------------------------------------------------------------------------------------
outputDetails("Making table for matching TP2 clusters to integers (for metrics process) ...", TRUE)
processed_tp2 <- intClusters(time2); rm(time2)

th <- params$th[2]

if (as.logical(params$nsTP2[2])) {
  remove_strains <- strainsInSingletons(processed_tp2$new_cols, th)
  processed_tp2$new_cols %<>% filter(!(Strain %in% remove_strains))
  strain_data %<>% filter(!(Strain %in% remove_strains))
}

# INTERVAL TYPE ------------------------------------------------------------------------------------------------
save_to <- file.path(paste0("intermediate_data/cgms/", tolower(params$int_type[2])))
dir.create(save_to, showWarnings = FALSE)

metadata <- strain_data %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
  mutate(YearMonth = format(Date, "%Y-%m")) %>% mutate(Week = strftime(Date, format = "%V")) %>% 
  arrange(Week) %>% as.data.table()
rm(strain_data)

f2 <- processed_tp2$new_cols %>% as.data.table() %>% select(Strain, all_of(params$th[2])) %>% 
  rename(isolate = Strain) %>% arrange(isolate)

if (params$int_type[2] == "weekly") {
  interval <- "Week"
  interval_list <- metadata$Week %>% unique() %>% sort()
  
  interval_clusters <- metadata %>% select(c(Strain, all_of(interval))) %>% 
    inner_join(., f2, by = c("Strain" = "isolate")) %>% 
    set_colnames(c("isolate", "ivl", "heightx"))
  
}else if (params$int_type[2] == "monthly") {
  interval <- "YearMonth"
  interval_list <- metadata %>% pull(interval) %>% unique() %>% sort()
  
  interval_clusters <- metadata %>% select(c(Strain, all_of(interval))) %>% 
    inner_join(., f2, by = c("Strain" = "isolate")) %>% 
    set_colnames(c("isolate", "ivl", "heightx"))
  
  
}else if (params$int_type[2] == "multiset") {
  setdivider <- strsplit(params$divs[2], split = ",") %>% unlist() %>% as.Date(., format = "%Y-%m-%d")
  
  assertion1 <- lapply(setdivider, function(x) nchar(as.character(x)) == 10) %>% unlist()
  assert("Provided input(s) has/have 10 characters", all(assertion1))
  assert("Provided input(s) is/are date type, correct format", !any(is.na(setdivider)))
  assertion3 <- lapply(setdivider, function(x) x >= min(metadata$Date) & x <= max(metadata$Date)) %>% unlist()
  assert("Provided input(s) is/are between min date and max date", all(assertion3))
  
  sdiv <- c(min(metadata$Date), setdivider, max(metadata$Date))
  fdivs <- tibble(start = sdiv[1:(length(sdiv)-1)], end = sdiv[2:length(sdiv)]) %>% 
    mutate(ivl = paste0("set", 1:nrow(.)))
  
  metadata <- metadata %>% add_column(Multiset = "")
  for (i in 1:nrow(fdivs)) {
    metadata[metadata$Date >= fdivs$start[i] & metadata$Date < fdivs$end[i]]$Multiset <- fdivs$ivl[i]
  }
  metadata[metadata$Date == fdivs$end[i]]$Multiset <- fdivs$ivl[i]
  rm(fdivs)
  
  interval <- "Multiset"
  interval_list <- metadata %>% pull(interval) %>% unique() %>% sort()
  
  interval_clusters <- metadata %>% select(c(Strain, Multiset)) %>%
    inner_join(., f2, by = c("Strain" = "isolate")) %>%
    set_colnames(c("isolate", "ivl", "heightx"))
}
rm(f2)

clustersets <- vector(mode = "list", length = length(interval_list)) %>% set_names(interval_list)

for (xj in interval_list) {
  # cluster assignments for clusters that changed when interval i strains were added
  int_j <- interval_clusters[heightx %in% interval_clusters[ivl == xj]$heightx]
  sofar <- interval_clusters[heightx %in% interval_clusters[ivl <= xj]$heightx]
  clustersets[[xj]] <- list(int_j, sofar) %>% set_names(c("ivl", "sofar"))
}
rm(sofar); rm(int_j); rm(interval_clusters); rm(interval_list)

# SAVING RESULTS -----------------------------------------------------------------------------------------------
saveRDS(clustersets, file.path("inputs", "processed", "clustersets.Rds")); rm(clustersets)
writeData(metadata, file.path("inputs", "processed", "strain_info.txt")); rm(metadata)
writeData(processed_tp2$new_cols, file.path("inputs", "processed", "tp2_clusters.txt"))
saveRDS(processed_tp2, file.path("inputs", "processed", "allTP2.Rds")); rm(processed_tp2)

outputDetails(paste0("\nFinished process at: ", Sys.time(), "\n||", paste0(rep("-", 14), collapse = ""), 
                     " Saved formatted inputs to 'processed' in the ",
                     "inputs", " directory", paste0(rep("-", 15), collapse = ""), "||"), TRUE)
