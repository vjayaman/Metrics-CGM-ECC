#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/dist_matrices.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y))); rm(y); rm(libs)

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source)); rm(files)

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", 
              help = "Time point 2 file name (TP2)"), 
  make_option(c("-f", "--intervalfile"), metavar = "file", default = "inputs/processed/clustersets.Rds"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))
arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)

cat(paste0("\n||", paste0(rep("-", 23), collapse = ""), 
           " Generating non-redundant pairwise distances ", 
           paste0(rep("-", 23), collapse = ""), "||\nStarted process at: ", Sys.time()))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

# COLLECT dist matrices using TPN clusters ---------------------------------------------------------
outputMessages(paste0("\nCollecting and saving distances for groups at TPN"))

params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

hx <- strsplit(as.character(params$th[2]), split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))
tp2 <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)
m <- read_tsv(arg$metadata) %>% processedStrains()
basedir <- "intermediate_data/TPN/"

# Extremes, for scaling ----------------------------------------------------------------------
ext_geo_dists <- m$assignments %>% select(Longitude, Latitude) %>% unique() %>% 
  rownames_to_column("id") %>% distMatrix(., "geo", c("Longitude", "Latitude"))
ext_temp_dists <- m$assignments %>% select(Date) %>% unique() %>% 
  rownames_to_column("id") %>% distMatrix(., "temp", "Date")

extremes <- list(maxt = max(ext_temp_dists), mint = min(ext_temp_dists), 
                 maxg = max(ext_geo_dists), ming = min(ext_geo_dists))
saveRDS(extremes, paste0(basedir, "extreme_dists.Rds"))
rm(ext_geo_dists); rm(ext_temp_dists); rm(extremes)

# Pairwise distances within clusters at TPN --------------------------------------------------
metadata <- m$strain_data %>% as.data.table()

clustersets <- readRDS(arg$intervalfile)
interval_list <- names(clustersets)
rm(clustersets)

k <- last(interval_list)
save_to <- paste0("intermediate_data/TPN/dists/")
dir.create(save_to,  recursive = TRUE, showWarnings = FALSE)

if (params$int_type[2] == "multiset") {
  interval <- "Multiset"
}else if (params$int_type[2] == "monthly") {
  interval <- "YearMonth"
}else if (params$int_type[2] == "weekly") {
  interval <- "Week"
}

typing_data <- lapply(1:length(interval_list), function(i) {
  n1 <- as.character(interval_list[i])
  tpkstrains <- metadata[get(interval) <= n1]$Strain
  dfz <- tp2$filedata %>% rownames_to_column("isolate") %>%
    select(isolate, all_of(hx$h)) %>%
    filter(isolate %in% tpkstrains) %>% column_to_rownames("isolate")
  dfz[,hx$h[1],drop=FALSE] %>% set_colnames(hx$th[1])
}) %>% set_names(as.character(interval_list))

td <- typing_data[[length(typing_data)]] %>% rownames_to_column("Strain") %>% as.data.table()
rm(typing_data)

parts <- m$dr_matches %>% filter(Strain %in% td$Strain) %>% 
  left_join(td, ., by = "Strain") %>% sectionClusters(.)
saveRDS(parts, "intermediate_data/TPN/parts.Rds")

outputMessages(paste0("Collecting and saving distances for cluster groups at TP", k, ":\n"))

tpkstrains <- metadata[get(interval) <= k]$Strain
collectDistances(parts$drs, parts$results, m$dr_matches, m$assignments, tpkstrains, save_to)
rm(m); rm(parts)

assert("Distances were collected and saved", file.exists(paste0(basedir, "extreme_dists.Rds")))

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
timeTaken(pt = "distances collection", stopwatch) %>% outputDetails(., newcat = TRUE)
cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), 
           " End of distances collection ", paste0(rep("-", 31), collapse = ""), "||\n"))