#! /usr/bin/env Rscript

msg <- file("logs/logfile_epiquant1.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "Rcpp", "RcppArmadillo")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y))); rm(libs); rm(y)

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)
sourceCpp("scripts/epicohversions.cpp")

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source)); rm(files)

assert("Distances were collected and saved", file.exists("intermediate_data/TPN/extreme_dists.Rds"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ",
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", 
              help = "Time point 2 file name (TP2)"), 
  make_option(c("f", "--intervalfile"), metavar = "file", default = "inputs/processed/clustersets.Rds"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))
# option_list <- list(
#   make_option(c("-m", "--metadata"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
#   make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
#   make_option(c("-x", "--heights"), metavar = "character", default = "0",
#               help = "Comma-delimited string of heights to collect ECCs for"),
#   make_option(c("-t", "--trio"), metavar = "character", default = "0-1-0",
#               help = "source, temporal, geographic coefficients"), 
#   make_option(c("-i", "--intervaltype"), metavar = "char", default = readLines("scripts/date.txt")[1], 
#               help = "Type of intervals, choices are: weekly, monthly, multiset. If multiset, provide a time to split the dataset at."))

arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)

params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

hx <- strsplit(as.character(params$th[2]), split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))
tp2 <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)
m <- read_tsv(arg$metadata) %>% processedStrains()
metadata <- m$strain_data %>% as.data.table()
extremes <- readRDS("intermediate_data/TPN/extreme_dists.Rds")

if (params$int_type[2] == "multiset") {
  interval <- "Multiset"
}else if (params$int_type[2] == "monthly") {
  interval <- "YearMonth"
}else if (params$int_type[2] == "weekly") {
  interval <- "Week"
}

# source("scripts/interval_prep.R")
clustersets <- readRDS(arg$intervalfile)
interval_list <- names(clustersets)
rm(clustersets)

for (xj in interval_list) {
  paste0("intermediate_data/TP", xj, "/eccs/") %>% dir.create(., showWarnings = FALSE, recursive = TRUE)
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
df <- parts$drs
cx <- setdiff(colnames(df), c("Strain", "dr"))
results <- parts$results

dfx <- params$coeffs[2] %>% strsplit(., ",") %>% unlist() %>% 
  expand.grid(x = ., k = interval_list) %>% as.data.frame()

k <- last(interval_list)
tpkstrains <- metadata[get(interval) <= k]$Strain
key_cls <- parts$drs[Strain %in% tpkstrains] %>% select(-Strain, -dr) %>% pull() %>% unique()
y <- lapply(parts$results, function(x) any(key_cls %in% pull(x, 1))) %>% unlist()
fnames <- names(y[y])

for (j in 1:length(fnames)) {
  cat(paste0("\nStep ", j, " / ", length(fnames), ":\n"))
  f <- fnames[j]
  dms <- readRDS(paste0("intermediate_data/TP", as.character(last(dfx$k)), "/dists/group", f, ".Rds"))
  tr_dists <- collectTransforms2(dms, extremes)
  
  qb <- txtProgressBar(min = 0, max = nrow(dfx), initial = 0, style = 3)
  for (i in 1:nrow(dfx)) {
    
    k <- as.character(dfx$k[i])
    tpkstrains <- metadata[get(interval) <= k]$Strain
    key_cls <- parts$drs[Strain %in% tpkstrains] %>% select(-Strain, -dr) %>% pull() %>% unique()
    k_drs <- m$dr_matches %>% filter(Strain %in% tpkstrains) %>% pull(dr)
    
    c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "-") %>% unlist() %>% as.numeric()
    tau <- c1[2]
    gamma <- c1[3]
    
    cluster_x <- df[df[[cx]] %in% pull(results[[f]], cx),-"Strain"]
    selected_tp <- m$strain_data %>% filter(Strain %in% tpkstrains)
    
    epiCollectionByClusterV2(selected_tp, tau, gamma, tr_dists, k, cluster_x[dr %in% k_drs]) %>% 
      saveRDS(., paste0("intermediate_data/TP", k, "/eccs/group", f, ".Rds"))
    
    setTxtProgressBar(qb, i)
  }
  close(qb)
}

ecc_results <- lapply(1:nrow(dfx), function(j) {
  fnames <- list.files(paste0("intermediate_data/TP", dfx$k[j], "/eccs"), full.names = TRUE)
  tpk <- lapply(fnames, function(f) {readRDS(f)}) %>% bind_rows() %>% 
    mutate(across(.cols = everything(), as.double))
  
  cname <- strsplit(colnames(tpk), split = "_") %>% sapply(., "[[", 1) %>% unique()
  colnames(tpk) <- gsub(paste0(cname, "_"), "", colnames(tpk))
  
  tpk %>% add_column(TP = gsub("TP", "", cname), .before = 1)
}) %>% bind_rows()

saveRDS(ecc_results, "results/ecc_results.Rds")

# Generating ECC results file ----------------------------------------------------
assert("No -Inf ECC results", !any(is.infinite(abs(pull(ecc_results[,4])))))

# assert("Average distances were collected and saved", 
#        file.exists("intermediate_data/average_dists.Rds"))
# all_avg_dists <- readRDS("intermediate_data/average_dists.Rds")
# 
# cat(paste0("\n\nStep ", nrow(dfx) + 2, " / ", nrow(dfx) + 2, ":"))
# outputMessages("   Merging collected ECCs ...\n")
# inner_join(all_avg_dists, all_eccs) %>% 
#   write.table(., file = "results/ECCs.tsv", col.names = TRUE, 
#               row.names = FALSE, quote = FALSE, sep = "\t")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ",
           paste0(rep("-", 31), collapse = ""), "||\n"))
