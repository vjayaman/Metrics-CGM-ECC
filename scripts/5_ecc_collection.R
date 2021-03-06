#! /usr/bin/env Rscript

msg <- file("logs/ecc_collection.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "Rcpp", "RcppArmadillo")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y))); rm(libs); rm(y)

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)
sourceCpp("scripts/epicohversions.cpp")

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", 
           "scripts/ECC/ecc_functions.R", 
           "scripts/ECC/dist_functions.R")
invisible(sapply(files, source)); rm(files)

assert("Distances were collected and saved", file.exists("intermediate_data/TPN/extreme_dists.Rds"))

outputDetails(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ",
                     paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", 
              help = "Time point 2 file name (TP2)"), 
  make_option(c("-f", "--intervalfile"), metavar = "file", default = "inputs/processed/clustersets.Rds"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))

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

clustersets <- readRDS(arg$intervalfile)
interval_list <- names(clustersets)
rm(clustersets)

basedir <- file.path("intermediate_data", params$int_type[2], "eccs")
for (xj in interval_list) {
  file.path(basedir, paste0("TP", xj)) %>% 
    dir.create(., showWarnings = FALSE, recursive = TRUE)
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

parts <- m$dr_matches %>% filter(Strain %in% td$Strain) %>% 
  left_join(td, ., by = "Strain") %>% sectionClusters(.)
df <- parts$drs
cx <- setdiff(colnames(df), c("Strain", "dr"))
results <- parts$results

dfx <- params$coeffs[2] %>% strsplit(., ",") %>% unlist() %>% 
  expand.grid(x = ., k = interval_list) %>% as.data.frame() %>% as.data.table()

k <- last(interval_list)
tpkstrains <- metadata[get(interval) <= k]$Strain
key_cls <- parts$drs[Strain %in% tpkstrains] %>% select(-Strain, -dr) %>% pull() %>% unique()
y <- lapply(parts$results, function(x) any(key_cls %in% pull(x, 1))) %>% unlist()
fnames <- names(y[y]) # groups of clusters, the pairwise distances for each group

m$dr_matches <- m$dr_matches %>% as.data.table()

dates <- as.character(unique(dfx$k))

for (j in 1:length(fnames)) {
  f <- fnames[j]
  outputDetails(paste0("\n\nGroup of clusters ", f, ", ", j, " / ", length(fnames), " ... "))
  outputDetails(paste0("\n     TP             tau             gamma"))
  dms <- readRDS(paste0("intermediate_data/TPN/dists/group", f, ".Rds"))
  tr_dists <- collectTransforms2(dms, extremes)
  
  for (k_i in dates) {
    # print(k_i)
    tpkstrains <- metadata[get(interval) <= k_i]$Strain
    key_cls <- parts$drs[Strain %in% tpkstrains] %>% select(-Strain, -dr) %>% pull() %>% unique()
    k_drs <- m$dr_matches[Strain %in% tpkstrains] %>% pull(dr)
    cluster_x <- df[df[[cx]] %in% pull(results[[f]], cx),-"Strain"]
    selected_tp <- m$strain_data %>% filter(Strain %in% tpkstrains)
    
    eccs <- lapply(as.character(unique(dfx$x)), function(x_i) {
      coeffs <- unlist(strsplit(x_i, split = "-"))
      tau <- coeffs[2]
      gamma <- coeffs[3]
      outputDetails(paste0("\n   ", k_i, "           ", tau, "                ", gamma))
      cluster_y <- cluster_x[dr %in% k_drs]
      if (nrow(cluster_y) > 0) {
        epiCollectionByClusterV2(selected_tp, tau, gamma, tr_dists, k_i, cluster_y)  
      }
    })
    
    if (length(eccs[!sapply(eccs, is.null)]) > 0) {
      save_to <- file.path(basedir, paste0("TP", k_i), paste0("group", f, ".Rds"))
      suppressMessages(Reduce(inner_join, eccs)) %>% saveRDS(., save_to)  
    }
  }
}

datafiles <- lapply(interval_list, function(tp) {
  data.table(tp, fname = list.files(file.path(basedir, paste0("TP", tp)), full.names = TRUE))
}) %>% bind_rows() %>% arrange(tp)

ecc_results <- lapply(1:nrow(datafiles), function(i) {
  dfy <- readRDS(datafiles$fname[i]) %>% mutate(across(.cols = everything(), as.double))
  colnames(dfy) <- gsub(paste0("TP", datafiles$tp[i], "_"), "", colnames(dfy))
  dfy %>% add_column(TP = datafiles$tp[i], .before = 1)
}) %>% bind_rows()

if (params$int_type[2] == "multiset") {
  res_file <- gsub("-", "", params$divs[2]) %>% gsub(",", "-", .) %>% 
    paste0("results/ECC-", ., "-midpoints.Rds")  
}else {
  res_file <- paste0("results/ECC-", params$int_type[2], "-intervals.Rds")
}

assert("No -Inf ECC results", !any(is.infinite(abs(pull(ecc_results[,4])))))
saveRDS(ecc_results, res_file)

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
timeTaken(pt = "ECC data collection", stopwatch) %>% outputDetails(., newcat = TRUE)
cat(paste0("||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ",
           paste0(rep("-", 31), collapse = ""), "||\n"))
