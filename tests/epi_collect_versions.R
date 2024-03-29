
libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "Rcpp", "RcppArmadillo")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

sourceCpp("scripts/epicohversions.cpp")

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
# "scripts/ECC/ecc_functions2.R")
invisible(sapply(files, source))

assert("Distances were collected and saved", file.exists("intermediate_data/TPN/extreme_dists.Rds"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ",
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-t", "--trio"), metavar = "character", default = "0-1-0",
              help = "source, temporal, geographic coefficients"), 
  make_option(c("-i", "--intervaltype"), metavar = "char", default = readLines("scripts/date.txt")[1], 
              help = "Type of intervals, choices are: weekly, monthly, multiset. If multiset, provide a time to split the dataset at."))

params <- parse_args(OptionParser(option_list=option_list))

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))
m <- read_tsv(params$metadata) %>% processedStrains()

metadata <- m$strain_data %>% 
  mutate(YearMonth = format(Date, "%Y-%m")) %>% 
  mutate(Week = strftime(Date, format = "%V")) %>% 
  select(-TP1, -TP2) %>% arrange(Week) %>% as.data.table()

tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
f2 <- tp2$filedata %>% rownames_to_column("isolate") %>% as.data.table() %>% 
  select(isolate, all_of(params$heights)) %>% arrange(isolate)

extremes <- readRDS("intermediate_data/TPN/extreme_dists.Rds")

source("scripts/interval_prep.R")

clusters <- vector(mode = "list", length = length(interval_list)) %>% set_names(interval_list)

for (xj in interval_list) {
  # cluster assignments for clusters that changed when interval i strains were added
  int_j <- interval_clusters[heightx %in% interval_clusters[ivl == xj]$heightx]
  sofar <- interval_clusters[heightx %in% interval_clusters[ivl <= xj]$heightx]
  clusters[[xj]] <- list(int_j, sofar) %>% set_names(c("ivl", "sofar"))
}

typing_data <- lapply(1:length(interval_list), function(i) {
  n1 <- as.character(interval_list[i])
  tpkstrains <- metadata[get(interval) <= n1]$Strain
  dfz <- tp2$filedata %>% rownames_to_column("isolate") %>%
    select(isolate, all_of(params$heights)) %>%
    filter(isolate %in% tpkstrains) %>% column_to_rownames("isolate")
  dfz[,hx$h[1],drop=FALSE] %>% set_colnames(hx$th[1])
}) %>% set_names(as.character(interval_list))

td <- typing_data[[length(typing_data)]] %>% rownames_to_column("Strain") %>% as.data.table()
parts <- m$dr_matches %>% filter(Strain %in% td$Strain) %>% 
  left_join(td, ., by = "Strain") %>% sectionClusters(.)
df <- parts$drs
cx <- setdiff(colnames(df), c("Strain", "dr"))
results <- parts$results

dfx <- params$trio %>% strsplit(., ",") %>% unlist() %>% 
  expand.grid(x = ., k = interval_list) %>% as.data.frame()

i <- 2

k <- as.character(dfx$k[i])
tpkstrains <- metadata[get(interval) <= k]$Strain
key_cls <- parts$drs[Strain %in% tpkstrains] %>% select(-Strain, -dr) %>% pull() %>% unique()
y <- lapply(parts$results, function(x) any(key_cls %in% pull(x, 1))) %>% unlist()

k_drs <- m$dr_matches %>% filter(Strain %in% tpkstrains) %>% pull(dr)

c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "-") %>% unlist() %>% as.numeric()
tau <- c1[2]
gamma <- c1[3]

fnames <- names(y[y])

eccresults <- lapply(fnames, function(f) {
  cluster_x <- df[df[[cx]] %in% pull(results[[f]], cx),-"Strain"]
  dms <- readRDS(paste0("intermediate_data/TP", k, "/dists/group", f, ".Rds"))
  selected_tp <- m$strain_data %>% filter(Strain %in% tpkstrains)
  
  tr_dists <- collectTransforms2(dms, extremes)
  eccnew <- epiCollectionByClusterV2(selected_tp, tau, gamma, tr_dists, k, cluster_x[dr %in% k_drs])
  
  transformed_dists <- collectTransforms(dms, extremes)
  eccold <- epiCollectionByClusterV1(selected_tp, tau, gamma, transformed_dists, k, cluster_x[dr %in% k_drs])
  
  identical(eccnew, eccold)
}) %>% unlist()

assert("V4 and V3 have same results", all(eccresults))  
