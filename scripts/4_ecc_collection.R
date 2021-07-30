#! /usr/bin/env Rscript

msg <- file("logs/logfile_epiquant.txt", open="wt")
sink(msg, type="message")

source("scripts/ECC/ecc_opener.R")
source("scripts/ECC/dist_functions.R")

assert("Distances were collected and saved", file.exists("intermediate_data/dist_extremes.Rds"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ",
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))

extremes <- readRDS("intermediate_data/dist_extremes.Rds")

dfx <- expand.grid(x = combos, k = c(1,2)) %>% as.data.frame()

for (k in unique(dfx$k)) {
  dir.create(paste0("intermediate_data/TP", k, "/ecc_groups/"), showWarnings = FALSE)
}

for (i in 1:nrow(dfx)) {
  cat(paste0("\n\nStep ", i, " / ", nrow(dfx) + 2, ":"))
  c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "-") %>% 
    unlist() %>% as.numeric()
  
  paste0("intermediate_data/TP", dfx$k[i], "/ecc_groups/", as.character(dfx$x[i])) %>% 
    dir.create(., showWarnings = FALSE)
  
  outputMessages(paste0("\nCollecting and saving ECCs for groups of clusters at TP", 
                        dfx$k[i], ", for ECC coefficient triple ", as.character(dfx$x[i])))
    
  k <- dfx$k[i]
  parts <- sectionClusters(k, typing_data, m)
  read_from <- paste0("intermediate_data/TP", k, "/dists/")
  save_to <- paste0("intermediate_data/TP", k, "/ecc_groups/", as.character(dfx$x[i]), "/")
  
  collectECCs(k, m, parts, extremes, c1, read_from, save_to)
  
  c1 %>% paste0(., collapse = ", ") %>% 
    paste0("\nMerging ECC files at TP", dfx$k[i], ", for ECC parameters ", .) %>% 
    outputMessages()
}

cat(paste0("\n\nStep ", nrow(dfx) + 1, " / ", nrow(dfx) + 2, ":"))
for (i in 1:nrow(dfx)) {
  fnames <- list.files(paste0("intermediate_data/TP", dfx$k[i], 
                              "/ecc_groups/", dfx$x[i], "/"), full.names = TRUE)
  
  tpk <- lapply(fnames, function(f) {readRDS(f)}) %>% bind_rows()
  
  saveRDS(tpk, paste0("intermediate_data/TP", dfx$k[i], "/", dfx$x[i], "-eccs.Rds"))
}

# Generating ECC results file ----------------------------------------------------
outputMessages("   Putting together saved ECCs, replacing blanks with NAs")
a1 <- lapply(combos, function(x) {
  readRDS(paste0("intermediate_data/TP1/", x, "-eccs.Rds"))
}) %>% Reduce(inner_join, .) %>% unique()

b1 <- lapply(combos, function(x) {
  readRDS(paste0("intermediate_data/TP2/", x, "-eccs.Rds"))
}) %>% Reduce(inner_join, .) %>% unique()

tp1_eccs <- tp1$proc[,c("Strain",colnames(a1)[1])] %>% left_join(., a1)
tp2_eccs <- tp2$proc[,c("Strain",colnames(b1)[1])] %>% left_join(., b1)

all_eccs <- right_join(tp1_eccs, tp2_eccs, by = "Strain")

# The -Infs are because of TP1 singletons
tp1size <- colnames(all_eccs) %>% 
  grep("Size", ., value = TRUE) %>% grep("TP1", ., value = TRUE)
inf_inds <- which(all_eccs == -Inf, arr.ind = TRUE) %>% as.data.frame() %>% pull("row")
num_others <- all_eccs[inf_inds,] %>% filter(!is.na(!!as.symbol(tp1size))) %>% 
  pull(tp1size) %>% unique() %>% setdiff(., 1) %>% length()
assert("The -Infs are because of TP1 singletons", num_others == 0)
all_eccs[all_eccs == -Inf] <- NA

assert("Average distances were collected and saved", 
       file.exists("intermediate_data/average_dists.Rds"))
all_avg_dists <- readRDS("intermediate_data/average_dists.Rds")

cat(paste0("\n\nStep ", nrow(dfx) + 2, " / ", nrow(dfx) + 2, ":"))
outputMessages("   Merging collected ECCs ...\n")
inner_join(all_avg_dists, all_eccs) %>% 
  write.table(., file = "results/ECCs.tsv", col.names = TRUE, 
              row.names = FALSE, quote = FALSE, sep = "\t")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ",
           paste0(rep("-", 31), collapse = ""), "||\n"))
