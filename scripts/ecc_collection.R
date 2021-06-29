#! /usr/bin/env Rscript

msg <- file("logs/logfile_epiquant.txt", open="wt")
sink(msg, type="message")

source("scripts/ecc_opener.R")
assert("Distances were collected and saved", file.exists("results/TP1/dists/group1.Rds"))

# TIMEPOINT 2 ANALYSIS -----------------------------------------------------------------------
extremes <- readRDS("results/dist_extremes.Rds")

dfx <- expand.grid(x = combos, k = c(1,2)) %>% as.data.frame()

for (k in unique(dfx$k)) {dir.create(paste0("results/TP", k, "/eccs"), showWarnings = FALSE)}

for (i in 1:nrow(dfx)) {
  c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "") %>% unlist() %>% as.numeric()
  
  dir.create(paste0("results/TP", dfx$k[i], "/eccs/", as.character(dfx$x[i])), showWarnings = FALSE)
  
  outputMessages(paste0("\nCollecting and saving ECCs for groups of clusters at TP", dfx$k[i], 
                 ", for ecc triple ", as.character(dfx$x[i])))
  sectionClusters(dfx$k[i], typing_data, m) %>% 
    collectECCs(dfx$k[i], m, ., extremes, c1, 
                paste0("results/TP", dfx$k[i], "/dists/"), 
                paste0("results/TP", dfx$k[i], "/eccs/", as.character(dfx$x[i]), "/"))
  
  outputMessages(paste0("\nMerging ECC files at TP", dfx$k[i]))
}


for (i in 1:nrow(dfx)) {
  dir_k <- paste0("results/TP", dfx$k[i], "/eccs/", dfx$x[i], "/")
  fnames <- list.files(dir_k)
  
  tpk <- lapply(fnames, function(f) {
    ecc_j <- readRDS(paste0(dir_k, f))
    file.remove(paste0(dir_k, f))
    return(ecc_j)
  }) %>% bind_rows()
  
  saveRDS(tpk, paste0(dir_k, "TP", dfx$k[i], "-", dfx$x[i], "-eccs.Rds"))
}

# ---------------------------------------------------------------------------------------------
# Generating ECC results file

a1 <- lapply(combos, function(x) {
  readRDS(paste0("results/TP1/eccs/", x, "/TP1-", x, "-eccs.Rds"))
}) %>% Reduce(inner_join, .)

b1 <- lapply(combos, function(x) {
  readRDS(paste0("results/TP2/eccs/", x, "/TP2-", x, "-eccs.Rds"))
}) %>% Reduce(inner_join, .)

tp1_eccs <- tp1$proc[,c("Strain",colnames(a1)[1])] %>% left_join(., a1)
tp2_eccs <- tp2$proc[,c("Strain",colnames(b1)[1])] %>% left_join(., b1)

all_eccs <- right_join(tp1_eccs, tp2_eccs, by = "Strain")

# The -Infs are because of TP1 singletons
tp1size <- grep("Size", colnames(all_eccs), value = TRUE) %>% grep("TP1", ., value = TRUE)
inf_inds <- which(all_eccs == -Inf, arr.ind = TRUE) %>% as.data.frame() %>% pull("row")
num_others <- all_eccs[inf_inds,] %>% filter(!is.na(!!as.symbol(tp1size))) %>% 
  pull(tp1size) %>% unique() %>% setdiff(., 1) %>% length()
assert("The -Infs are because of TP1 singletons", num_others == 0)
all_eccs[all_eccs == -Inf] <- NA

cat(paste0("\n\nStep ", length(combos) + 2, ":"))
outputMessages("   Merging collected ECCs ...\n")
write.table(all_eccs, file = "results/ECCs.tsv", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ",
           paste0(rep("-", 31), collapse = ""), "||\n"))
