#! /usr/bin/env Rscript

msg <- file("logs/merging_ecc.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source))

assert("Distances were collected and saved", file.exists("intermediate_data/dist_extremes.Rds"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ",
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))

option_list <- list(
  make_option(c("-m", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-a", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-t", "--trio"), metavar = "character", default = "0-1-0",
              help = "source, temporal, geographic coefficients"))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$trio %>% strsplit(., ",") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

dfx <- expand.grid(x = combos, k = c(1,2)) %>% as.data.frame()

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