#! /usr/bin/env Rscript

msg <- file("logs/logfile_alldists.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "progress")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", 
           "scripts/ECC/dist_functions.R")
invisible(sapply(files, source))

option_list <- list(
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"), 
  make_option(c("-t", "--timepoint"), metavar = "numeric", default = 2))

cat(paste0("\n||", paste0(rep("-", 22), collapse = ""), 
           " Saving full set of (strain) pairwise distances ", 
           paste0(rep("-", 21), collapse = ""), "||\nStarted process at: ", Sys.time()))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

m <- read_tsv(params$strains) %>% processedStrains()

k <- params$timepoint

intra_cl <- paste0("intermediate_data/TP", k, "/dists/") %>% 
  list.files(., full.names = TRUE, pattern = "group") %>% 
  lapply(., function(fname) {readRDS(fname)})

parts <- sectionClusters(k, typing_data, m)
rm(tp1)
rm(tp2)

dir.create(paste0("results/", "dists"), showWarnings = FALSE)
dir.create(paste0("results/dists/", "TP", k), showWarnings = FALSE)

outputMessages(paste0("\nFormatting and merging distances matrix for each cluster ...\n\n"))

intra_dists <- lapply(1:length(intra_cl), function(i) {
  a1 <- intra_cl[[i]][["temp"]] %>% as.data.frame() %>% rownames_to_column("dr") %>% 
    as.data.table() %>% melt.data.table(id.vars = "dr") %>% 
    set_colnames(c("dr1", "dr2", "Temp.Dist"))
  
  a2 <- intra_cl[[i]][["geo"]] %>% as.data.frame() %>% rownames_to_column("dr") %>% 
    as.data.table() %>% melt.data.table(id.vars = "dr") %>% 
    set_colnames(c("dr1", "dr2", "Geog.Dist"))
  
  merge.data.table(a1, a2) %>% return()
}) %>% bind_rows() %>% unique() %>% 
  mutate(across(dr2, as.character)) %>% as.data.table()

rm(intra_cl)
gc()
outputMessages("\n")

dr_matches <- parts$drs %>% as.data.table() %>% set_colnames(c("Strain", "TP_cl", "dr"))
clusters <- dr_matches %>% select(TP_cl) %>% unique() %>% pull()

pb <- txtProgressBar(min = 0, max = length(clusters), initial = 0, style = 3)

for (i in 1:length(clusters)) {
  x <- clusters[i]
  cl1 <- dr_matches[TP_cl %in% x]

  cl1_dists <- inner_join(cl1, intra_dists, by = c("dr"="dr1")) %>%
    select(-TP_cl, -dr) %>% rename(Strain1 = Strain) %>%
    inner_join(cl1, ., by = c("dr" = "dr2")) %>%
    select(-TP_cl, -dr) %>% rename(Strain2 = Strain)

  saveRDS(cl1_dists, paste0("results/dists/TP", k, "/cluster", 
                            formatC(x, width = nchar(max(clusters)), format = "d", flag = "0"), 
                            ".Rds"))
  setTxtProgressBar(pb, x)
}
close(pb)

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())

outputDetails(paste0("\nSuccessfully collected data for all heights."), newcat = TRUE)
timeTaken(pt = "Distance data collection", stopwatch) %>% outputDetails(., newcat = TRUE)

cat(paste0("\n||", paste0(rep("-", 27), collapse = ""),
           " End of processing full distances set ",
           paste0(rep("-", 26), collapse = ""), "||\n"))
