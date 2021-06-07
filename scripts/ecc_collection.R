#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_epiquant.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

files <- paste0("scripts/ECC") %>% list.files(., full.names = TRUE)
invisible(sapply(files, source))

# Title: "EpiQuant - Salmonella Enteritidis Project (2019-2020)"
# Authors of original work and initial modifications: Ben Hetman, Elissa Giang, Dillon Barker
# Responsible for changes made during and after merge with CGM process: Vasena Jayamanna

option_list <- list(
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"),
  make_option(c("-t", "--trio"), metavar = "character", default = "010-001",
              help = "temporal, geographic coefficients"))

cat(paste0("\n||", hyps(34), " ECC metric generation ", hyps(34), "||\nStarted process at: ", Sys.time()))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))
combos <- params$trio %>% strsplit(., "-") %>% unlist()

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

strain_data <- read_tsv(params$strains) %>% 
  mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")))

cat(paste0("\n\nStep 1:\n   Note that the source coefficent is always 0 in this version"))

loc_cols <- length(intersect(c("Country", "Province", "City"), colnames(strain_data)))

# Note: dr stands for data representative
# in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
outputMessages("   Removing redundancy (comparing date, lat-long, etc. - not every pair of strains)")

outputMessages("   Identifying which strains match which non-redundant 'data representatives'")

assignments <- strain_data %>% select(Date, Latitude, Longitude) %>% 
  unique() %>% rownames_to_column("dr")

dr_matches <- assignments %>% 
  left_join(strain_data, ., by = byNames(strain_data, .)) %>% select(Strain, dr)

dists <- nonRedundantDists(assignments)

avgdistvals <- lapply(1:length(typing_data), function(i) {
  g_cuts <- drCounts(typing_data[[i]], dr_matches)
  
  outputMessages(paste0("      Calculating average (not transformed) distances for timepoint ", i))
  a2 <- avgDists(g_cuts, dists$temp, "Temp.Dist", paste0("TP", i, "_", colnames(g_cuts)[1]))
  b2 <- avgDists(g_cuts, dists$geo, "Geog.Dist", paste0("TP", i, "_", colnames(g_cuts)[1]))
  return(list(temp = a2, geo = b2))
}) %>% set_names(c("TP1", "TP2"))

collected_eccs <- lapply(1:length(combos), function(j) {
  c1 <- unlist(strsplit(combos[j], split = "")) %>% as.numeric()
  cat(paste0("\n\nStep ", j + 1, ":"))
  tau <- c1[2]
  gamma <- c1[3]
  epitable <- prepEpiTable(dists$tr_dists, tau, gamma)
  sim_matrix <- prepEpiMatrix(epitable)
  x1 <- epiCollection(sim_matrix, typing_data, tau, gamma, dr_matches, avgdistvals)
})

cat(paste0("\n\nStep ", length(combos) + 2, ":"))
outputMessages("   Merging collected ECCs ...\n")
full_set <- mergeECCs(collected_eccs, 1, tp1$proc) %>%
  merge.data.table(., mergeECCs(collected_eccs, 2, tp2$proc), by = "Strain", all.y = TRUE) %>%
  mutate(TP1 = ifelse(is.na(TP1), 0, TP1))

write.table(full_set, file = "results/ECCs.tsv", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", hyps(30), " End of ECC metric generation ", hyps(31), "||\n"))
