#! /usr/bin/env Rscript

# msg <- file("logs/logfile_epiquant.txt", open="wt")
# sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- paste0("scripts/ECC") %>% list.files(., full.names = TRUE)
invisible(sapply(files, source))

# Title: "EpiQuant - Salmonella Enteritidis Project (2019-2020)"
# Authors of original work and initial modifications: Ben Hetman, Elissa Giang, Dillon Barker
# Responsible for changes made during and after merge with CGM process: Vasena Jayamanna

option_list <- list(
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"),
  make_option(c("-t", "--trio"), metavar = "character", default = "010",
              help = "temporal, geographic coefficients"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ", 
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$trio %>% strsplit(., "-") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

m <- read_tsv(params$strains) %>% processedStrains()

cat(paste0("\n\nStep 1:"))
cat(paste0("\n   Note that the source coefficent is always 0 in this version"))

# Note: dr stands for data representative
# in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
# outputMessages("   Removing redundancy (comparing date, lat-long, etc. - not every pair of strains)")
# outputMessages("   Identifying which strains match which non-redundant 'data representatives'")


# TIMEPOINT 2 ANALYSIS -----------------------------------------------------------------------

k <- 2
df <- typing_data[[k]] %>% rownames_to_column("Strain") %>% as.data.table() %>% 
  left_join(., m$dr_matches, by = "Strain")

final_steps <- lapply(1:p, function(j) {
  
  outputMessages(paste0("   Collecting ECCs for group of clusters ", j, " / ", p))
  cluster_x <- df[df[[cx]] %in% pull(results[[j]], cx),-"Strain"]
  dms <- readRDS(paste0("results/tmp/dists-", j, ".Rds"))
  
  b1 <- transformData2(dms$temp, "temp", min_temp, max_temp)
  transformed_temp <- b1 %>% formatData(., c("dr1","dr2","Temp.Dist"))

  transformed_geo <- dms$geo %>%
    transformData2(., "geo", min_geo, max_geo) %>%
    formatData(., c("dr1","dr2","Geog.Dist"))
  
  rm(dms); gc()
  
  transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
  
  outputMessages("   Clearing up some memory")
  rm(transformed_temp); rm(transformed_geo); gc()
  
  a1 <- epiCollectionByCluster(m$strain_data, tau, gamma, transformed_dists, k, cluster_x)
  saveRDS(a1, paste0("results/tmp/TP", k, "-", j, "-eccs.Rds"))
})

new_eccs <- lapply(1:p, function(j) {
  readRDS(paste0("results/tmp/TP", k, "-", j, "-eccs.Rds"))
}) %>% bind_rows()

# saveRDS(new_eccs, "results/tmp/new_eccs.Rds")
# # new_eccs[which(new_eccs[,3] == -Inf),3] <- NA
# ---------------------------------------------------------------------------------------------

# Distance matrix step
a1 <- m$assignments %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
identical(a1[rownames(dms$temp), colnames(dms$temp)], dms$temp)

# Transformation step
a2 <- a1 %>% transformData(., "temp")
identical(a2[rownames(b1),colnames(b1)], b1)

tr_temp <- a2 %>% formatData(., c("dr1", "dr2", "Temp.Dist"))

tr_geo <- m$assignments %>% select(dr, Latitude, Longitude) %>%
  distMatrix(., "geo", c("Latitude", "Longitude")) %>%
  pairwiseDists(., "geo", c("dr1", "dr2", "Geog.Dist"))

tr_dists <- merge.data.table(tr_temp, tr_geo)
original_eccs <- tr_dists %>%
  epiCollection(m$strain_data, tau, gamma, typing_data, ., dm_temp, dm_geo, m$dr_matches)

saveRDS(tr_temp, "results/tmp/tr_temp.Rds")
saveRDS(tr_geo, "results/tmp/tr_geo.Rds")
saveRDS(original_eccs, "results/tmp/original_eccs.Rds")

# original_eccs <- readRDS("results/tmp/original_eccs.Rds")[[k]]
colnames(new_eccs)[grepl("ECC", colnames(new_eccs))] %<>% paste0("New_", .)
eccs <- inner_join(new_eccs, original_eccs)
eccs$dif <- abs(pull(eccs, 3) - pull(eccs, 4))
assert("Original ECCs match the new ones (cluster sectioning)", nrow(eccs[!is.na(dif)][dif > 1e-10]) == 0)
  
eccs[!is.na(dif)][dif > 1e-10] %>% as.data.table()



# # a1 <- readRDS("results/collected_eccs.Rds")
# collected_eccs <- lapply(1:length(combos), function(j) {
#   c1 <- unlist(strsplit(combos[j], split = "")) %>% as.numeric()
#   cat(paste0("\n\nStep ", j + 1, ":"))
#   tau <- c1[2]
#   gamma <- c1[3]
#   epiCollectionByCluster(m$strain_data, tau, gamma, transformed_dists, 2, cluster_x)
#   # epiCollection(m$strain_data, tau, gamma, typing_data, 
#   #               transformed_dists, dm_temp, dm_geo, m$dr_matches, j)
# })
# all_equal(a1, original_eccs)
# 
# all_equal(a1[[1]][[2]], collected_eccs[[1]])
# all_equal(a1[[2]][[2]], collected_eccs[[2]])

# cat(paste0("\n\nStep ", length(combos) + 2, ":"))
# outputMessages("   Merging collected ECCs ...\n")
# full_set <- mergeECCs(collected_eccs, 1, tp1$proc) %>%
#   merge.data.table(., mergeECCs(collected_eccs, 2, tp2$proc), by = "Strain", all.y = TRUE) %>%
#   mutate(TP1 = ifelse(is.na(TP1), 0, TP1))
# 
# write.table(full_set, file = "results/ECCs.tsv", col.names = TRUE, 
#             row.names = FALSE, quote = FALSE, sep = "\t")
# 
# stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
# cat(timeTaken(pt = "ECC data collection", stopwatch))
# 
# cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ", 
#            paste0(rep("-", 31), collapse = ""), "||\n"))
