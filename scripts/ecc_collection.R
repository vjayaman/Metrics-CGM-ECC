#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

# msg <- file("logs/logfile_epiquant.txt", open="wt")
# sink(msg, type="message")

source("scripts/tmp.R")

k <- 2
df <- typing_data[[k]] %>% rownames_to_column("Strain") %>%
  as_tibble() %>% left_join(., dr_matches, by = "Strain")

gc()

results <- sectionTypingData(df, 8000)
assert("No clusters overlooked", length(setdiff(pull(df,2), pull(rbindlist(results),1))) == 0)

cx <- setdiff(colnames(df), c("Strain", "dr"))

c1 <- unlist(strsplit(combos[1], split = "")) %>% as.numeric()
tau <- c1[2]
gamma <- c1[3]

# p <- 5 # length(results)
# dir.create("results/tmp", showWarnings = FALSE)
# new_eccs <- lapply(1:p, function(j) {
  j <- 1
  x1 <- df %>% filter(!!as.symbol(cx) %in% pull(results[[j]], cx))
  cluster_asmts <- assignments %>% filter(dr %in% pull(x1, dr))
  cluster_x <- x1 %>% select(-Strain)
  
  outputMessages("   Generating all possible date pair distances ...")
  transformed_temp <- cluster_asmts %>% select(dr, Date) %>% 
    distMatrix(., "temp", "Date") %>% 
    pairwiseDists(., "temp", c("dr1", "dr2", "Temp.Dist"))

  outputMessages("   Generating all possible lat-long pair distances ...")
  transformed_geo <- cluster_asmts %>% select(dr, Latitude, Longitude) %>% 
    distMatrix(., "geo", c("Latitude", "Longitude")) %>% 
    pairwiseDists(., "geo", c("dr1", "dr2", "Geog.Dist"))

  gc()
  
  transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
  
  outputMessages("   Clearing up some memory")
  rm(transformed_temp)
  rm(transformed_geo)
  gc()
  
  a1 <- epiCollectionByCluster(strain_data, tau, gamma, transformed_dists, k, cluster_x)
  saveRDS(a1, paste0("results/tmp/all-TP", k, "-", j, ".Rds"))
  
  rm(cluster_x)
  rm(transformed_dists)
  # rm(a1)
  gc()
  
  
# }) %>% bind_rows()

# new_eccs <- lapply(1:length(results), function(j) {
#   readRDS(paste0("results/tmp/TP", k, "-", j, ".Rds"))
# }) %>% bind_rows()
# new_eccs[which(new_eccs[,3] == -Inf),3] <- NA
  

tr_temp <- assignments %>% select(dr, Date) %>%
  distMatrix(., "temp", "Date") %>%
  pairwiseDists(., "temp", c("dr1", "dr2", "Temp.Dist"))
tr_geo <- assignments %>% select(dr, Latitude, Longitude) %>%
  distMatrix(., "geo", c("Latitude", "Longitude")) %>%
  pairwiseDists(., "geo", c("dr1", "dr2", "Geog.Dist"))
tr_dists <- merge.data.table(tr_temp, tr_geo)
# original_eccs <- tr_dists %>% 
#   epiCollection(strain_data, tau, gamma, typing_data, ., dm_temp, dm_geo, dr_matches)

original_eccs <- readRDS("results/tmp/original_eccs.Rds")[[k]]
all_equal(a1, original_eccs)
  
  

# # a1 <- readRDS("results/collected_eccs.Rds")
# collected_eccs <- lapply(1:length(combos), function(j) {
#   c1 <- unlist(strsplit(combos[j], split = "")) %>% as.numeric()
#   cat(paste0("\n\nStep ", j + 1, ":"))
#   tau <- c1[2]
#   gamma <- c1[3]
#   epiCollectionByCluster(strain_data, tau, gamma, transformed_dists, 2, cluster_x)
#   # epiCollection(strain_data, tau, gamma, typing_data, 
#   #               transformed_dists, dm_temp, dm_geo, dr_matches, j)
# })

# all_equal(a1[[1]][[2]], collected_eccs[[1]])
# all_equal(a1[[2]][[2]], collected_eccs[[2]])

cat(paste0("\n\nStep ", length(combos) + 2, ":"))
outputMessages("   Merging collected ECCs ...\n")
full_set <- mergeECCs(collected_eccs, 1, tp1$proc) %>%
  merge.data.table(., mergeECCs(collected_eccs, 2, tp2$proc), by = "Strain", all.y = TRUE) %>%
  mutate(TP1 = ifelse(is.na(TP1), 0, TP1))

write.table(full_set, file = "results/ECCs.tsv", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))
