# temporal_dists <- strain_data %>% 
#   select(dr, Date) %>%
#   distMatrix(., "temp", "Date") %>% 
#   transformData(., "temp") %>%
#   formatData(., c("dr1", "dr2", "Temp.Dist"))

# # Distance matrix, then transformation step, then formatting data
# tr_temp <- m$assignments %>% 
#   select(dr, Date) %>%
#   distMatrix(., "temp", "Date") %>% 
#   transformData(., "temp") %>%
#   formatData(., c("dr1", "dr2", "Temp.Dist"))
# 
# tr_geo <- m$assignments %>% select(dr, Latitude, Longitude) %>%
#   distMatrix(., "geo", c("Latitude", "Longitude")) %>%
#   pairwiseDists(., "geo", c("dr1", "dr2", "Geog.Dist"))
# 
# tr_dists <- merge.data.table(tr_temp, tr_geo)
# original_eccs <- tr_dists %>%
#   epiCollection(m$strain_data, tau, gamma, typing_data, ., dm_temp, dm_geo, m$dr_matches)
# # saveRDS(original_eccs, "results/tmp/original_eccs.Rds")
# 
# tp1_eccs <- readRDS("results/tmp/TP1-dists/all_eccs.Rds")
# colnames(tp1_eccs)[grepl("ECC", colnames(tp1_eccs))] %<>% paste0("New_", .)
# 
# inner_join(original_eccs[[1]], tp1_eccs) %>% 
#   mutate(dif = abs(New_TP1_T0_ECC.0.1.0 - TP1_T0_ECC.0.1.0)) %>% 
#   filter(!is.na(dif)) %>% 
#   filter(dif > 0)
# 
# # original_eccs <- readRDS("results/tmp/original_eccs.Rds")[[k]]
# colnames(new_eccs)[grepl("ECC", colnames(new_eccs))] %<>% paste0("New_", .)
# eccs <- inner_join(new_eccs, original_eccs)
# eccs$dif <- abs(pull(eccs, 3) - pull(eccs, 4))
# assert("Original ECCs match the new ones (cluster sectioning)", nrow(eccs[!is.na(dif)][dif > 1e-10]) == 0)
#   
# eccs[!is.na(dif)][dif > 1e-10] %>% as.data.table()
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
