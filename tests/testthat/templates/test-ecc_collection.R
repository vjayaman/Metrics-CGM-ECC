libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/ecc_opener.R")

test_results <- vector(length = 12) %>% 
  setNames(c("Part 1", "Part 2", "Part 3", "Part 4", "Part 5", "Part 6"))

# Part 1
# assert("Distances were collected and saved", file.exists("results/TP1/dists/group1.Rds"))
# extremes <- readRDS("results/dist_extremes.Rds")
# dfx <- expand.grid(x = combos, k = c(1,2)) %>% as.data.frame()
# for (k in unique(dfx$k)) {
#   dir.create(paste0("results/TP", k, "/eccs"), showWarnings = FALSE)
# }
test_that("Part 1", {})

# Part 2
# for (i in 1:nrow(dfx)) {
#   c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "") %>% 
#     unlist() %>% as.numeric()
#   paste0("results/TP", dfx$k[i], "/eccs/", as.character(dfx$x[i])) %>% 
#     dir.create(., showWarnings = FALSE)
#   outputMessages(paste0("\nCollecting and saving ECCs for groups of clusters at TP", 
#                         dfx$k[i], ", for ECC coefficient triple ", as.character(dfx$x[i])))
#   sectionClusters(dfx$k[i], typing_data, m) %>% 
#     collectECCs(dfx$k[i], m, ., extremes, c1, 
#                 paste0("results/TP", dfx$k[i], "/dists/"), 
#                 paste0("results/TP", dfx$k[i], "/eccs/", 
#                        as.character(dfx$x[i]), "/"))
#   outputMessages(paste0("\nMerging ECC files at TP", dfx$k[i]))
# }
test_that("Part 2", {})

# Part 3
# for (i in 1:nrow(dfx)) {
#   dir_k <- paste0("results/TP", dfx$k[i], "/eccs/", dfx$x[i], "/")
#   fnames <- list.files(dir_k)
#   tpk <- lapply(fnames, function(f) {
#     ecc_j <- readRDS(paste0(dir_k, f))
#     file.remove(paste0(dir_k, f))
#     return(ecc_j)
#   }) %>% bind_rows()
#   saveRDS(tpk, paste0(dir_k, "TP", dfx$k[i], "-", dfx$x[i], "-eccs.Rds"))
# }
test_that("Part 3", {})

# Part 4
# a1 <- lapply(combos, function(x) {
#   readRDS(paste0("results/TP1/eccs/", x, "/TP1-", x, "-eccs.Rds"))
# }) %>% Reduce(inner_join, .) %>% unique()
# b1 <- lapply(combos, function(x) {
#   readRDS(paste0("results/TP2/eccs/", x, "/TP2-", x, "-eccs.Rds"))
# }) %>% Reduce(inner_join, .) %>% unique()
# tp1_eccs <- tp1$proc[,c("Strain",colnames(a1)[1])] %>% left_join(., a1)
# tp2_eccs <- tp2$proc[,c("Strain",colnames(b1)[1])] %>% left_join(., b1)
# all_eccs <- right_join(tp1_eccs, tp2_eccs, by = "Strain")

# tp1_eccs <- tp1$proc[,c("Strain", colnames(a1[[1]])[1])] %>% 
#   left_join(., a1[[1]], by = intersect(colnames(.), colnames(a1[[1]]))) %>% 
#   left_join(., a1[[2]], by = intersect(colnames(.), colnames(a1[[2]])))
test_that("Part 4", {
  
})

# Part 5
# tp1size <- colnames(all_eccs) %>% 
#   grep("Size", ., value = TRUE) %>% grep("TP1", ., value = TRUE)
# inf_inds <- which(all_eccs == -Inf, arr.ind = TRUE) %>% as.data.frame() %>% pull("row")
# num_others <- all_eccs[inf_inds,] %>% filter(!is.na(!!as.symbol(tp1size))) %>% 
#   pull(tp1size) %>% unique() %>% setdiff(., 1) %>% length()
# assert("The -Infs are because of TP1 singletons", num_others == 0)
# all_eccs[all_eccs == -Inf] <- NA
test_that("Part 5", {})

# Part 6
# assert("Average distances were collected and saved", file.exists("results/average_dists.Rds"))
# all_avg_dists <- readRDS("results/average_dists.Rds")
# cat(paste0("\n\nStep ", length(combos) + 2, ":"))
# outputMessages("   Merging collected ECCs ...\n")
# inner_join(all_avg_dists, all_eccs) %>% 
#   write.table(., file = "results/ECCs.tsv", col.names = TRUE, 
#               row.names = FALSE, quote = FALSE, sep = "\t")
test_that("Part 6", {})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input ECC collection tests have passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}
