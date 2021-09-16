libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "testthat")
y <- suppressWarnings(suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE)))
assert("All packages loaded correctly", all(unlist(y)))

source("tests/testthat/merged_results_functions.R")

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-b", "--tp2"), metavar = "file", 
              default = "inputs/processed/allTP2.Rds", help = "TP2 data"), 
  make_option(c("-c", "--cgm_results"), metavar = "file", 
              default = "results/Merged_cluster_results.tsv", help = "CGM results"))

arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)

test_results <- vector(length = 27) %>% 
  setNames(c("strains", "intervals", "actual_tp1_cluster and tp1_cluster", "cluster sizes", 
             "cluster growth", "sneakers", "novels", "average tp1 geo distance", 
             "average tp1 temporal distance", "average tp2 geo distance", "average tp2 temporal distance", 
             "average tp1 latitude (raw)", "average tp2 latitude (raw)", "average tp1 longitude (raw)", 
             "average tp4 longitude (raw)", "average tp1 date", "average tp2 date", "first tp1 flag", 
             "last tp1 flag", "first tp2 flag", "last tp2 flag", "actual tp2 cluster", 
             "actual growth rate", "number of novels", "type", "eccs", "delta eccs"))

cgm_results <- suppressMessages(read_tsv(arg$cgm_results)) %>% as.data.table() %>% 
  set_colnames(tolower(colnames(.))) %>% 
  rename("tp1_cl_size1" = "tp1 cluster size (1)", 
         "tp2_cl_size1" = "tp2 cluster size (1)", 
         "first_tp1_flag" = "first time this cluster was seen in tp1", 
         "last_tp1_flag" = "last time this cluster was seen in tp1", 
         "first_tp2_flag" = "first time this cluster was seen in tp2", 
         "last_tp2_flag" = "last time this cluster was seen in tp2", 
         "tp1_cl_size2" = "tp1 cluster size + 1 (2)", 
         "tp2_cl_size2" = "tp2 cluster size + 1 (2)", 
         "actual_cl_growth" = "actual cluster growth (tp2 size - tp1 size)", 
         "add_tp1" = "number of additional tp1 strains in the tp2 match", 
         "num_novs" = "number of novels in the tp2 match", 
         "actual_growth_rate" = "actual growth rate = (tp2 size - tp1 size) / (tp1 size)", 
         "novel_growth" = "novel growth = (tp2 size) / (tp2 size - number of novels)")
colnames(cgm_results) <- colnames(cgm_results) %>% gsub(" ", "_", .)

tp2 <- readRDS(arg$tp2)
tp2clusters <- tp2$new_cols %>% melt.data.table(id.vars = "Strain")

metadata <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()
strain_data <- metadata$strain_data %>% as.data.table()


# Randomly select a cluster to track (TP1)
starter_set <- tp2$lookup_table[new_h == '0']
random_cluster_row <- sample.int(n = nrow(starter_set), size = 1, replace = FALSE)
# random_cluster_row <- 134
# print(random_cluster_row)
cluster_id <- starter_set[random_cluster_row]

# # the strains within the randomly selected cluster (new cluster number)
a1 <- tp2$new_cols %>% select(Strain, as.character(cluster_id$new_h)) %>% 
  set_colnames(c("Strain", "Hx")) %>% filter(Hx == cluster_id$new_cl)

# cgm cluster results
cluster_results <- cgm_results[tp1_cluster == cluster_id$new_cl]
rowx <- cluster_results[1,]

start_date <- substr(rowx$interval, 1, 7)
end_date <- substr(rowx$interval, 9, 15)
end_of_tp2 <- max(strain_data[YearMonth == end_date]$Date)
end_of_tp1 <- max(strain_data[YearMonth == start_date]$Date)
tp1_cl <- strain_data[Strain %in% a1$Strain][Date <= end_of_tp1]
tp2_cl <- strain_data[Strain %in% a1$Strain][Date <= end_of_tp2]

tpx1 <- tp2$new_cols[Strain %in% strain_data[Date <= end_of_tp1]$Strain]
tpx2 <- tp2$new_cols[Strain %in% strain_data[Date <= end_of_tp2]$Strain]
ph <- max(nchar(colnames(tpx1)[-1]), nchar(colnames(tpx2)[-1]))
pc <- tpx2 %>% select(-Strain) %>% max(., tpx2 %>% select(-Strain)) %>% nchar()



test_results[[1]] <- test_that("strains", {
  x1 <- tp2$original %>% select(Strain, cluster_id$old_h) %>%
    set_colnames(c("Strain", "Hx")) %>%
    filter(Hx %in% cluster_id$old_cl)
  # the strains within the randomly selected cluster (actual cluster name): x1

  expect_identical(sort(a1$Strain), sort(x1$Strain))
})

test_results[[2]] <- test_that("intervals", {
  # want to be sure when the cluster first formed and when the last novel strain was added
  intervals <- cluster_results %>% arrange(interval) %>% slice(1,n()) %>% pull(interval)
  start_ivl <- substr(intervals[1], 1, 7)
  end_ivl <- substr(intervals[2], 9, 15)
  
  actual_ivls <- strain_data[Strain %in% a1$Strain] %>% arrange(Date) %>%
    slice(1,n()) %>% pull(Date) %>% substr(., 1, 7)
  
  assert("Last strain was introduced to cluster before end date", actual_ivls[2] <= end_ivl)
  expect_equal(actual_ivls[1], start_ivl)
  expect_true(actual_ivls[2] <= end_ivl)
})


test_results[[3]] <- test_that("actual_tp1_cluster and tp1_cluster", {
  expect_identical(rowx$actual_tp1_cluster, cluster_id$old_cl)
  expect_identical(as.double(rowx$tp1_cluster), as.double(cluster_id$new_cl))
})


test_results[[4]] <- test_that("cluster sizes", {
  expect_equal(rowx$tp1_cl_size1, nrow(tp1_cl))
  expect_equal(rowx$tp1_cl_size2, nrow(tp1_cl) + 1)
  expect_equal(rowx$tp2_cl_size1, nrow(tp2_cl))
  expect_equal(rowx$tp2_cl_size2, nrow(tp2_cl) + 1)
})

test_results[[5]] <- test_that("cluster growth", {
  expect_equal(rowx$actual_cl_growth, nrow(tp2_cl) - nrow(tp1_cl))
})

test_results[[6]] <- test_that("sneakers", {
  # TP1 strains that are not found in tp1_cl
  tp1s_outside_cluster <- strain_data[YearMonth <= start_date][!(Strain %in% tp1_cl$Strain)]
  intersect(tp2_cl$Strain, tp1s_outside_cluster$Strain) %>% length() %>%
    expect_equal(rowx$add_tp1, .)
})

present_at_tp2 <- strain_data[Date <= end_of_tp2]
originals <- strain_data[Date <= end_of_tp1]
number_of_novels <- length(intersect(a1$Strain, present_at_tp2$Strain)) - length(intersect(a1$Strain, originals$Strain))

test_results[[7]] <- test_that("novels", {
  # strains found at TP2 but not at TP1
  expect_equal(rowx$num_novs, number_of_novels)
})

test_results[[8]] <- test_that("average tp1 geo distance", {
  tp1gdists <- tp1_cl %>% select(Longitude, Latitude) %>% 
    as.data.frame() %>% earth.dist(dist = TRUE) %>% as.matrix()
  denom_pri <- length(tp1gdists[upper.tri(tp1gdists)])
  denom <- ifelse(denom_pri == 0, 1, denom_pri)
  tp1_geo_avg_dist <- sum(tp1gdists[upper.tri(tp1gdists)]) / denom
  expect_equal(rowx$tp1_geo.avg.dist, tp1_geo_avg_dist)
})

test_results[[9]] <- test_that("average tp1 temporal distance", {
  tp1tdists <- tp1_cl %>% pull(Date) %>% 
    dist(diag = FALSE, upper = FALSE, method = "euclidean") %>% as.matrix()
  tp1tdists <- tp1tdists[upper.tri(tp1tdists)]
  denom_pri <- length(tp1tdists)
  denom <- ifelse(denom_pri == 0, 1, denom_pri)
  tp1_temp_avg_dist <- sum(tp1tdists) / denom
  expect_equal(rowx$tp1_temp.avg.dist, tp1_temp_avg_dist)
})

test_results[[10]] <- test_that("average tp2 geo distance", {
  tp2gdists <- tp2_cl %>% select(Longitude, Latitude) %>% 
    as.data.frame() %>% earth.dist(dist = TRUE) %>% as.matrix()
  denom_pri <- length(tp2gdists[upper.tri(tp2gdists)])
  denom <- ifelse(denom_pri == 0, 1, denom_pri)
  tp2_geo_avg_dist <- sum(tp2gdists[upper.tri(tp2gdists)]) / denom
  expect_equal(rowx$tp2_geo.avg.dist, tp2_geo_avg_dist)
})

test_results[[11]] <- test_that("average tp1 temporal distance", {
  tp2tdists <- tp2_cl %>% pull(Date) %>% 
    dist(diag = FALSE, upper = FALSE, method = "euclidean") %>% as.matrix()
  denom_pri <- length(tp2tdists[upper.tri(tp2tdists)])
  denom <- ifelse(denom_pri == 0, 1, denom_pri)
  tp2_temp_avg_dist <- sum(tp2tdists[upper.tri(tp2tdists)]) / denom
  expect_equal(rowx$tp2_temp.avg.dist, tp2_temp_avg_dist)
})

test_results[[12]] <- test_that("average tp1 latitude (raw)", {
  expect_equal(rowx$average_tp1_latitude, mean(tp1_cl$Latitude))
})

test_results[[13]] <- test_that("average tp2 latitude (raw)", {
  expect_equal(rowx$average_tp2_latitude, mean(tp2_cl$Latitude))
})

test_results[[14]] <- test_that("average tp1 longitude (raw)", {
  expect_equal(rowx$average_tp1_longitude, mean(tp1_cl$Longitude))
})

test_results[[15]] <- test_that("average tp4 longitude (raw)", {
  expect_equal(rowx$average_tp2_longitude, mean(tp2_cl$Longitude))
})

test_results[[16]] <- test_that("average tp1 date", {
  expect_identical(as.character(rowx$average_tp1_date), as.character(mean(tp1_cl$Date)))
})

test_results[[17]] <- test_that("average tp2 date", {
  expect_identical(as.character(rowx$average_tp2_date), as.character(mean(tp2_cl$Date)))
})

test_results[[18]] <- test_that("first tp1 flag", {
  # first (by threshold) tp2 cluster to contain all of the tp1 strains
  # note that these are the same if we are looking at the same cluster sets
  # but if we someday want to revisit separately clustered time points again, 
  # this will become highly relevant
  
  # only strains available at TP1
  alltp1data <- tp2clusters[Strain %in% strain_data[Date <= end_of_tp1]$Strain]
  tp1sizes <- alltp1data %>% group_by(variable, value) %>% summarise(n = n(), .groups = "drop") %>% as.data.table()
  
  # of those, which were the first and last clusters to contain only the tp1_cl strains?
  tracked_tp1s <- alltp1data[Strain %in% tp1_cl$Strain] %>% select(-Strain) %>% unique() %>% 
    merge.data.table(., tp1sizes)
  
  first_tp1 <- tracked_tp1s %>% slice(1)
  first_tp1 %>% newID(., "tp1", "variable", "value", ph, pc) %>% pull(id) %>% 
    expect_identical(rowx$first_tp1_flag, .)
})


test_results[[19]] <- test_that("last tp1 flag", {
  # only strains available at TP1
  alltp1data <- tp2clusters[Strain %in% strain_data[Date <= end_of_tp1]$Strain]
  tp1sizes <- alltp1data %>% group_by(variable, value) %>% summarise(n = n(), .groups = "drop") %>% as.data.table()
  
  # of those, which were the first and last clusters to contain only the tp1_cl strains?
  tracked_tp1s <- alltp1data[Strain %in% tp1_cl$Strain] %>% select(-Strain) %>% unique() %>% 
    merge.data.table(., tp1sizes)
  
  first_tp1 <- tracked_tp1s %>% slice(1)
  
  tracked_tp1s[n == first_tp1$n] %>% slice(n()) %>% 
    newID(., "tp1", "variable", "value", ph, pc) %>% pull(id) %>% 
    expect_identical(rowx$last_tp1_flag, .)
})

# only strains available at TP2
alltp2data <- tp2clusters[Strain %in% strain_data[Date <= end_of_tp2]$Strain]
tp2sizes <- alltp2data %>% group_by(variable, value) %>% summarise(n = n(), .groups = "drop") %>% as.data.table()

tracked_tp2s <- alltp2data[Strain %in% tp2_cl$Strain] %>% select(-Strain) %>% unique() %>% 
  merge.data.table(., tp2sizes)

# first TP2 cluster to contain the tp1_cl strains
first_tp2 <- tracked_tp2s %>% slice(1)

test_results[[20]] <- test_that("first tp2 flag", {
  first_tp2 %>% newID(., "tp2", "variable", "value", ph, pc) %>% pull(id) %>% 
    expect_identical(rowx$first_tp2_flag, .)
})


test_results[[21]] <- test_that("last tp2 flag", {
  tracked_tp2s[n == first_tp2$n] %>% slice(n()) %>% 
    newID(., "tp2", "variable", "value", ph, pc) %>% pull(id) %>% 
    expect_identical(rowx$last_tp2_flag, .)
  expect_equal(rowx$tp2_cluster, first_tp2$value)
})

test_results[[22]] <- test_that("actual tp2 cluster", {
  starter_set[new_h == as.double(as.character(first_tp2$variable)) & 
                new_cl == first_tp2$value]$old_cl %>% 
    expect_equal(rowx$actual_tp2_cluster, .)
})

test_results[[23]] <- test_that("actual growth rate", {
  ((nrow(tp2_cl) - nrow(tp1_cl)) / (nrow(tp1_cl)+1)) %>% round(., digits = 3) %>% 
    expect_equal(rowx$actual_growth_rate, .)
})

test_results[[24]] <- test_that("number of novels", {
  ((nrow(tp2_cl) + 1) / (nrow(tp2_cl) + 1 - number_of_novels)) %>% round(., digits = 3) %>% 
    expect_equal(rowx$novel_growth, number_of_novels)
})

test_results[[25]] <- test_that("type", {
  if (nrow(tp1_cl) > 1 & nrow(tp2_cl) > 1 & nrow(tp1_cl) == nrow(tp2_cl)) {
    expect_equal(rowx$type, "Type1")
  }else if (nrow(tp1_cl) > 1 & nrow(tp2_cl) > 1 & nrow(tp2_cl) > nrow(tp1_cl)) {
    expect_equal(rowx$type, "Type2")
  }else if (nrow(tp1_cl) < 2 & nrow(tp2_cl) > 1) {
    expect_equal(rowx$type, "Type3")
  }else {
    expect_equal(rowx$type, "Type4")
  }
})

ecc_cols <- grep("ecc", colnames(cgm_results), value = TRUE) %>% 
  grep("delta", ., value = TRUE, invert = TRUE)

test_results[[26]] <- test_that("eccs", {
  extremes <- readRDS("intermediate_data/TPN/extreme_dists.Rds")
  tpcls <- list(tp1 = tp1_cl, tp2 = tp2_cl)
  tr_dists <- list(tp1 = manualEpiMelt(tp1_cl, extremes), 
                   tp2 = manualEpiMelt(tp2_cl, extremes))
  
  for (tp_ecc in ecc_cols) {
    tpx <- strsplit(tp_ecc, "_") %>% unlist() %>% extract2(1)
    tau <- strsplit(tp_ecc, "\\.") %>% unlist() %>% extract2(3) %>% as.double()
    gamma <- strsplit(tp_ecc, "\\.") %>% unlist() %>% extract2(4) %>% as.double()
    
    epi_melt <- tr_dists[[tpx]] %>% 
      mutate(Total.Distinv = 1 - sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
      select(Strain.1, Strain.2, Total.Distinv) %>% 
      set_colnames(c("Var1", "Var2", "value"))
    
    cluster_size <- nrow(tpcls[[tpx]])
    ecc_result <- pull(rowx, paste(paste0(tpx, "_ecc.0"), tau, gamma, sep = "."))
    if (cluster_size == 1) {
      if (is.na(ecc_result)) {
        expect_true(TRUE)
      }else {
        expect_equal(ecc_result, 1)
      }
    }else {
      ((sum(epi_melt$value) - cluster_size) / (cluster_size * (cluster_size - 1))) %>% 
        expect_equal(ecc_result, .)  
    }
  }
})

test_results[[27]] <- test_that("delta eccs", {
  params <- gsub("tp1_ecc.|tp2_ecc.", "", ecc_cols) %>% unique()
  coeffs <- data.table(tau = as.double(substr(params, 3, 3)), 
                       gamma = as.double(substr(params, 5, 5)))
  
  for (i in 1:nrow(coeffs)) {
    cnames <- list(tp1 = paste0("tp1_ecc.0.", coeffs$tau[i], ".", coeffs$gamma[i]), 
                   tp2 = paste0("tp2_ecc.0.", coeffs$tau[i], ".", coeffs$gamma[i]), 
                   delta = paste0("delta_ecc_0.", coeffs$tau[i], ".", coeffs$gamma[i]))
    expect_equal(pull(rowx, cnames$tp2) - pull(rowx, cnames$tp1), pull(rowx, cnames$delta))
  }
})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input cluster result tests pass, for randomly sampled cluster:\n")
  print(cluster_id)
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}


