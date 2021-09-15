libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "testthat")
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
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

test_results <- vector(length = 1) %>% 
  setNames(c("cluster_results"))

cgm_results <- read_tsv(arg$cgm_results) %>% as.data.table() %>% 
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
metadata <- read_tsv(arg$metadata) %>% processedStrains()
strain_data <- metadata$strain_data %>% as.data.table()

test_results[[1]] <- test_that("cluster_results", {

  # # Randomly select a cluster to track (TP1)
  # starter_set <- tp2$lookup_table[new_h == '0']
  # random_cluster_row <- sample.int(n = nrow(starter_set), size = 1, replace = FALSE)
  # random_cluster_row <- 402
  # cluster_id <- starter_set[random_cluster_row]
  # 
  # x1 <- tp2$original %>% select(Strain, cluster_id$old_h) %>%
  #   set_colnames(c("Strain", "Hx")) %>%
  #   filter(Hx %in% cluster_id$old_cl)
  # # the strains within the randomly selected cluster (actual cluster name)
  # x1

  # # the strains within the randomly selected cluster (new cluster number)
  a1 <- tp2$new_cols %>% select(Strain, as.character(cluster_id$new_h)) %>% 
    set_colnames(c("Strain", "Hx")) %>% filter(Hx == cluster_id$new_cl)
  
  expect_identical(sort(a1$Strain), sort(x1$Strain))
  
  # cgm cluster results
  cluster_results <- cgm_results[tp1_cluster == cluster_id$new_cl]
  
  # want to be sure when the cluster first formed and when the last novel strain was added
  intervals <- cluster_results %>% arrange(interval) %>% slice(1,n()) %>% pull(interval)
  start_ivl <- substr(intervals[1], 1, 7)
  end_ivl <- substr(intervals[2], 9, 15)

  actual_ivls <- strain_data[Strain %in% a1$Strain] %>% arrange(Date) %>%
    slice(1,n()) %>% pull(Date) %>% substr(., 1, 7)
  
  expect_equal(actual_ivls[1], start_ivl)
  expect_true(actual_ivls[2] <= end_ivl)
  assert("Last strain was introduced to cluster before end date", actual_ivls[2] <= end_ivl)

  rowx <- cluster_results[1,]
  expect_identical(rowx$actual_tp1_cluster, cluster_id$old_cl)
  expect_identical(as.double(rowx$tp1_cluster), as.double(cluster_id$new_cl))

  start_date <- substr(rowx$interval, 1, 7)
  end_date <- substr(rowx$interval, 9, 15)
  end_of_tp2 <- max(strain_data[YearMonth == end_date]$Date)
  end_of_tp1 <- max(strain_data[YearMonth == start_date]$Date)
  tp1_cl <- strain_data[Strain %in% a1$Strain][Date <= end_of_tp1]
  tp2_cl <- strain_data[Strain %in% a1$Strain][Date <= end_of_tp2]
  
  expect_equal(rowx$tp1_cl_size1, nrow(tp1_cl))
  expect_equal(rowx$tp1_cl_size2, nrow(tp1_cl) + 1)
  expect_equal(rowx$tp2_cl_size1, nrow(tp2_cl))
  expect_equal(rowx$tp2_cl_size2, nrow(tp2_cl) + 1)
  
  expect_equal(rowx$actual_cl_growth, nrow(tp2_cl) - nrow(tp1_cl))
  
  # TP1 strains that are not found in tp1_cl
  tp1s_outside_cluster <- strain_data[YearMonth <= start_date][!(Strain %in% tp1_cl$Strain)]
  intersect(tp2_cl$Strain, tp1s_outside_cluster$Strain) %>% length() %>%
    expect_equal(rowx$add_tp1, .)

  # strains found at TP2 but not at TP1
  present_at_tp2 <- strain_data[Date <= end_of_tp2]
  originals <- strain_data[Date <= end_of_tp1]
  number_of_novels <- length(intersect(a1$Strain, novels$Strain)) - length(intersect(a1$Strain, originals$Strain))
  expect_equal(rowx$num_novs, number_of_novels)
  
  tp1gdists <- tp1_cl %>% select(Longitude, Latitude) %>% 
    as.data.frame() %>% earth.dist(dist = TRUE) %>% as.matrix()
  tp1_geo_avg_dist <- sum(tp1gdists[upper.tri(tp1gdists)]) / length(tp1gdists[upper.tri(tp1gdists)])
  expect_equal(rowx$tp1_geo.avg.dist, tp1_geo_avg_dist)
  
  tp1tdists <- tp1_cl %>% pull(Date) %>% 
    dist(diag = FALSE, upper = FALSE, method = "euclidean") %>% as.matrix()
  tp1_temp_avg_dist <- sum(tp1tdists[upper.tri(tp1tdists)]) / length(tp1tdists[upper.tri(tp1tdists)])
  expect_equal(rowx$tp1_temp.avg.dist, tp1_temp_avg_dist)
  
  tp2gdists <- tp2_cl %>% select(Longitude, Latitude) %>% 
    as.data.frame() %>% earth.dist(dist = TRUE) %>% as.matrix()
  tp2_geo_avg_dist <- sum(tp2gdists[upper.tri(tp2gdists)]) / length(tp2gdists[upper.tri(tp2gdists)])
  expect_equal(rowx$tp2_geo.avg.dist, tp2_geo_avg_dist)

  # FAILS 
  # tp2tdists <- tp2_cl %>% pull(Date) %>% 
  #   dist(diag = FALSE, upper = FALSE, method = "euclidean") %>% as.matrix()
  # tp2_temp_avg_dist <- sum(tp2tdists[upper.tri(tp2tdists)]) / length(tp2tdists[upper.tri(tp2tdists)])
  # expect_equal(rowx$tp2_temp.avg.dist, tp2_temp_avg_dist)
  
  expect_equal(rowx$average_tp1_latitude, mean(tp1_cl$Latitude))
  expect_equal(rowx$average_tp2_latitude, mean(tp2_cl$Latitude))
  
  expect_equal(rowx$average_tp1_longitude, mean(tp1_cl$Longitude))
  expect_equal(rowx$average_tp2_longitude, mean(tp2_cl$Longitude))
  
  expect_identical(as.character(rowx$average_tp1_date), as.character(mean(tp1_cl$Date)))
  expect_identical(as.character(rowx$average_tp2_date), as.character(mean(tp2_cl$Date)))
  
  # first (by threshold) tp2 cluster to contain all of the tp1 strains
  # note that these are the same if we are looking at the same cluster sets
  # but if we someday want to revisit separately clustered time points again, 
  # this will become highly relevant
  tp2_cluster <- tp2$new_cols[Strain %in% tp2_cl$Strain][Strain %in% tp1_cl$Strain] %>% 
    melt.data.table(id.vars = "Strain") %>% group_by(variable, value) %>% 
    summarise(n = n(), .groups = "drop") %>% filter(n >= nrow(tp1_cl)) %>% slice(1)
  expect_equal(rowx$tp2_cluster, tp2_cluster$value)
  
  starter_set[new_h == as.double(as.character(tp2_cluster$variable)) & 
                new_cl == tp2_cluster$value]$old_cl %>% 
    expect_equal(rowx$actual_tp2_cluster, .)
  
  tpx1 <- tp2$new_cols[Strain %in% strain_data[Date <= end_of_tp1]$Strain]
  tpx2 <- tp2$new_cols[Strain %in% strain_data[Date <= end_of_tp2]$Strain]
  ph <- max(nchar(colnames(tpx1)[-1]), nchar(colnames(tpx2)[-1]))
  pc <- tpx2 %>% select(-Strain) %>% max(., tpx2 %>% select(-Strain)) %>% nchar()
  
  c1 <- "TP1"
  c2 <- tp2_cluster %>% pull(variable) %>% as.character() %>% as.integer() %>% 
    formatC(., width = max(3, ph), format = "d", flag = "0") %>% paste0("h", .)
  c3 <- tp2_cluster %>% pull(value) %>% as.character() %>% as.integer() %>% 
    formatC(., width = max(3, pc), format = "d", flag = "0") %>% paste0("c", .)
  expect_identical(rowx$first_tp1_flag, paste(c1, c2, c3, sep = "_"))
})



"interval"
"actual_tp1_cluster"
"tp1_cluster"
"tp1_cl_size1"
# "tp1_ecc.0.1.0"
# "tp1_ecc.0.0.1"
"actual_tp2_cluster"
"tp2_cluster"
"tp2_cl_size1"
# "tp2_ecc.0.1.0"
# "tp2_ecc.0.0.1"
# "delta_ecc_0.1.0"
# "delta_ecc_0.0.1"
"average_tp1_date"
"tp1_temp.avg.dist"
"average_tp1_latitude"
"average_tp1_longitude"
"tp1_geo.avg.dist"
"average_tp2_date"
# "tp2_temp.avg.dist"
"average_tp2_latitude"
"average_tp2_longitude"
"tp2_geo.avg.dist"
# "first_tp1_flag"
# "last_tp1_flag"
# "first_tp2_flag"
# "last_tp2_flag"
"tp1_cl_size2"
"tp2_cl_size2"
# "actual_cl_growth"
"add_tp1"
"num_novs"
# "actual_growth_rate"
# "novel_growth"
# "type"




# (4) outputMessages() function
test_results[[3]] <- test_that("outputMessages()", {
  expect_identical(outputMessages(""), cat(paste0("\n","")))
  expect_output(outputMessages("Test"), "\nTest")
  expect_silent(outputMessages())
})


test_results[[4]] <- test_that("transformData2() - temp", {
  fake_obj <- fakeSecCluster(sample.int(50,1), sample.int(200,1)*2, sample.int(20,1))
  results <- fake_obj$b$results
  inds <- fake_obj$b$drs[["th"]] %in% pull(results[[1]], "th")
  cluster_x <- fake_obj$b$drs[inds,-"Strain"]
  cluster_asmts <- fake_obj$a[dr %in% pull(cluster_x, dr)]
  
  dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
  
  returned <- transformData2(dm_temp, "temp", min(dm_temp), max(dm_temp))
  
  logtr <- dm_temp %>% add(10) %>% log10()
  logtr[logtr == -Inf] <- 0
  min_dm <- min(dm_temp) %>% add(10) %>% log10()
  max_dm <- max(dm_temp) %>% add(10) %>% log10()
  
  pairs <- expand.grid(cluster_asmts$dr, cluster_asmts$dr)
  
  for (i in 1:nrow(pairs)) {
    ri <- as.character(pairs$Var1[i])
    ci <- as.character(pairs$Var2[i])
    
    x1 <- (logtr[ri, ci] - min_dm) / (max_dm - min_dm)
    expect_equal(returned[ri, ci], x1)  
  }
})




test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input ECC tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}

# "interval"
# "actual_tp1_cluster"
# "tp1_cluster"
# "tp1_cl_size1"
# "tp1_ecc.0.1.0"
# "tp1_ecc.0.0.1"
# "actual_tp2_cluster"
# "tp2_cluster"
# "tp2_cl_size1"
# "tp2_ecc.0.1.0"
# "tp2_ecc.0.0.1"
# "delta_ecc_0.1.0"
# "delta_ecc_0.0.1"
# "average_tp1_date"
# "tp1_temp.avg.dist"
# "average_tp1_latitude"
# "average_tp1_longitude"
# "tp1_geo.avg.dist"
# "average_tp2_date"
# "tp2_temp.avg.dist"
# "average_tp2_latitude"
# "average_tp2_longitude"
# "tp2_geo.avg.dist"
# "first_tp1_flag"
# "last_tp1_flag"
# "first_tp2_flag"
# "last_tp2_flag"
# "tp1_cl_size2"
# "tp2_cl_size2"
# "actual_cl_growth"
# "add_tp1"
# "num_novs"
# "actual_growth_rate"
# "novel_growth"
# "type"
