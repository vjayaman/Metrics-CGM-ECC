libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "testthat")
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/ECC/functions/ecc_functions.R")
source("scripts/ECC/functions/dist_functions.R")
source("tests/testthat/helper_functions.R")

test_results <- vector(length = 9) %>% 
  setNames(c("timeTaken() - null inputs", "timeTaken() - time inputs", 
             "outputMessages()", "transformData2() - temp", "transformData2() - geo", 
             "processedStrains() - requirements test", "processedStrains() - additional columns test", 
             "epiCohesion()", "epiCollectionByCluster()"))


# (3) timeTaken() function
# - Indicates length of a process in hours, minutes, and seconds, when given a name of the process 
# ("pt") and a two-element named vector with Sys.time() values named "start_time" and "end_time"
test_results[[1]] <- test_that("timeTaken() - null inputs", {
  stopwatch <- list("start_time" = NULL, "end_time" = NULL)
  expect_identical(timeTaken("test", stopwatch), "Neither start nor end time were collected")
  
  stopwatch$end_time <- as.character.POSIXt(Sys.time())
  expect_identical(timeTaken("test", stopwatch), "Start time was not collected.")
  
  stopwatch$start_time <- as.character.POSIXt(Sys.time())
  stopwatch$end_time <- NULL
  expect_identical(timeTaken("test", stopwatch), "End time was not collected.")
})

test_results[[2]] <- test_that("timeTaken() - time inputs", {
  x <- as.character.POSIXt(Sys.time())
  stopwatch <- list("start_time" = x, "end_time" = NULL)
  stopwatch$end_time <- as.character.POSIXt(as.POSIXct(x) + 30)
  expect_identical(timeTaken("test", stopwatch), "\nThe test process took 30 second(s).")
  
  stopwatch$end_time <- as.character.POSIXt(as.POSIXct(x) + 2*60 + 30)
  expect_identical(timeTaken("test", stopwatch), "\nThe test process took 2 minute(s) and 30 second(s).")
  
  stopwatch$end_time <- as.character.POSIXt(as.POSIXct(x) + 60*60 + 2*60 + 30)
  expect_identical(timeTaken("test", stopwatch), 
                   "\nThe test process took 1 hour(s), 2 minute(s), and 30 second(s).")
})


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

test_results[[5]] <- test_that("transformData2() - geo", {
  fake_obj <- fakeSecCluster(sample.int(50,1), sample.int(200,1)*2, sample.int(20,1))
  results <- fake_obj$b$results
  inds <- fake_obj$b$drs[["th"]] %in% pull(results[[1]], "th")
  cluster_x <- fake_obj$b$drs[inds,-"Strain"]
  cluster_asmts <- fake_obj$a[dr %in% pull(cluster_x, dr)]
  
  dm_geo <- cluster_asmts %>% select(dr, Longitude, Latitude) %>% 
    distMatrix(., "geo", c("Longitude", "Latitude"))
  
  returned <- transformData2(dm_geo, "geo", min(dm_geo), max(dm_geo))
  
  logtr <- dm_geo %>% add(10) %>% log10()
  min_dm <- min(dm_geo) %>% add(10) %>% log10()
  max_dm <- max(dm_geo) %>% add(10) %>% log10()
  
  pairs <- expand.grid(cluster_asmts$dr, cluster_asmts$dr)
  
  if (max_dm == 1) {
    expect_true(all(unique(returned) == 0))
  }else {
    for (i in 1:nrow(pairs)) {
      ri <- as.character(pairs$Var1[i])
      ci <- as.character(pairs$Var2[i])
      
      x1 <- (logtr[ri, ci] - min_dm) / (max_dm - min_dm)
      expect_equal(returned[ri, ci], x1)  
    }  
  }
})

# checkEncoding - can check by inspection, don't necessarily need a test for this?
#   - unless the inputs are invalid, in which case will do this one in the next
#   iteration of test development (malicious/incorrect inputs)
# test_that("checkEncoding()", {})


# formatData - can check by inspection, don't necessarily need a test for this?
#   - unless the inputs are invalid, in which case will do this one in the next 
#   iteration of test development (malicious/incorrect inputs)
# test_that("formatData()", {formatData(dm, c("dr1","dr2","Temp.Dist"))})

test_results[[6]] <- test_that("processedStrains() - requirements test", {
  base_strains <- sample.int(200, 1) %>% fakeStrains()
  returned <- processedStrains(base_strains)
  
  strain_data <- base_strains %>% 
    mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
    select(colnames(returned$strain_data))
  
  assignments <- strain_data %>% select(Date, Latitude, Longitude) %>% 
    unique() %>% rownames_to_column("dr") %>% as.data.table()
  
  dr_matches <- left_join(strain_data, assignments, by = c("Latitude", "Longitude", "Date")) %>% 
    select(Strain, dr)
  
  expect_identical(strain_data, returned$strain_data)
  expect_identical(assignments, returned$assignments)
  expect_identical(dr_matches, returned$dr_matches)
})

test_results[[7]] <- test_that("processedStrains() - additional columns test", {
  additional_cols <- sample(c("Country","Province","City"), sample.int(3, 1)) %>% sort()
  base_strains <- sample.int(200, 1) %>% fakeStrains(., additional_cols)
  returned <- processedStrains(base_strains)
  
  strain_data <- base_strains %>% 
    mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")), 
           Location = do.call(paste, c(base_strains[additional_cols], sep = "_"))) %>% 
    select(colnames(returned$strain_data))
  
  assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
    unique() %>% rownames_to_column("dr") %>% as.data.table()
  
  dr_matches <- left_join(strain_data, assignments, by = c("Latitude", "Longitude", "Date")) %>% 
    select(Strain, dr)
  
  expect_identical(strain_data, returned$strain_data)
  expect_identical(assignments, returned$assignments)
  expect_identical(dr_matches, returned$dr_matches)
})


test_results[[8]] <- test_that("epiCohesion()", {
  ndrs <- 10
  alldrs <- sample(seq(1,200), ndrs)
  nstrains <- ceiling(sample(25:75, 1))
  nclusters <- 15
  
  drs <- sample(alldrs, nstrains, replace = TRUE)
  
  dr_vals <- data.table(dr = unique(drs), v = sample(seq(0, 1, by = 0.0001), length(unique(drs))))
  
  dr_matches <- data.table(
    th = sample(1:nclusters, nstrains, replace = TRUE),
    strain = 1:nstrains %>% as.character(), 
    dr = drs
    ) %>% left_join(., dr_vals, by = "dr") %>% 
    as_tibble() %>% arrange(th) %>% 
    mutate(both = paste(strain, dr, sep = "-"))
  
  df2 <- dr_matches %>% select(strain, dr)
  epi_melt_interim <- strainByStrainEpiMelt(dr_matches)
  
  epi_melt <- left_join(epi_melt_interim, df2, by = c("strain1" = "strain")) %>% 
    select(-strain1) %>% rename(dr_1 = dr) %>% 
    left_join(., df2, by = c("strain2" = "strain")) %>% rename(dr_2 = dr, value = Total.Dist) %>% 
    select(dr_1, dr_2, value) %>% unique()
  
  g_cuts <- dr_matches %>% group_by(th) %>% count(dr) %>% ungroup() %>% as_tibble()
  
  returned <- epiCohesion(g_cuts, epi_melt)
  g_cuts <- dr_matches %>% select(-dr, -both)
  actual <- epi_melt_interim %>% rename(value = Total.Dist) %>% strainByStrainECC(g_cuts, .)
  
  expect_equal(returned, actual)
})


test_results[[9]] <- test_that("epiCollectionByCluster()", {
  fake_obj <- fakeSecCluster(sample.int(50,1), sample.int(200,1)*2, sample.int(20,1))
  strain_raws <- left_join(fake_obj$c, fake_obj$a, by = "dr")
  
  cluster_asmts <- strain_raws %>% select(-Strain, -th) %>% unique()
  
  returned <- runDrECBC(cluster_asmts, strain_data = strain_raws, tau = 1, gamma = 0)
  actual <- runStrainECBC(strain_raws, tau = 1, gamma = 0)
  cat("\n")
  expect_equal(returned, actual)
})


# collectECCs - can check by inspection, don't necessarily need a test for this?
#   - unless the inputs are invalid, in which case will do this one in the next 
#   iteration of test development (malicious/incorrect inputs)
#   - it's basically just running epiCollectionByCluster()
# test_that("collectECCs()", {})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input ECC tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}