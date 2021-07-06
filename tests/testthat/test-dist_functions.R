libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/ECC/functions/dist_functions.R")

# INITIALLY, TESTING FUNCTIONS FOR CORRECT INPUTS ONLY
# WILL ADD TESTS FOR MALICIOUS/INCORRECT INPUTS AFTERWARDS

# FUNCTIONS for checking sizes and making fakes, stubs, mocks, etc.
checkSizes <- function(res, sz) {
  lapply(1:length(res), function(i) sum(res[[i]]$n) <= sz) %>% unlist() %>% all() %>% return()  
}

fakeRaw <- function(ndrs, n, nc) {
  data.table(th = sample(1:nc, n, replace = TRUE) %>% sort(), 
             dr = rep(1:ndrs, ceiling(n/ndrs))[1:n]) %>% 
    rownames_to_column("Strain") %>% return()
}

fakeSecCluster <- function(ndrs, n, nc) {
  testcase <- fakeRaw(ndrs, n, nc)
  nunique <- length(unique(testcase$dr))
  
  assignments <- data.table(
    dr = testcase$dr %>% unique(), 
    Date = sample(seq(as.Date('2019/01/01'), as.Date('2021/01/01'), by="day"), nunique), 
    Latitude = sample(seq(25.00000, 65.00000, by = 0.00001), nunique), 
    Longitude = sample(seq(-100.00000, 125.00000, by = 0.00001), nunique))
  
  # fake sectionClusters() result object:
  parts <- formatForSectioning(testcase, 5) %>% sectionTypingData() %>% 
    list("drs" = testcase, "results" = .)
  return(list("a" = assignments, "b" = parts, "c" = testcase))
}

# FAKES - used in the following tests
# In this fake object, there are 12 drs, 60 strains, and 10 clusters
fake_obj <- fakeSecCluster(12, 60, 10)
assignments <- fake_obj$a
results <- fake_obj$b$results

# TESTS - part 1/3 - checking results with correct inputs -----------------------------------

test_that("formatForSectioning()", {
  testcase <- data.table(th = sample(1:5, 30, replace = TRUE) %>% sort(), 
                         dr = sample(1:10, 30, replace = TRUE) %>% sort()) %>% 
    rownames_to_column("Strain")
  
  x <- formatForSectioning(testcase, 1)
  maxdrs <- max(x[[1]]$n)
  
  x1 <- formatForSectioning(testcase, maxdrs - 1)
  expect_equal(names(x1[2]), "dr")
  
  x2 <- formatForSectioning(testcase, maxdrs + 1)
  expect_equal(names(x2[2]), "ms")
})


test_that("sectionTypingData()", {
  testcase <- data.table(th = sample(1:5, 30, replace = TRUE) %>% sort(), 
                         dr = sample(1:10, 30, replace = TRUE) %>% sort()) %>% 
    rownames_to_column("Strain")
  
  a1 <- testcase %>% formatForSectioning(., 1)
  sectionTypingData(a1) %>% checkSizes(., a1[2]) %>% expect_true()
  
  a2 <- testcase %>% formatForSectioning(., 5)
  sectionTypingData(a2) %>% checkSizes(., a2[2]) %>% expect_true()
  
  a3 <- testcase %>% formatForSectioning(., 30)
  sectionTypingData(a3) %>% checkSizes(., a3[2]) %>% expect_true()
})


# no test included for sectionClusters(), just a particular use case of sectionTypingData()


test_that("distMatrix() - temp", {
  
  inds <- fake_obj$b$drs[["th"]] %in% pull(results[[1]], "th")
  cluster_x <- fake_obj$b$drs[inds,-"Strain"]
  cluster_asmts <- assignments[dr %in% pull(cluster_x, dr)]
  
  dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
  
  # check the distances manually - can do this since we have a relatively small distance matrix
  pairs <- expand.grid(cluster_asmts$dr, cluster_asmts$dr)
  for (i in 1:nrow(pairs)) {
    dr1 <- cluster_asmts[dr == pairs$Var1[i]]$Date
    dr2 <- cluster_asmts[dr == pairs$Var2[i]]$Date
    expect_equal(dm_temp[as.character(pairs$Var1[i]), as.character(pairs$Var2[i])], 
                 abs(as.vector(difftime(dr1, dr2, units = "days"))))
  }
})


test_that("distMatrix() - geo", {
  
  inds <- fake_obj$b$drs[["th"]] %in% pull(results[[1]], "th")
  cluster_x <- fake_obj$b$drs[inds,-"Strain"]
  cluster_asmts <- assignments[dr %in% pull(cluster_x, dr)]
  
  dm_geo <- cluster_asmts %>% select(dr, Longitude, Latitude) %>% 
    distMatrix(., "geo", c("Longitude", "Latitude"))
  
  # check the distances manually - can do this since we have a relatively small distance matrix
  pairs <- expand.grid(cluster_asmts$dr, cluster_asmts$dr)
  for (i in 1:nrow(pairs)) {
    x1 <- as.character(pairs$Var1[i])
    x2 <- as.character(pairs$Var2[i])
    dr1 <- cluster_asmts[dr == x1]
    dr2 <- cluster_asmts[dr == x2]
    
    expect_equal(dm_geo[x1, x2], 
                 deg.dist(dr1$Longitude, dr1$Latitude, dr2$Longitude, dr2$Latitude))
  }
})


test_that("outputMessages()", {
  expect_identical(outputMessages(""), cat(paste0("\n","")))
  expect_output(outputMessages("Test"), "\nTest")
  expect_silent(outputMessages())
})


test_that("collectDistances()", {
  extremes <- collectDistances(save_extremes = FALSE, assignments, fake_obj$b)
  
  min_geo <- min_temp <- Inf
  max_geo <- max_temp <- -Inf
  
  for (j in 1:length(results)) {
    inds <- fake_obj$b$drs[["th"]] %in% pull(results[[j]], "th")
    cluster_x <- fake_obj$b$drs[inds,-"Strain"]
    cluster_asmts <- assignments[dr %in% pull(cluster_x, dr)]
    
    dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
    dm_geo <- cluster_asmts %>% select(dr, Longitude, Latitude) %>% 
      distMatrix(., "geo", c("Longitude", "Latitude"))
    
    min_temp <- min(min_temp, min(dm_temp))
    max_temp <- max(max_temp, max(dm_temp))
    
    # min_geo <- min(min_geo, min(dm_geo + 10)) # necessary for the transformation step?
    min_geo <- min(min_geo, min(dm_geo)) # necessary for the transformation step?
    max_geo <- max(max_geo, max(dm_geo))
  }
  
  cat("\n")
  expect_equal(extremes$temp$max, max_temp)
  expect_equal(extremes$temp$min, min_temp)
  
  expect_equal(extremes$geo$max, max_geo)
  expect_equal(extremes$geo$min, min_geo)
})


test_that("countDataReps()", {
  dr_matches <- fake_obj$c
  dr_td1 <- dr_matches %>% select(-Strain)
  g_cuts <- countDataReps(dr_td1)
  
  clusters <- unique(pull(dr_matches, "th"))
  
  for (cl in clusters) {
    a1 <- dr_matches[th == cl] %>% group_by(dr) %>% summarise(n = n()) %>% 
      as.data.table() %>% mutate(across(dr, as.character)) %>% arrange(dr)
    b1 <- g_cuts[th == cl] %>% select(-th) %>% arrange(dr)
    expect_identical(a1, b1)
  }
})


test_that("avgDists()", {
  every_cluster <- bind_rows(results)
  for (j in 1:nrow(every_cluster)) {
    cluster_x <- fake_obj$b$drs[,-"Strain"][th %in% every_cluster$th[j]]
    cluster_asmts <- assignments[dr %in% pull(cluster_x, dr)]
    g_cuts <- countDataReps(fake_obj$c[,-"Strain"])[th %in% unique(cluster_x$th)]
    
    dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
    returned_temp <- avgDists(g_cuts, dm_temp, "Temp.Dist", 
                         paste0("TP_test_", colnames(g_cuts)[1])) %>% pull(2)
    actual_temp <- sum(dm_temp)/(nrow(dm_temp)*ncol(dm_temp))
    
    dm_geo <- cluster_asmts %>% select(dr, Longitude, Latitude) %>% 
      distMatrix(., "geo", c("Longitude", "Latitude"))
    returned_geo <- avgDists(g_cuts, dm_geo, "Geog.Dist", 
                              paste0("TP_test_", colnames(g_cuts)[1])) %>% pull(2)
    actual_geo <- sum(dm_geo)/(nrow(dm_geo)*ncol(dm_geo))
    
    expect_equal(returned_temp, actual_temp)
    expect_equal(returned_geo, actual_geo)
  }
})

# TESTS - part 2/3 - checking results with malicious/incorrect inputs -----------------------
# TESTS - part 3/3 - checking results with border cases? ------------------------------------