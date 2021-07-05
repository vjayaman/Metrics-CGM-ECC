libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/ECC/functions/ecc_functions.R")


# (3) timeTaken() function
# - Indicates length of a process in hours, minutes, and seconds, when given a name of the process 
# ("pt") and a two-element named vector with Sys.time() values named "start_time" and "end_time"

test_that("timeTaken() - null inputs", {
  stopwatch <- list("start_time" = NULL, "end_time" = NULL)
  expect_identical(timeTaken("test", stopwatch), "Neither start nor end time were collected")
  
  stopwatch$end_time <- as.character.POSIXt(Sys.time())
  expect_identical(timeTaken("test", stopwatch), "Start time was not collected.")
  
  stopwatch$start_time <- as.character.POSIXt(Sys.time())
  stopwatch$end_time <- NULL
  expect_identical(timeTaken("test", stopwatch), "End time was not collected.")
})

test_that("timeTaken() - time inputs", {
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
test_that("outputMessages()", {
  expect_identical(outputMessages(""), cat(paste0("\n","")))
  expect_output(outputMessages("Test"), "\nTest")
  expect_silent(outputMessages())
})


# (10)
# checkEncoding(fp)
# readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()

# (1) Incorporating the allele data with the epidemiological data 
# epiCollectionByCluster(strain_data, tau, gamma, transformed_dists, tpx, cluster_x)

# (2) Note that we sum the values for both directions e.g. (192, 346) and (346, 192) and then divide by the total number of pairs (counting both directions). This gives the same result as if we had used only one direction - if we did this, we would divide by 2 in the numerator and the denominator --> would cancel
# avgDists(g_cuts, dm, cname, newname)

# (5) maxd <- max(logdata) # mind <- min(logdata)
# transformData2(dm, dtype, min_dist, max_dist)

# (6)
# transformData(dm, dtype)

# (7)
# formatData(dm, newnames)

# (8)
# epi_cohesion_new(g_cuts, epi_melt)

# (9) Incorporating the allele data with the epidemiological data
# epiCollection(strain_data, tau, gamma, typing_data, transformed_dists, dm_temp, dm_geo, dr_matches)

# (11)
# mergeECCs(eccs, tpx, typing_data)

# (12)
# processedStrains(base_strains)

# (13)
# collectECCs(k, m, parts, extremes, c1, fpath1, fpath2)
