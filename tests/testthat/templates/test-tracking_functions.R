libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/CGM/tracking_functions.R")

test_results <- vector(length = 7) %>% 
  setNames(c("compsSet()", "noChange()", "checkEachIsolate()", "trackSingletons()", 
             "trackClusters()", "oneHeight()", "findingSneakers()"))

test_results[[1]] <- test_that("compsSet()", {
  expect_true(FALSE)
})

test_results[[2]] <- test_that("noChange()", {
  expect_true(FALSE)
})

test_results[[3]] <- test_that("checkEachIsolate()", {
  
  expect_true(FALSE)
})

test_results[[4]] <- test_that("trackSingletons()", {
  expect_true(FALSE)
})

test_results[[5]] <- test_that("trackClusters()", {
  expect_true(FALSE)
})

test_results[[6]] <- test_that("oneHeight()", {
  expect_true(FALSE)
})

test_results[[7]] <- test_that("findingSneakers()", {
  expect_true(FALSE)
})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input CGM tracking functions tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}