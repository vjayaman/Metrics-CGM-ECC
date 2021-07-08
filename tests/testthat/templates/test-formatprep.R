libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/formatprep.R")

test_results <- vector(length = 7) %>% 
  setNames(c("checkEncoding()", "outputDetails()", "readData()", "writeData()", 
             "intClusters()", "updateStrains()", "strainsInSingletons()"))

test_results[[1]] <- test_that("checkEncoding()", {
  expect_true(FALSE)
})

test_results[[2]] <- test_that("outputDetails()", {
  expect_true(FALSE)
})

test_results[[3]] <- test_that("readData()", {
  
  expect_true(FALSE)
})

test_results[[4]] <- test_that("writeData()", {
  expect_true(FALSE)
})

test_results[[5]] <- test_that("intClusters()", {
  expect_true(FALSE)
})

test_results[[6]] <- test_that("updateStrains()", {
  expect_true(FALSE)
})

test_results[[7]] <- test_that("strainsInSingletons()", {
  expect_true(FALSE)
})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input formatprep functions tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}