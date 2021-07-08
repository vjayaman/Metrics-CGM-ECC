libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/CGM/formatting_cgm.R")

test_results <- vector(length = 14) %>% 
  setNames(c("addingType()", "flaggingClusters()", "timeTaken()", 
             "outputDetails()", "newID()", "meltData()", "factorToInt()", 
             "codeIsolates()", "meltedSizing()", "checkEncoding()", 
             "readBaseData()", "padCol()", "meltedIDs()", "convertAndSave()"))

test_results[[1]] <- test_that("addingType()", {
  expect_true(FALSE)
})

test_results[[2]] <- test_that("flaggingClusters()", {
  expect_true(FALSE)
})

test_results[[3]] <- test_that("timeTaken()", {

  expect_true(FALSE)
})

test_results[[4]] <- test_that("outputDetails()", {
  expect_true(FALSE)
})

test_results[[5]] <- test_that("newID()", {
  expect_true(FALSE)
})

test_results[[6]] <- test_that("meltData()", {
  expect_true(FALSE)
})

test_results[[7]] <- test_that("factorToInt()", {
  expect_true(FALSE)
})

test_results[[8]] <- test_that("codeIsolates()", {
  expect_true(FALSE)
})

test_results[[9]] <- test_that("meltedSizing()", {
  expect_true(FALSE)
})

test_results[[10]] <- test_that("checkEncoding()", {
  expect_true(FALSE)
})

test_results[[11]] <- test_that("readBaseData()", {
  expect_true(FALSE)
})

test_results[[11]] <- test_that("padCol()", {
  expect_true(FALSE)
})

test_results[[11]] <- test_that("meltedIDs()", {
  expect_true(FALSE)
})

test_results[[11]] <- test_that("convertAndSave()", {
  expect_true(FALSE)
})



test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input CGM formatting function tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}