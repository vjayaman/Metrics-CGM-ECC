libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/MergeFunctions/type_handling.R")

test_results <- vector(length = 11) %>% 
  setNames(c("checkTypes()", "colsToInherit()", "inheritTP1Data()", 
             "type2Inheritance()", "type3Inheritance()", "writeData()", 
             "checkEncoding()", "readData()", "getAverage()", 
             "getDistCols()", "replaceDistName()"))

test_results[[1]] <- test_that("checkTypes()", {
  expect_true(FALSE)
})

test_results[[2]] <- test_that("colsToInherit()", {
  expect_true(FALSE)
})

test_results[[3]] <- test_that("inheritTP1Data()", {
  
  expect_true(FALSE)
})

test_results[[4]] <- test_that("type2Inheritance()", {
  expect_true(FALSE)
})

test_results[[5]] <- test_that("type3Inheritance()", {
  expect_true(FALSE)
})

test_results[[6]] <- test_that("writeData()", {
  expect_true(FALSE)
})

test_results[[7]] <- test_that("checkEncoding()", {
  expect_true(FALSE)
})

test_results[[8]] <- test_that("readData()", {
  expect_true(FALSE)
})

test_results[[9]] <- test_that("getAverage()", {
  expect_true(FALSE)
})

test_results[[10]] <- test_that("getDistCols()", {
  expect_true(FALSE)
})

test_results[[11]] <- test_that("replaceDistName()", {
  expect_true(FALSE)
})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input type handling function (from the merging section) tests have passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}