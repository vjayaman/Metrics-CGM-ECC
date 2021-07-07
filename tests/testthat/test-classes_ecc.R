libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/ECC/functions/classes_ecc.R")

test_results <- vector(length = 4) %>% 
  setNames(c("Timepoint - new()", "Timepoint - readTyping()", 
             "Timepoint - Process()", "Timepoint - listHeights()"))

test_results[[1]] <- test_that("Timepoint - new()", {
  # self$filepath <- fpath
  # self$name <- toupper(name)
  # self$readTyping()
  expect_true(FALSE)
})

test_results[[2]] <- test_that("Timepoint - readTyping()", {
  # self$filedata <- read.table(self$filepath, header = TRUE, sep = "\t", row.names = 1, 
  #                             check.names = FALSE, quote = "", stringsAsFactors = FALSE, 
  #                             fileEncoding = checkEncoding(self$filepath))
  expect_true(FALSE)
})

test_results[[3]] <- test_that("Timepoint - Process()", {
  # self$proc <- self$filedata %>% 
  #   select(hx$h) %>% 
  #   set_colnames(paste0(self$name, "_T", colnames(.))) %>% 
  #   rownames_to_column("Strain") %>% as_tibble() %>% 
  #   add_column(tpx = 1)
  # colnames(self$proc)[colnames(self$proc) == "tpx"] <- self$name
  expect_true(FALSE)
})

test_results[[4]] <- test_that("Timepoint - listHeights()", {
  # self$height_list <- lapply(1:nrow(hx), function(i) {
  #   self$filedata[,hx$h[i],drop=FALSE] %>% set_colnames(hx$th[i])
  # }) %>% set_names(paste0(self$name, "_", hx$th))
  expect_true(FALSE)
})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input ECC classes tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}

