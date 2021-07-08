libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

files <- paste0("scripts/ECC/functions/") %>% list.files(., full.names = TRUE)
invisible(sapply(files, source))

test_results <- vector(length = 12) %>% 
  setNames(c("Part 1", "Part 2", "Part 3", "Part 4", "Part 5", "Part 6"))

# Part 1
# combos <- params$trio %>% strsplit(., "-") %>% unlist()
# z <- vector("list", length = length(combos)) %>% set_names(combos)
test_that("Part 1", {})

# Part 2
# hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))
test_that("Part 2", {})

# Part 3
# tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
test_that("Part 3", {})

# Part 4
# tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
test_that("Part 4", {})

# Part 5
typing_data <- tp1$height_list %>% append(tp2$height_list)
test_that("Part 5", {})

# Part 6
# m <- read_tsv(params$strains) %>% processedStrains()
test_that("Part 6", {})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input CGM collection tests have passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}