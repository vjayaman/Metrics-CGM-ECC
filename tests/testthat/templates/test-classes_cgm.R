libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/CGM/classes_cgm.R")

test_results <- vector(length = 11) %>% 
  setNames(c("Timedata - new()", "Timedata - set_comps()", "Timedata - flag_clusters()", 
             "Timedata - coded_status()", "Timedata - set_cnames()", "Heightdata - new()", 
             "Heightdata - clust_tracking()", "Heightdata - post_data()", 
             "Heightdata - unchanged()", "Heightdata - update_iteration()", 
             "Heightdata - reset_values()"))

test_results[[1]] <- test_that("Timedata - new()", {
  # tp1 <- Timedata$new("tp1", raw = f1, all_isolates, pad_height = ph, pad_cluster = pc)
  # self$name <- name
  # self$start()
  # self$raw <- raw
  # self$coded <- codeIsolates(raw, name, isos, pad_height, pad_cluster)
  # self$melted <- meltedIDs(raw, name, pad_height, pad_cluster)  
  expect_true(FALSE)
})

test_results[[2]] <- test_that("Timedata - set_comps()", {
  # compsSet(., toupper(self$name), indicate_progress = TRUE)  
  expect_true(FALSE)
})

test_results[[3]] <- test_that("Timedata - flag_clusters()", {
  # self$flagged <- flaggingClusters(self$comps, self$name)
  expect_true(FALSE)
})

test_results[[4]] <- test_that("Timedata - coded_status()", {
  expect_true(FALSE)
})

test_results[[5]] <- test_that("Timedata - set_cnames()", {
  expect_true(FALSE)
})


test_results[[6]] <- test_that("Heightdata - new()", {
  # self$h_after <- self$h_before <- starter
  # self$comps <- t1_comps %>% filter(tp1_h == starter) %>% arrange(tp1_h, tp1_cl)
  # self$bef <- t1_comps %>% filter(tp1_h == self$h_before) %>% 
  #   set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
  # self$results <- vector(mode = "list", length = length(hvals)) %>% set_names(hvals)  
  expect_true(FALSE)
})

test_results[[7]] <- test_that("Heightdata - clust_tracking()", {
  # self$changed <- trackClusters(self$comps, t2_comps, t2_cnames, t1_coded, t2_coded, indp)
  expect_true(FALSE)
})

test_results[[8]] <- test_that("Heightdata - post_data()", {
  # self$aft <- t1_comps %>% filter(tp1_h == self$h_after) %>% 
  #   set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))
  expect_true(FALSE)
})

test_results[[9]] <- test_that("Heightdata - unchanged()", {
  # self$same <- noChange(self$aft, self$bef, self$tracked)
  expect_true(FALSE)
})

test_results[[10]] <- test_that("Heightdata - update_iteration()", {
  # self$tracked <- bind_rows(self$changed, self$same) %>% arrange(tp1_h, tp1_cl)
  # self$results[[self$h_after]] <- self$tracked
  expect_true(FALSE)
})

test_results[[11]] <- test_that("Heightdata - reset_values()", {
  # self$h_before <- self$h_after
  # self$h_after <- NULL
  # self$bef <- self$aft %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
  # self$aft <- tibble()
  expect_true(FALSE)
})






test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input CGM classes tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}

