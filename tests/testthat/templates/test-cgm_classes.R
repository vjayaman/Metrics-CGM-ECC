libs <- c("R6", "tibble", "testit", "testthat", "reader", "magrittr", "dplyr", 
          "data.table")
# libs <- c("optparse","readr", "reshape2","fossil","tidyr","purrr")
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/CGM/cgm_functions.R")
source("tests/testthat/helper_functions.R")

test_results <- vector(length = 11) %>% 
  setNames(c("Timedata - new()", "Timedata - set_comps()", "Timedata - flag_clusters()", 
             "Timedata - coded_status()", "Timedata - set_cnames()", "Heightdata - new()", 
             "Heightdata - clust_tracking()", "Heightdata - post_data()", 
             "Heightdata - unchanged()", "Heightdata - update_iteration()", 
             "Heightdata - reset_values()"))

# Fakes ----------------------------------------------------------------------
t1isos <- sample(1:10000, 1)
tp1 <- simulatedAssignments(ntph = 5, nisos = t1isos, nclusts = 100) %>% 
  set_colnames(c("isolate", as.character(0:(ncol(.) - 2)))) %>% as_tibble()

tp2 <- simulatedAssignments(ntph = 5, 
                            nisos = t1isos + sample(1:5000, 1), nclusts = 300) %>% 
  set_colnames(c("isolate", as.character(0:(ncol(.) - 2)))) %>% as_tibble()

all_isolates <- unique(c(tp1$isolate, tp2$isolate)) %>% as_tibble() %>% 
  set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

test_cases <- list("tp1" = tp1, "tp2" = tp2)
m <- sample(1:2, 1)

tp <- test_cases[[m]]
ntp <- names(test_cases)[m]
ph <- max(nchar(colnames(tp1)[-1]), nchar(colnames(tp2)[-1]))
pc <- tp2 %>% select(-isolate) %>% max(., tp2 %>% select(-isolate)) %>% nchar()
# ----------------------------------------------------------------------------

test_results[[1]] <- test_that("Timedata - new()", {
  tpx <- Timedata$new(name = ntp, raw = tp, isos = all_isolates, 
                      pad_height = ph, pad_cluster = pc)
  
  expect_identical(tpx$name, ntp)
  expect_identical(tpx$raw, tp)
  
  coded <- codeIsolates(tp, ntp, all_isolates, ph, pc)
  expect_identical(tpx$coded, coded)
  
  melted <- meltedIDs(tp, ntp, ph, pc)
  expect_identical(tpx$melted, melted)
  
  expect_true(all(is.null(tpx$flagged), is.null(tpx$cnames), 
                  is.null(tpx$comps), is.null(tpx$status)))
})

test_results[[2]] <- test_that("Timedata - set_comps()", {
  tpx <- Timedata$new(name = ntp, raw = tp, isos = all_isolates, 
                      pad_height = ph, pad_cluster = pc)
  tpx$set_comps()
  
  composition <- tpx$coded %>% 
    set_colnames(gsub(tpx$name, "tp", colnames(.))) %>% 
    compsSet(., toupper(tpx$name), indicate_progress = TRUE)
  
  expect_identical(tpx$comps, composition)
})

test_results[[3]] <- test_that("Timedata - flag_clusters()", {
  tpx$flag_clusters()
  flagged <- flaggingClusters(tpx$comps, tpx$name)
  expect_identical(tpx$flagged, flagged)
})

test_results[[4]] <- test_that("Timedata - coded_status()", {
  tp1_df <- Timedata$new(name = "tp1", raw = tp1, isos = all_isolates, 
                         pad_height = ph, pad_cluster = pc)
  
  tp2_df <- Timedata$new(name = "tp2", raw = tp2, isos = all_isolates, 
                         pad_height = ph, pad_cluster = pc)
  
  expect_null(tp2_df$status)
  
  novels <- setdiff(tp2_df$coded$isolate, tp1_df$coded$isolate)
  tp2_df$coded_status(novels)
  
  tp2_df$status %>% 
    filter(isolate %in% tp2_df$coded$isolate) %>% 
    filter(!(isolate %in% tp1_df$coded$isolate)) %>% 
    pull(status) %>% unique() %>% expect_equal(., "novs")
  
  tp2_df$status %>% 
    filter(isolate %in% tp1_df$coded$isolate) %>% 
    pull(status) %>% unique() %>% is.na(.) %>% expect_true()
})

test_results[[5]] <- test_that("Timedata - set_cnames()", {
  m <- sample(1:2, 1)
  tp <- test_cases[[m]]
  ntp <- names(test_cases)[m]
  
  tpx <- Timedata$new(name = ntp, raw = tp, isos = all_isolates, 
                      pad_height = ph, pad_cluster = pc)
  expect_null(tpx$cnames)
  
  tpx$set_cnames()
  x <- tpx$coded %>% pull(grep("h", colnames(tpx$coded), value = TRUE)) %>% unique() %>% sort()
  expect_identical(tpx$cnames, x)
})

# starter height: 0
test_results[[6]] <- test_that("Heightdata - new()", {
  tp1_df <- Timedata$new("tp1", raw = tp1, all_isolates, pad_height = ph, pad_cluster = pc)$
    set_comps()$flag_clusters()
  
  tp2_df <- Timedata$new("tp2", raw = tp2, all_isolates, pad_height = ph, pad_cluster = pc)$
    set_comps()$set_cnames()

  novels <- setdiff(tp2_df$coded$isolate, tp1_df$coded$isolate)
  
  tp2_df$comps <- tp2_df$coded %>% filter(isolate %in% novels) %>% 
    group_by(tp2_id) %>% 
    summarise(num_novs = n(), .groups = "drop") %>% 
    left_join(tp2_df$comps, ., by = "tp2_id") %>% 
    mutate(num_novs = ifelse(is.na(num_novs), 0, num_novs))
  
  tp2_df$flag_clusters()$coded_status(novels)
  tp1_df$coded_status(novels)
  
  heights <- "0"
  
  hx <- Heightdata$new(starter = heights[1], t1_comps = tp1_df$comps, hvals = heights)

  expect_equal(hx$h_after, heights[1])
  expect_equal(hx$h_before, heights[1])
  
  tp1_df$comps %>% filter(tp1_h == heights[1]) %>% arrange(tp1_h, tp1_cl) %>% 
    expect_identical(hx$comps, .)
  
  tp1_df$comps %>% filter(tp1_h == heights[1]) %>% 
    set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
  
  expect_null(unlist(hx$results))
  expect_identical(names(hx$results), heights[1])
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