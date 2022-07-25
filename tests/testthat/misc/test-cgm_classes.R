libs <- c("R6", "tibble", "testit", "testthat", "reader", "magrittr", "dplyr", 
          "data.table", "progress")
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/CGM/cgm_functions.R")
source("tests/testthat/helper_functions.R")

test_results <- vector(length = 12) %>% 
  setNames(c("Timedata - new()", "Timedata - set_comps()", "Timedata - flag_clusters()", 
             "Timedata - coded_status()", "Timedata - set_cnames()", "Heightdata - new()", 
             "Heightdata - clust_tracking() - setup - base case", 
             "Heightdata - clust_tracking() - check tracking - base case", 
             "Heightdata - update_iteration() - base case", 
             "Heightdata - post_data() - first iteration", 
             "Heightdata - unchanged() - first iteration", "Heightdata - reset_values()"))

# Fakes ----------------------------------------------------------------------
tp2_assignments <- simulatedAssignments(ntph = 5, 
                            nisos = sample(3000:10000, 1), 
                            nclusts = sample(200:600, 1)) %>% 
  set_colnames(c("isolate", as.character(0:(ncol(.) - 2)))) %>% as_tibble()

tp1_assignments <- tp2_assignments %>% filter(!(isolate %in% sample(tp2_assignments$isolate, sample(1000:2500, 1))))

all_isolates <- unique(c(tp1_assignments$isolate, tp2_assignments$isolate)) %>% as_tibble() %>% 
  set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

test_cases <- list("tp1" = tp1_assignments, "tp2" = tp2_assignments)
m <- sample(1:2, 1)

tp <- test_cases[[m]]
ntp <- names(test_cases)[m]
ph <- max(nchar(colnames(tp1_assignments)[-1]), nchar(colnames(tp2_assignments)[-1]))
pc <- tp2_assignments %>% select(-isolate) %>% max(., tp2_assignments %>% select(-isolate)) %>% nchar()
# ----------------------------------------------------------------------------

tpx <- Timedata$new(name = ntp, raw = tp, isos = all_isolates, pad_height = ph, 
                      pad_cluster = pc, msg = FALSE, ind_prog = FALSE)

test_results[[1]] <- test_that("Timedata - new()", {
  expect_identical(tpx$name, ntp)
  expect_identical(tpx$raw, tp)
  
  coded <- codeIsolates(tp, ntp, all_isolates, ph, pc)
  expect_identical(tpx$coded, coded)
  
  melted <- meltedIDs(tp, ntp, ph, pc)
  expect_identical(tpx$melted, melted)
  
  expect_true(all(is.null(tpx$flagged), is.null(tpx$cnames), 
                  is.null(tpx$comps), is.null(tpx$status)))
})

tpx$set_comps()

test_results[[2]] <- test_that("Timedata - set_comps()", {
  composition <- tpx$coded %>% 
    set_colnames(gsub(tpx$name, "tp", colnames(.))) %>% 
    compsSet(., toupper(tpx$name), indicate_progress = FALSE)
  expect_identical(tpx$comps, composition)
})

tpx$flag_clusters()

test_results[[3]] <- test_that("Timedata - flag_clusters()", {
  flagged <- flaggingClusters(tpx$comps, tpx$name)
  expect_identical(tpx$flagged, flagged)
})

tp1 <- Timedata$new(name = "tp1", raw = tp1_assignments, isos = all_isolates, pad_height = ph, 
                    pad_cluster = pc, msg = FALSE, ind_prog = FALSE)
tp2 <- Timedata$new(name = "tp2", raw = tp2_assignments, isos = all_isolates, pad_height = ph, 
                    pad_cluster = pc, msg = FALSE, ind_prog = FALSE)
  
novels <- setdiff(tp2$coded$isolate, tp1$coded$isolate)

test_results[[4]] <- test_that("Timedata - coded_status()", {
  expect_null(tp2$status)
  tp2$coded_status(novels)
  
  tp2$status %>% 
    filter(isolate %in% tp2$coded$isolate) %>% 
    filter(!(isolate %in% tp1$coded$isolate)) %>% 
    pull(status) %>% unique() %>% expect_equal(., "novs")
  
  tp2$status %>% 
    filter(isolate %in% tp1$coded$isolate) %>% 
    pull(status) %>% unique() %>% is.na(.) %>% expect_true()
})

m <- sample(1:2, 1)
tp <- test_cases[[m]]
ntp <- names(test_cases)[m]

test_results[[5]] <- test_that("Timedata - set_cnames()", {
  tpx <- Timedata$new(name = ntp, raw = tp, isos = all_isolates, pad_height = ph, 
                      pad_cluster = pc, msg = FALSE, ind_prog = FALSE)
  expect_null(tpx$cnames)
  
  tpx$set_cnames()
  x <- tpx$coded %>% pull(grep("h", colnames(tpx$coded), value = TRUE)) %>% unique() %>% sort()
  expect_identical(tpx$cnames, x)
})

tp1 <- tp1$set_comps()$flag_clusters()
tp2 <- tp2$set_comps()$set_cnames()

tp2$comps <- tp2$coded %>% filter(isolate %in% novels) %>% group_by(tp2_id) %>% 
  summarise(num_novs = n(), .groups = "drop") %>% left_join(tp2$comps, ., by = "tp2_id") %>% 
  mutate(num_novs = ifelse(is.na(num_novs), 0, num_novs))

heights <- colnames(tp1$raw)[-1] %>% paste0(., collapse = ",") %>% 
    strsplit(., split = ",") %>% unlist()

# starter height: 0
test_results[[6]] <- test_that("Heightdata - new()", {
  tp2$flag_clusters()$coded_status(novels)
  tp1$coded_status(novels)
  
  hx <- Heightdata$new(starter = heights[1], t1_comps = tp1$comps, hvals = heights)

  expect_equal(hx$h_after, heights[1])
  expect_equal(hx$h_before, heights[1])
  
  tp1$comps %>% filter(tp1_h == heights[1]) %>% arrange(tp1_h, tp1_cl) %>% 
    expect_identical(hx$comps, .)
  
  tp1$comps %>% filter(tp1_h == heights[1]) %>% 
    set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef")) %>% 
    expect_identical(hx$bef, .)
  
  expect_null(unlist(hx$results))
  expect_identical(names(hx$results), heights)
  expect_equal(nrow(hx$aft), 0)
  expect_equal(ncol(hx$aft), 0)
})
  
tp2$flag_clusters()$coded_status(novels)
tp1$coded_status(novels)

hx <- Heightdata$new(starter = heights[1], t1_comps = tp1$comps, hvals = heights)
hx$clust_tracking(tp2$comps, tp2$cnames, tp1$coded, tp2$coded, ind_prog = FALSE)

changed <- hx$changed %>% 
    mutate(both1 = paste0(tp1_h, "-", tp1_cl), 
           both2 = paste0(tp2_h, "-", tp2_cl))


test_results[[7]] <- test_that("Heightdata - clust_tracking() - setup - base case", {
  changed %>% pull(tp1_id) %>% checkIDs(., "h", 2) %>% 
    expect_identical(., changed$tp1_h)
  changed %>% pull(tp1_id) %>% checkIDs(., "c", 3) %>% as.double() %>% 
    expect_identical(., changed$tp1_cl)
  changed %>% pull(tp1_id) %>% strsplit(., "_") %>% 
    sapply(., '[[', 1) %>% unique() %>% expect_equal(., toupper("tp1"))
  
  
  changed %>% pull(tp2_id) %>% checkIDs(., "h", 2) %>% 
    expect_identical(., changed$tp2_h)
  changed %>% pull(tp2_id) %>% checkIDs(., "c", 3) %>% as.double() %>% 
    expect_identical(., changed$tp2_cl)
  changed %>% pull(tp2_id) %>% strsplit(., "_") %>% 
    sapply(., '[[', 1) %>% unique() %>% expect_equal(., toupper("tp2"))
})

test_results[[8]] <- test_that("Heightdata - clust_tracking() - check tracking - base case", {
  tp2$raw <- tp2$raw %>% as.data.table()
  mtp2 <- tp2$raw %>% melt.data.table(id.vars = "isolate", 
                                      variable.name = "tp2_h", value.name = "tp2_cl") %>% 
    mutate(both2 = paste0(tp2_h, "-", tp2_cl)) %>% 
    arrange(tp2_h, tp2_cl)
  
  initial_pairs <- changed %>% select(tp1_h, tp1_cl, both1) %>% unique()
  
  for (i in 1:nrow(initial_pairs)) {
    h1 <- initial_pairs$tp1_h[i]
    c1 <- initial_pairs$tp1_cl[i]
    tp1_isos <- tp1$raw %>% filter(!!as.symbol(as.character(h1)) == c1) %>% pull(isolate)
    t2_clusters <- mtp2[isolate %in% tp1_isos] %>% select(-isolate) %>% unique()
    
    actual_tracked <- mtp2[both2 %in% t2_clusters$both2]
    
    number_novels <- actual_tracked[!(isolate %in% tp1$raw$isolate)] %>% 
      group_by(both2) %>% summarise(num_novs = n()) %>% 
      mutate(across(num_novs, as.double))
    
    tracked_clusters <- actual_tracked %>% group_by(both2) %>% 
      summarise(tp2_cl_size = n(), .groups = "drop") %>% 
      filter(tp2_cl_size >= length(tp1_isos)) %>% 
      left_join(., number_novels, by = "both2") %>% 
      mutate(num_novs = ifelse(is.na(num_novs), 0, num_novs)) %>% 
      as.data.table()
    
    changed %>% filter(both1 %in% initial_pairs$both1[i]) %>% 
      select(both2, tp2_cl_size, num_novs) %>% 
      as.data.table() %>% 
      expect_equal(., tracked_clusters)
  }
})

test_results[[9]] <- test_that("Heightdata - update_iteration() - base case", {
  expect_null(hx$results[[heights[1]]])
  expect_null(hx$tracked)
  
  hx$update_iteration()
  
  bind_rows(hx$changed, hx$same) %>% arrange(tp1_h, tp1_cl) %>% 
    expect_identical(., hx$tracked)
  
  expect_identical(hx$results[[heights[1]]], hx$tracked)
})

test_results[[10]] <- test_that("Heightdata - post_data() - first iteration", {
  h2 <- heights[-1][1]
  hx$h_after <- h2
  
  expect_equal(nrow(hx$aft), 0)
  expect_equal(ncol(hx$aft), 0)
  
  hx$post_data(tp1$comps)
  
  tp1$comps %>% filter(tp1_h == h2) %>% 
    set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft")) %>% 
    expect_identical(., hx$aft)
})

test_results[[11]] <- test_that("Heightdata - unchanged() - first iteration", {
  expect_equal(nrow(hx$same), 0)
  expect_equal(ncol(hx$same), 0)
  
  hx$unchanged()
  
  part1 <- tp1$raw %>% select(isolate, hx$h_before) %>% as.data.table() %>% 
      arrange(!!as.symbol(hx$h_before), isolate) %>% set_colnames(c("isolate", "hbef"))
  
  part2 <- tp1$raw %>% select(isolate, hx$h_after) %>% as.data.table() %>% 
      arrange(!!as.symbol(hx$h_after), isolate) %>% set_colnames(c("isolate", "haft"))
  
  a1 <- tibble()
  for (initial_cluster in unique(part1$hbef)) {
    isos1 <- part1[hbef == initial_cluster] %>% pull(isolate) %>% sort()
    cl2 <- part2[isolate %in% isos1] %>% pull(haft) %>% unique()
    if (length(cl2) == 1) { # cluster at TP1 found in only one cluster at TP2
      isos2 <- part2[haft == cl2] %>% pull(isolate) %>% sort()
      if (identical(isos1, isos2)){
        a1 <- tibble(h_bef = as.integer(hx$h_before), cl_bef = initial_cluster, 
                     h_aft = as.integer(hx$h_after), cl_aft = cl2) %>% bind_rows(a1, .)
      }
    }
  }
  actual_nochange <- hx$tracked %>% 
    left_join(a1, ., by = c("h_bef" = "tp1_h", "cl_bef" = "tp1_cl")) %>% 
    select(-h_bef, -cl_bef) %>% rename(tp1_h = h_aft, tp1_cl = cl_aft) %>% 
    select(-tp1_id) %>% 
    newID(., "tp1", "tp1_h", "tp1_cl", ph, pc) %>% rename(tp1_id = id) %>% 
    select(colnames(hx$same)) %>% arrange(tp1_h, tp1_cl)
  
  hx$same %>% arrange(tp1_h, tp1_cl) %>% expect_identical(actual_nochange, .)
})


test_results[[12]] <- test_that("Heightdata - reset_values()", {
  expect_false(hx$h_after == hx$h_before)
  expect_false(identical(hx$bef, hx$aft))
  
  hafter <- hx$h_after
  haft <- hx$aft %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
  
  hx$bef$h_bef %>% unique() %>% expect_equal(., as.integer(hx$h_before))
  expect_gt(nrow(hx$aft), 0)
  
  hx$reset_values()
  
  expect_equal(hx$h_before, hafter)
  hx$bef$h_bef %>% unique() %>% expect_equal(., as.integer(hafter))
  expect_null(hx$h_after)
  expect_equal(nrow(hx$aft), 0)
  expect_identical(hx$bef, haft)
})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input CGM classes tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
  Sys.time() %>% gsub(" ", "-", .) %>% gsub(":", ".", .) %>% 
    paste0("tests/testthat/logs/tp1-", ., ".Rds") %>% 
    saveRDS(tp1, .)
  
  Sys.time() %>% gsub(" ", "-", .) %>% gsub(":", ".", .) %>% 
    paste0("tests/testthat/logs/tp2-", ., ".Rds") %>% 
    saveRDS(tp2, .)
}