libs <- c("R6", "testit", "magrittr", "tibble", "dplyr", "data.table", "testthat", "tidyr")
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/CGM/cgm_functions.R")
source("tests/testthat/helper_functions.R")

test_results <- vector(length = 8) %>% 
  setNames(c("addingType()", "timeTaken() - null inputs", "timeTaken() - time inputs", 
             "newID()", "codeIsolates()", "meltedIDs()", "compsSet()", "flaggingClusters()"))

# requires as input a dataset with no column "type", with columns "tp1_cl_size" and "tp2_cl_size"
test_results[[1]] <- test_that("addingType()", {
  initial_sizes <- c(sampleSizes(0, 200, 60), sampleSizes(0, 3, 40))
  df1 <- syntheticTypeSet(sample(1:4, sample(1:4, 1))) %>% 
    addingType() %>% 
    rename(TP2 = tp2_cl_size, TP1 = tp1_cl_size)

  p1 <- df1 %>% filter(type == "Type1")
  # print(paste0("Randomized set:"))
  # print(paste0("  Number of type 1 cases: ", nrow(p1)))
  if (nrow(p1) > 0) {
    expect_true(all(p1$TP1 > 2 & p1$TP2 > 2 & p1$TP1 == p1$TP2))  
  }else {
    numT1s <- df1 %>% filter(TP1 > 2 & TP2 > 2 & TP1 == TP2) %>% nrow()
    expect_equal(numT1s, 0)
  }
  
  p2 <- df1 %>% filter(type == "Type2")
  # print(paste0("  Number of type 2 cases: ", nrow(p2)))
  if (nrow(p2) > 0) {
    expect_true(all(p2$TP1 > 2 & p2$TP2 > 2 & p2$TP2 > p2$TP1))  
  }else {
    numT2s <- df1 %>% filter(TP1 > 2 & TP2 > 2 & TP2 > TP1) %>% nrow()
    expect_equal(numT2s, 0)
  }
  
  p3 <- df1 %>% filter(type == "Type3")
  # print(paste0("  Number of type 3 cases: ", nrow(p3)))
  if (nrow(p3) > 0) {
    expect_true(all(p3$TP1 < 3 & p3$TP2 > 2))
  }else {
    numT3s <- df1 %>% filter(TP1 < 3 & TP2 > 2) %>% nrow()
    expect_equal(numT3s, 0)
  }
  
  p4 <- df1 %>% filter(type == "Type4")
  # print(paste0("  Number of type 4 cases: ", nrow(p4)))
  if (nrow(p4) > 0) {
    expect_true(all(p4$TP1 < 3 & p4$TP2 < 3))
  }else {
    numT4s <- df1 %>% filter(TP1 < 3 & TP2 < 3) %>% nrow()
    expect_equal(numT4s, 0)
  }
  
  df1 %>% filter(!(type %in% c("Type1", "Type2", "Type3", "Type4"))) %>% 
    nrow() %>% expect_equal(., 0)
})

# - Indicates length of a process in hours, minutes, and seconds, when given a name of the process 
# ("pt") and a two-element named vector with Sys.time() values named "start_time" and "end_time"
test_results[[2]] <- test_that("timeTaken() - null inputs", {
  stopwatch <- list("start_time" = NULL, "end_time" = NULL)
  expect_identical(timeTaken("test", stopwatch), "Neither start nor end time were collected")
  
  stopwatch$end_time <- as.character.POSIXt(Sys.time())
  expect_identical(timeTaken("test", stopwatch), "Start time was not collected.")
  
  stopwatch$start_time <- as.character.POSIXt(Sys.time())
  stopwatch$end_time <- NULL
  expect_identical(timeTaken("test", stopwatch), "End time was not collected.")
})

test_results[[3]] <- test_that("timeTaken() - time inputs", {
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


# outputDetails() - can be checked by inspection
#   - will need to add tests for incorrect/malicious inputs later on
# test_results[[]] <- test_that("outputDetails()", {expect_true(FALSE)})

# checkEncoding() - can be checked by inspection
#   - will need to add tests for incorrect/malicious inputs later on
# test_results[[]] <- test_that("checkEncoding()", {expect_true(FALSE)})

# readBaseData() - can be checked by inspection
#   - will need to add tests for incorrect/malicious inputs later on
# test_results[[]] <- test_that("readBaseData()", {expect_true(FALSE)})

# meltData() - can be checked by inspection
#   - will need to add tests for incorrect/malicious inputs later on
# test_results[[]] <- test_that("meltData()", {expect_true(FALSE)})

# factorToInt() - can be checked by inspection
#   - will need to add tests for incorrect/malicious inputs later on
# test_results[[]] <- test_that("factorToInt()", {expect_true(FALSE)})

# padCol() - can be checked by inspection - very similar to newID()
#   - will need to add tests for incorrect/malicious inputs later on
# For padding height and cluster columns with h0..0.., and c0..0.., respectively
# test_results[[]] <- test_that("padCol()", {expect_true(FALSE)})


# Given a dataframe df, two column names c1 and c2 (height and cluster respectively) and a new
# ID prefix tpx (e.g. "tp1"), creates an ID column and adds to df before returning df
test_results[[4]] <- test_that("newID()", {
  df <- tibble(height1 = sample(seq(0,100), 5), clusters1 = sample(seq(0,100), 5))
  
  tpx <- sample(c("tp1", "tp2"), 1)
  ph <- sample(0:100,1)
  pc <- sample(0:100,1)
  with_id <- newID(df, tpx, "height1", "clusters1", ph, pc)
  
  sapply(strsplit(with_id$id, "_"), "[[", 1) %>% unique() %>% expect_match(., toupper(tpx))
  
  returned_heights <- sapply(strsplit(with_id$id, "_"), "[[", 2)
  actual_heights <- paste0("%0", ph, "d") %>% sprintf(., df$height1) %>% paste0("h", .)
  expect_equal(returned_heights, actual_heights)
  
  returned_clusters <- sapply(strsplit(with_id$id, "_"), "[[", 3)
  actual_clusters <- paste0("%0", pc, "d") %>% sprintf(., df$clusters1) %>% paste0("c", .)
  expect_equal(returned_clusters, actual_clusters)
})


# -------------------------------------------------------------------------------------------------
# Data setup --------------------------------------------------------------------------------------
t1isos <- sample(1:10000, 1)
tp1 <- simulatedAssignments(ntph = 5, nisos = t1isos, nclusts = 100) %>% 
  set_colnames(c("isolate", as.character(0:(ncol(.) - 2)))) %>% as_tibble()
  
tp2 <- simulatedAssignments(ntph = 5, nisos = t1isos + sample(1:5000, 1), nclusts = 300) %>% 
  set_colnames(c("isolate", as.character(0:(ncol(.) - 2)))) %>% as_tibble()

all_isolates <- unique(c(tp1$isolate, tp2$isolate)) %>% as_tibble() %>% 
  set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

ph <- max(nchar(colnames(tp1)[-1]), nchar(colnames(tp2)[-1]))
pc <- tp2 %>% select(-isolate) %>% max(., tp2 %>% select(-isolate)) %>% nchar()

test_cases <- list("tp1" = tp1, "tp2" = tp2)
# -------------------------------------------------------------------------------------------------  
# -------------------------------------------------------------------------------------------------


test_results[[5]] <- test_that("codeIsolates()", {
  m <- sample(1:2, 1)
  heights <- "0"
  tp <- test_cases[[m]]
  ntp <- names(test_cases)[m]
  
  coded_isos <- codeIsolates(df = tp, tpx = ntp, all_iso = all_isolates, ph, pc) %>% 
    set_colnames(gsub(ntp, "tp", colnames(.)))
  
  coded_isos$tp_id %>% strsplit(., "_") %>% sapply(., '[[', 2) %>% 
    gsub("h", "", .) %>% as.integer() %>% 
    expect_identical(., coded_isos$tp_h)
  
  coded_isos$tp_id %>% strsplit(., "_") %>% sapply(., '[[', 3) %>% 
    gsub("c", "", .) %>% as.double() %>% 
    expect_identical(., coded_isos$tp_cl)
  
  all_isolates %<>% mutate(isolate = paste0("-", num_isolate, "-"))
  coded_isos %>% select(-tp_id) %>% left_join(., all_isolates, by = "isolate") %>% 
    select(char_isolate, tp_h, tp_cl) %>% rename(isolate = char_isolate) %>% 
    dcast(., isolate ~ tp_h, value.var = "tp_cl") %>% as_tibble() %>% 
    expect_identical(., tp)
})

test_results[[6]] <- test_that("meltedIDs()", {
  m <- sample(1:2, 1)
  
  melted_results <- meltedIDs(test_cases[[m]], names(test_cases)[m], ph, pc) %>% 
    set_colnames(c("isolate", "tp_h", "tp_cl", "tp_id"))

  melted_results$tp_id %>% strsplit(., split = "_") %>% 
    sapply(., '[[', 2) %>% gsub("h", "", .) %>% as.integer() %>% 
    as.character() %>% expect_identical(., melted_results$tp_h)
  
  melted_results$tp_id %>% strsplit(., split = "_") %>% 
    sapply(., '[[', 3) %>% gsub("c", "", .) %>% as.double() %>% 
    expect_identical(., melted_results$tp_cl)
  
  melted_results %>% as.data.table() %>% 
    dcast(., isolate ~ tp_h, value.var = "tp_cl") %>% 
    as_tibble() %>% expect_identical(., test_cases[[m]])
})


test_results[[7]] <- test_that("compsSet()", {
  m <- sample(1:2, 1)

  tp <- test_cases[[m]]
  ntp <- names(test_cases)[m]
  coded <- codeIsolates(tp, tpx = ntp, all_iso = all_isolates, ph, pc)

  composition <- coded %>% set_colnames(gsub(ntp, "tp", colnames(.))) %>% 
    compsSet(., toupper(ntp), indicate_progress = FALSE) %>% 
    set_colnames(gsub(ntp, "tp", colnames(.)))
  
  checkIDs(composition$tp_id, "h", 2) %>% expect_identical(., composition$tp_h)
  checkIDs(composition$tp_id, "c", 3) %>% as.double() %>% expect_identical(., composition$tp_cl)
  
  all_isolates %<>% rename(isolate = char_isolate) %>% mutate(code = paste0("-", num_isolate, "-"))
  
  a2 <- left_join(tp, all_isolates, by = "isolate") %>% select(-num_isolate) %>% as.data.table() %>% 
    melt.data.table(., id.vars = c("isolate", "code"), variable.name = "tp_h", value.name = "tp_cl")
  pairs <- a2 %>% select(tp_h, tp_cl) %>% unique()

  for (i in 1:nrow(pairs)) {
    coded_actual <- a2[tp_h %in% pairs$tp_h[i] & tp_cl %in% pairs$tp_cl[i]] %>% pull(code) %>% 
      gsub("-", "", .) %>% as.integer() %>% sort()
    cluster_x <- composition %>% filter(tp_h %in% pairs$tp_h[i] & tp_cl %in% pairs$tp_cl[i])
    coded_returned <- cluster_x %>% pull(composition) %>% strsplit(., ",") %>% 
      unlist() %>% gsub("-", "", .) %>% as.integer() %>% sort()
    expect_identical(coded_returned, coded_actual)
    expect_equal(cluster_x$tp_cl_size, length(coded_actual))
    expect_equal(cluster_x$tp_cl_size, length(coded_returned))
  }
})


test_results[[8]] <- test_that("flaggingClusters()", {
  m <- sample(1:2, 1)

  tp <- test_cases[[m]]
  ntp <- names(test_cases)[m]
  coded <- codeIsolates(tp, tpx = ntp, all_iso = all_isolates, ph, pc)

  composition <- coded %>% set_colnames(gsub(ntp, "tp", colnames(.))) %>% 
    compsSet(., toupper(ntp), indicate_progress = FALSE) %>% 
    set_colnames(gsub(ntp, "tp", colnames(.)))
  
  flagged <- flaggingClusters(composition, ntp)
  
  checkIDs(flagged$tp_id, "h", 2) %>% expect_identical(., flagged$tp_h)
  checkIDs(flagged$tp_id, "c", 3) %>% as.double() %>% expect_identical(., flagged$tp_cl)
  
  sizes <- tp %>% as.data.table() %>% melt.data.table(id.vars = "isolate") %>% 
    mutate(id = paste0(variable, "-", value)) %>% pull(id) %>% table() %>% 
    as.data.frame() %>% set_colnames(c("id", "tp_cl_size")) %>% as.data.table() %>% 
    separate(id, c("tp_h", "tp_cl"), "-") %>% 
    mutate(across(tp_h, as.integer)) %>% mutate(across(tp_cl, as.double)) %>% 
    arrange(tp_h, tp_cl)

  flagged %>% select(tp_h, tp_cl, tp_cl_size) %>% as.data.table() %>% 
    arrange(tp_h, tp_cl) %>% expect_equal(., sizes)
  
  mtp <- tp %>% as.data.table() %>% 
    melt.data.table(id.vars = "isolate", variable.factor = FALSE, variable.name = "tp_h", 
                    value.name = "tp_cl") %>% 
    mutate(across(tp_h, as.integer))
  mtp$isolate %<>% gsub("isolate", "", .) %>% as.integer()
  mtp %<>% arrange(tp_h, tp_cl, isolate)
  
  pairs <- mtp %>% select(-isolate) %>% unique() # ordered clusters
  for (i in 1:nrow(pairs)) {
    original <- mtp[tp_h %in% pairs$tp_h[i] & tp_cl %in% pairs$tp_cl[i]] %>% pull(isolate)

    clusters <- mtp[isolate %in% original] %>% select(tp_h, tp_cl) %>% unique()
    flagged_actual <- mtp %>% right_join(., clusters, by = c("tp_h","tp_cl")) %>% 
      select(-isolate) %>% group_by(tp_h, tp_cl) %>% summarise(n = n(), .groups = "drop") %>% 
      filter(n == length(original)) %>% slice(c(1, n())) %>% select(-n)
    
    a1 <- flagged %>% filter(tp_h %in% pairs$tp_h[i] & tp_cl %in% pairs$tp_cl[i]) %>% 
      select(first_tp2_flag, last_tp2_flag) %>% t()
    
    tibble(tp_h = checkIDs(a1, "h", 2), tp_cl = checkIDs(a1, "c", 3) %>% as.double()) %>% 
      expect_identical(., flagged_actual)
  }
})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input CGM formatting function tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}
