libs <- c("R6", "tibble", "testit", "testthat", "reader", "magrittr", "dplyr", 
          "data.table", "progress")
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/CGM/cgm_functions.R")
source("tests/testthat/helper_functions.R")

test_results <- vector(length = 4) %>% 
  setNames(c("checkEachIsolate()", "trackSingletons()", "oneHeight()", "findingSneakers()"))

# Fakes --------------------------------------------------------------------------------------
tp2_assignments <- simulatedAssignments(ntph = 5, 
                                        nisos = sample(3000:10000, 1), 
                                        nclusts = sample(200:600, 1)) %>% 
  set_colnames(c("isolate", as.character(0:(ncol(.) - 2)))) %>% as_tibble()

tp1_assignments <- tp2_assignments %>% filter(!(isolate %in% sample(tp2_assignments$isolate, sample(1000:2500, 1))))

all_isolates <- unique(c(tp1_assignments$isolate, tp2_assignments$isolate)) %>% as_tibble() %>% 
  set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

ph <- max(nchar(colnames(tp1_assignments)[-1]), nchar(colnames(tp2_assignments)[-1]))
pc <- tp2_assignments %>% select(-isolate) %>% max(., tp2_assignments %>% select(-isolate)) %>% nchar()


tp1 <- Timedata$new(name = "tp1", raw = tp1_assignments, isos = all_isolates, pad_height = ph, 
                    pad_cluster = pc, msg = FALSE, ind_prog = FALSE)
tp2 <- Timedata$new(name = "tp2", raw = tp2_assignments, isos = all_isolates, pad_height = ph, 
                    pad_cluster = pc, msg = FALSE, ind_prog = FALSE)

novels <- setdiff(tp2$coded$isolate, tp1$coded$isolate)

tp1 <- tp1$set_comps()$flag_clusters()
tp2 <- tp2$set_comps()$set_cnames()

tp2$comps <- tp2$coded %>% filter(isolate %in% novels) %>% group_by(tp2_id) %>% 
  summarise(num_novs = n(), .groups = "drop") %>% left_join(tp2$comps, ., by = "tp2_id") %>% 
  mutate(num_novs = ifelse(is.na(num_novs), 0, num_novs))

heights <- colnames(tp1$raw)[-1] %>% paste0(., collapse = ",") %>% 
  strsplit(., split = ",") %>% unlist()

tp2$flag_clusters()$coded_status(novels)
tp1$coded_status(novels)

hx <- Heightdata$new(starter = heights[1], t1_comps = tp1$comps, hvals = heights)
# --------------------------------------------------------------------------------------------


test_results[[1]] <- test_that("checkEachIsolate()", {
  multistrain <- hx$comps %>% filter(tp1_cl_size > 1)
  
  i <- sample(1:nrow(multistrain), 1)
  cluster_i <- multistrain[i,]
  
  results_i <- checkEachIsolate(cluster_i, tp2$coded, tp2$comps)
  
  h1 <- as.character(cluster_i$tp1_h)
  isos1 <- tp1$raw %>% select(isolate, all_of(h1)) %>% 
    filter(!!as.symbol(h1) == cluster_i$tp1_cl) %>% pull(isolate)
  
  mtp2 <- tp2$raw %>% as.data.table() %>% 
    melt.data.table(id.vars = "isolate", variable.name = "tp2_h", value.name = "tp2_cl") %>% 
    mutate(both2 = paste0(tp2_h, "-", tp2_cl))
  
  tp2_clusters <- mtp2 %>% filter(isolate %in% isos1) %>% select(-isolate) %>% unique()

  t2isos <- mtp2[both2 %in% tp2_clusters$both2]
  
  tp2_clusters <- t2isos %>% group_by(both2) %>% summarise(tp2_cl_size = n()) %>% 
    left_join(tp2_clusters, ., by = "both2")
  
  tp2_clusters <- t2isos[!(isolate %in% tp1$raw$isolate)] %>% group_by(both2) %>% 
    summarise(num_novs = n()) %>% 
    left_join(tp2_clusters, ., by = "both2") %>% 
    mutate(across(tp2_h, as.character)) %>% 
    mutate(across(tp2_h, as.integer)) %>% select(-both2)
  
  actual_results <- cluster_i %>% select(tp1_h, tp1_cl, tp1_cl_size) %>% bind_cols(tp2_clusters)
  expected_results <- results_i %>% select(colnames(actual_results)) %>% 
    mutate(across(num_novs, as.integer))
  expect_identical(actual_results, expected_results)
})

test_results[[2]] <- test_that("trackSingletons()", {
  t1set <- tibble()
  singletons <- hx$comps %>% filter(tp1_cl_size == 1) %>% select(tp1_id, tp1_cl_size)

  while (nrow(singletons) == 0) {
    tp2_assignments <- simulatedAssignments(ntph = 5, 
                                            nisos = sample(3000:10000, 1), 
                                            nclusts = sample(200:600, 1)) %>% 
      set_colnames(c("isolate", as.character(0:(ncol(.) - 2)))) %>% as_tibble()
    
    tp1_assignments <- tp2_assignments %>% filter(!(isolate %in% sample(tp2_assignments$isolate, sample(1000:2500, 1))))
    
    all_isolates <- unique(c(tp1_assignments$isolate, tp2_assignments$isolate)) %>% as_tibble() %>% 
      set_colnames("char_isolate") %>% rowid_to_column("num_isolate")
    
    ph <- max(nchar(colnames(tp1_assignments)[-1]), nchar(colnames(tp2_assignments)[-1]))
    pc <- tp2_assignments %>% select(-isolate) %>% max(., tp2_assignments %>% select(-isolate)) %>% nchar()
    
    
    tp1 <- Timedata$new(name = "tp1", raw = tp1_assignments, isos = all_isolates, pad_height = ph, 
                        pad_cluster = pc, msg = FALSE, ind_prog = FALSE)
    tp2 <- Timedata$new(name = "tp2", raw = tp2_assignments, isos = all_isolates, pad_height = ph, 
                        pad_cluster = pc, msg = FALSE, ind_prog = FALSE)
    
    novels <- setdiff(tp2$coded$isolate, tp1$coded$isolate)
    
    tp1 <- tp1$set_comps()$flag_clusters()
    tp2 <- tp2$set_comps()$set_cnames()
    
    tp2$comps <- tp2$coded %>% filter(isolate %in% novels) %>% group_by(tp2_id) %>% 
      summarise(num_novs = n(), .groups = "drop") %>% left_join(tp2$comps, ., by = "tp2_id") %>% 
      mutate(num_novs = ifelse(is.na(num_novs), 0, num_novs))
    
    heights <- colnames(tp1$raw)[-1] %>% paste0(., collapse = ",") %>% 
      strsplit(., split = ",") %>% unlist()
    
    tp2$flag_clusters()$coded_status(novels)
    tp1$coded_status(novels)
    
    hx <- Heightdata$new(starter = heights[1], t1_comps = tp1$comps, hvals = heights)
    
    singletons <- hx$comps %>% filter(tp1_cl_size == 1) %>% select(tp1_id, tp1_cl_size)
  }
  
  # at least one singleton cluster
  if (nrow(singletons) > 0) {
    t1set <- tp2$comps %>% select(tp2_id, tp2_cl_size, num_novs) %>%
      trackSingletons(singletons, tp1$coded, tp2$coded, .)
  }
  
  # checking: 
  i <- sample(1:nrow(singletons), 1)
  cluster_i <- hx$comps %>% filter(tp1_cl_size == 1) %>% slice(i)
  h1 <- as.character(cluster_i$tp1_h)
  isos1 <- tp1$raw %>% select(isolate, all_of(h1)) %>% 
    filter(!!as.symbol(h1) == cluster_i$tp1_cl) %>% pull(isolate)
  
  mtp2 <- tp2$raw %>% as.data.table() %>% 
    melt.data.table(id.vars = "isolate", variable.name = "tp2_h", value.name = "tp2_cl") %>% 
    mutate(both2 = paste0(tp2_h, "-", tp2_cl))
  
  tp2_clusters <- mtp2 %>% filter(isolate %in% isos1) %>% select(-isolate) %>% unique()
  
  t2isos <- mtp2[both2 %in% tp2_clusters$both2]
  
  tp2_clusters <- t2isos %>% group_by(both2) %>% summarise(tp2_cl_size = n()) %>% 
    left_join(tp2_clusters, ., by = "both2")
  
  tp2_clusters <- t2isos[!(isolate %in% tp1$raw$isolate)] %>% group_by(both2) %>% 
    summarise(num_novs = n()) %>% 
    left_join(tp2_clusters, ., by = "both2") %>% 
    mutate(across(tp2_h, as.character)) %>% 
    mutate(across(tp2_h, as.integer)) %>% select(-both2)
  
  actual_results <- cluster_i %>% select(tp1_h, tp1_cl, tp1_cl_size) %>% 
    bind_cols(tp2_clusters) %>% as.data.table()
  expected_results <- t1set %>% select(colnames(actual_results)) %>% 
    mutate(across(num_novs, as.integer))
  expect_identical(actual_results, expected_results)
})


test_results[[3]] <- test_that("oneHeight()", {
  multistrain <- hx$comps %>% filter(tp1_cl_size > 1)
  
  i <- sample(1:nrow(multistrain), 1)
  cluster_i <- multistrain[i,]
  
  toyset <- checkEachIsolate(cluster_i, tp2$coded, tp2$comps)
  
  df <- oneHeight(toyset)
  
  df %>% pull(actual_size_change) %>% 
    expect_equal(., df$tp2_cl_size - df$tp1_cl_size)
  
  df %>% pull(actual_growth_rate) %>% 
    expect_equal(., round((df$tp2_cl_size - df$tp1_cl_size)/df$tp1_cl_size, digits = 3))
  
  df %>% pull(new_growth) %>% 
    expect_equal(., round(df$tp2_cl_size / (df$tp2_cl_size - df$num_novs), digits = 3))
})

# findingSneakers <- function(novels, q1, q2, matched) {
#   compmatches <- matched[!duplicated(matched$tp1_id),]
#   
#   did_not_chg <- compmatches %>% filter(tp1_cl_size == tp2_cl_size) %>% add_column(add_TP1 = 0)
#   chg <- compmatches %>% filter(tp1_cl_size != tp2_cl_size)
#   
#   if (nrow(chg) > 0) {
#     # identifying the number of additional TP1 strains (sneakers) that show up in the TP2 cluster 
#     # that each TP1 cluster was first tracked to
#     sneakers <- lapply(1:nrow(chg), function(j) {
#       kc <- chg %>% slice(j) # key cluster j
#       tbl1 <- q1 %>% filter(tp1_id == kc$tp1_id) %>% mutate(status = "tp1_cl_size")
#       tbl2 <- q2 %>% filter(tp2_id == kc$tp2_id)
#       # number of novels in the TP2 cluster it kc was tracked to
#       x2 <- tbl2 %>% filter(status == "novs") %>% nrow()
#       # number of sneakers in the TP2 cluster it kc was tracked to
#       x3 <- tbl2 %>% filter(!(isolate %in% tbl1$isolate) & is.na(status)) %>% 
#         mutate(status = "additional_TP1") %>% nrow()
#       
#       c2_tally <- tibble(tp1_cl_size = nrow(tbl1), num_novs = x2, add_TP1 = x3, 
#                          tp2_cl_size = sum(nrow(tbl1), x2, x3))
#       
#       assert(paste0("oneHeight(): Novels check failed for ", kc$tp1_id), kc$num_novs == c2_tally$num_novs)
#       assert(paste0("Composition not calculated properly for ", kc$tp2_id, 
#                     " when tracking ", kc$tp1_id), c2_tally$tp2_cl_size == kc$tp2_cl_size)
#       
#       left_join(kc, c2_tally, by = c("tp1_cl_size", "tp2_cl_size", "num_novs")) %>% return()
#     }) %>% bind_rows()
#     
#     results <- bind_rows(did_not_chg, sneakers)
#   }else {
#     results <- did_not_chg
#   }
#   
#   results %>% return()
# }
test_results[[4]] <- test_that("findingSneakers()", {
  expect_true(FALSE)
})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input CGM tracking functions tests passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}