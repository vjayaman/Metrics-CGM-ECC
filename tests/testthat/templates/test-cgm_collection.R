libs <- c("R6", "tibble", "optparse", "magrittr", "dplyr", "reshape2", "progress", "testit")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

files <- paste0("scripts/CGM") %>% list.files(., full.names = TRUE)
invisible(sapply(files, source))

test_results <- vector(length = 12) %>% 
  setNames(c("Part 1", "Part 2", "Part 3", "Part 4", "Part 5", "Part 6", 
             "Part 7", "Part 8", "Part 9", "Part 10", "Part 11", "Part 12"))
# Part 1
# f1 <- readBaseData(arg$tp1, 1, reader::get.delim(arg$tp1))
# f2 <- readBaseData(arg$tp2, 2, reader::get.delim(arg$tp2))
# colnames(f1)[1] <- colnames(f2)[1] <- "isolate"
# 
# heights <- strsplit(as.character(arg$heights), split = ",") %>% unlist()
# 
# all_isolates <- unique(c(f1$isolate, f2$isolate)) %>% as_tibble() %>% 
#   set_colnames("char_isolate") %>% rowid_to_column("num_isolate")
# 
# ph <- max(nchar(colnames(f1)[-1]), nchar(colnames(f2)[-1]))
# pc <- f1 %>% select(-isolate) %>% max(., f1 %>% select(-isolate)) %>% nchar()

test_that("Part 1", {})

# tp1 <- Timedata$new("tp1", raw = f1, all_isolates, pad_height = ph, pad_cluster = pc)$
#   set_comps()$flag_clusters()
# tp2 <- Timedata$new("tp2", raw = f2, all_isolates, pad_height = ph, pad_cluster = pc)$
#   set_comps()$set_cnames()
test_that("Part 2", {})

# novels <- setdiff(tp2$coded$isolate, tp1$coded$isolate)
# tp2$comps <- tp2$coded %>% filter(isolate %in% novels) %>% 
#   group_by(tp2_id) %>% 
#   summarise(num_novs = n(), .groups = "drop") %>% 
#   left_join(tp2$comps, ., by = "tp2_id") %>% 
#   mutate(num_novs = ifelse(is.na(num_novs), 0, num_novs))
# tp2$flag_clusters()$coded_status(novels)
# tp1$coded_status(novels)
test_that("Part 3", {})

# hx <- Heightdata$new(starter = heights[1], t1_comps = tp1$comps, hvals = heights)$
#   clust_tracking(tp2$comps, tp2$cnames, tp1$coded, tp2$coded, TRUE)$
#   update_iteration()
test_that("Part 4", {})

# if (length(heights) > 1) {
#   outputDetails(paste0("\nStep 3 OF 3: Tracking and flagging clusters for the rest of the heights (",
#                        length(heights) - 1, " of them) ..........."), newcat = TRUE)
#   outputDetails(paste0("  This may take some time. \n  For a more detailed look at progress, ", 
#                        "see the logfile in the logs directory.\n"))
#   outputDetails("  Collecting data for other heights: ", newcat = TRUE)
#   
#   fcb <- txtProgressBar(min = 0, max = length(heights[-1])*2, initial = 0, style = 3)
#   for (j in 1:length(heights[-1])) {
#     hx$h_after <- heights[-1][j]
#     message(paste0("  Height ", j + 1, " / ", length(heights)))
#     
#     # Part 1: unchanged(): identifying clusters that have not changed from the previous height
#     # Part 2: takes comps for new height, previous height, and the tracked data for the previous height
#     #   --> identifies clusters that have not changed since the previous height, and reuses their tracking data
#     hx$post_data(tp1$comps)$unchanged()
#     
#     # Part 2: tracking clusters that changed, saving to results list, and prepping variable for next height
#     hx$comps <- hx$aft %>% filter(!(id_aft %in% hx$same$tp1_id)) %>% set_colnames(colnames(tp1$comps))
#     
#     setTxtProgressBar(fcb, j*2 - 1)
#     
#     hx$clust_tracking(tp2$comps, tp2$cnames, tp1$coded, tp2$coded, FALSE)$
#       update_iteration()$reset_values()
#     
#     setTxtProgressBar(fcb, j*2)
#   }
#   close(fcb)
# }else {
#   outputDetails(paste0("\nStep 3 OF 3: Only one threshold provided, so no further tracking necessary"), newcat = TRUE)
# }
test_that("Part 5", {})

# clusters_just_tp1 <- lapply(heights, function(h) {
#   hx$results[[h]] %>% left_join(., tp1$flagged) %>% 
#     left_join(., tp2$flagged) %>% 
#     arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
#     findingSneakers(novels, tp1$status, tp2$status, .) %>% return()
# }) %>% bind_rows()
test_that("Part 6", {})

# isolates_base <- tp1$melted %>% mutate(across(tp1_h, as.integer)) %>% 
#   filter(tp1_h %in% heights) %>% 
#   left_join(., clusters_just_tp1, by = c("tp1_id", "tp1_h", "tp1_cl")) %>% 
#   arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)
test_that("Part 7", {})

# first_nov_flag <- tp2$status %>% filter(!is.na(status)) %>% group_by(isolate) %>% slice(1) %>% ungroup() %>% pull(tp2_id)
# novels_only_tracking <- tp2$flagged %>% filter(tp2_id %in% first_nov_flag) %>% 
#   left_join(., tp2$comps[,c("tp2_id", "num_novs")]) %>% arrange(tp2_h, tp2_cl) %>% 
#   add_column(tp1_id = NA, tp1_h = NA, tp1_cl = NA, first_tp1_flag = NA, last_tp1_flag = NA)
# novel_asmts <- tp2$melted %>% mutate(across(tp2_h, as.integer)) %>% 
#   filter(isolate %in% setdiff(tp2$raw$isolate, tp1$raw$isolate))
test_that("Part 8", {})

# pure_novels <- novels_only_tracking %>% filter(num_novs == tp2_cl_size) %>% 
#   mutate(tp1_cl_size = 0, add_TP1 = 0) %>% 
#   right_join(novel_asmts, ., by = c("tp2_h", "tp2_cl", "tp2_id")) %>% select(colnames(isolates_base))
# mixed_novels <- novels_only_tracking %>% filter(num_novs != tp2_cl_size) %>% 
#   left_join(., novel_asmts, by = c("tp2_id", "tp2_h", "tp2_cl"))
# all_mixed <- isolates_base %>% 
#   select(-isolate, -tp1_id, -tp1_h, -tp1_cl, -first_tp1_flag, -last_tp1_flag) %>% 
#   unique() %>% left_join(mixed_novels, .) %>% select(colnames(isolates_base))
test_that("Part 9", {})

# if (nrow(pure_novels) > 0) {
#   isolates_file <- bind_rows(isolates_base, pure_novels)  
# }else {
#   isolates_file <- isolates_base
# }
# if (nrow(all_mixed) > 0) {
#   isolates_file %<>% bind_rows(., all_mixed)
# }
# isolates_file %<>% 
#   mutate(novel = ifelse(isolate %in% setdiff(tp2$raw$isolate, tp1$raw$isolate), 1, 0)) %>% 
#   rename(Strain = isolate)
# isolates_file[,c("tp1_h", "tp2_h")] %<>% apply(., 2, padCol, padval = ph, padchr = "h")
# isolates_file[,c("tp1_cl", "tp2_cl")] %<>% apply(., 2, padCol, padval = pc, padchr = "c")
test_that("Part 10", {})

# isolates_file %<>% 
#   mutate(tp1_cl_size = tp1_cl_size + 1, tp2_cl_size = tp2_cl_size + 1) %>% 
#   oneHeight() %>% 
#   addingType(.)
# isolates_file %>% 
#   select(Strain, novel, first_tp2_flag, tp2_h, tp2_cl, tp2_cl_size, last_tp2_flag, 
#          tp1_id, tp1_h, tp1_cl, tp1_cl_size, first_tp1_flag, last_tp1_flag, add_TP1, 
#          num_novs, actual_size_change, actual_growth_rate, new_growth, type)
test_that("Part 11", {})

# write.table(isolates_file, file.path("results","CGM_strain_results.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
test_that("Part 12", {})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input CGM collection tests have passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}