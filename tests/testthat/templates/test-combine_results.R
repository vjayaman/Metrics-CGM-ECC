libs <- c("optparse","magrittr","tibble", "dplyr", "readr", "testit", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/MergeFunctions/type_handling.R")

test_results <- vector(length = 14) %>% 
  setNames(c("Part 1", "Part 2", "Part 3", "Part 4", "Part 5", "Part 6", "Part 7", 
             "Part 8", "Part 9", "Part 10", "Part 11", "Part 12", "Part 13", "Part 14"))

# Part 1
# eccs <- readData(arg$ECCs)
# cgms <- readData(arg$CGMs)
# filtering_params <- readLines(arg$details, warn = FALSE) %>% 
#   strsplit(., split = ": ") %>%
#   set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
#               "th","nsTP1","nsTP2","temp_win","cnames"))
# cnames <- filtering_params$cnames[2] %>% strsplit(split = ",") %>% unlist()
# strain_data <- suppressMessages(read_tsv(arg$strains)) %>% 
#   mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
#   select(Strain, Latitude, Longitude, Day, Month, Year, Date, all_of(cnames), TP1, TP2) %>% 
#   as.data.table() %>% 
#   filter(TP2 == 1)
test_that("Part 1", {})

# Part 2
# step1 <- merge.data.table(cgms, eccs, by = "Strain") %>% 
#   merge.data.table(., strain_data, by = "Strain") %>% 
#   rename(found_in_TP1 = TP1, found_in_TP2 = TP2)
# ecccols <- grep("ECC", colnames(step1), value = TRUE) %>% sort(decreasing = TRUE)
# tp1eccs <- grep("TP1", ecccols, value = TRUE)
# tp2eccs <- grep("TP2", ecccols, value = TRUE)
# assert("No clusters with unassigned type", checkTypes(step1))
test_that("Part 2", {})

# Part 3
# cases2 <- step1 %>% filter(type == "Type2") %>% type2Inheritance(.)
# step1 <- step1[ type != "Type2" ] %>% bind_rows(cases2)
test_that("Part 3", {})

# Part 4
# cases3a <- step1 %>% filter(type == "Type3")
# if (nrow(cases3a) > 0) {
#   cases3b <- type3Inheritance(cases3a)
#   step1 <- step1[ type != "Type3" ] %>% bind_rows(cases3b)
# }
test_that("Part 4", {})

# Part 5
# index_eccs <- grep("ECC", colnames(step1))
# assert("All TP1 and TP2 ECCs for these should be blank (NA, for now)", 
#        all(is.na(step1[ type == "Type4", ..index_eccs])))
# cases4 <- step1 %>% filter(type == "Type4")
# assert("All singletons or nonexistent (at both TP1 and TP2)", 
#        all(c(cases4$tp1_cl_size, cases4$tp2_cl_size) < 3))
# assert("Only blank ECCs are for Type4 cases", identical(unique(step1[is.na(get(ecccols)),type]), "Type4"))
test_that("Part 5", {})

# Part 6
# step2 <- step1 %>% as_tibble()
# for (i in 1:(length(ecccols)/2)) {
#   a <- tp1eccs[i] %>% strsplit(., split = "ECC") %>% unlist() %>% extract2(2) %>% 
#     substr(., 2, nchar(.)) %>% paste0("delta_ECC_", .)
#   z <- pull(step2, tp2eccs[i]) - pull(step2, tp1eccs[i])
#   step2[a] <- z
# }
test_that("Part 6", {})

# Part 7
# step3 <- step2 %>% 
#   left_join(., getAverage(step2, tp1_cl, Latitude, "avg_lat_1"), by = "tp1_cl") %>% 
#   left_join(., getAverage(step2, tp1_cl, Longitude, "avg_long_1"), by = "tp1_cl") %>% 
#   left_join(., getAverage(step2, tp2_cl, Latitude, "avg_lat_2"), by = "tp2_cl") %>% 
#   left_join(., getAverage(step2, tp2_cl, Longitude, "avg_long_2"), by = "tp2_cl")
test_that("Part 7", {})

# Part 8
# step4 <- step3 %>% 
#   left_join(., getAverage(step3, tp1_cl, Date, "avg_date1"), by = "tp1_cl") %>% 
#   left_join(., getAverage(step3, tp2_cl, Date, "avg_date2"), by = "tp2_cl")
test_that("Part 8", {})

# Part 9
# dist_avgs <- grep("avg", colnames(step4), value = TRUE) %>% grep("dist", ., value = TRUE) %>% sort()
# time1 <- readRDS(arg$tp1)$lookup_table %>% 
#   select(1, as.double(filtering_params$th[2])+2) %>% 
#   set_colnames(c("Strain", "Actual TP1 cluster"))
# time2 <- readRDS(arg$tp2)$lookup_table %>% 
#   select(1, as.double(filtering_params$th[2])+2) %>% 
#   set_colnames(c("Strain", "Actual TP2 cluster"))
# step5 <- step4 %>% left_join(., time1) %>% left_join(., time2)
test_that("Part 9", {})

# Part 10
# step6 <- step5 %>% 
#   select(Strain, Country, Province, City, Latitude, Longitude, Day, Month, Year, found_in_TP1, 
#          `Actual TP1 cluster`, tp1_cl, TP1_T0_Size, all_of(tp1eccs), found_in_TP2, 
#          `Actual TP2 cluster`, tp2_cl, TP2_T0_Size, all_of(tp2eccs), 
#          grep("delta_ECC", colnames(step5), value = TRUE), avg_date1, 
#          getDistCols(dist_avgs, "TP1", "temp", TRUE), avg_lat_1, avg_long_1, 
#          getDistCols(dist_avgs, "TP1", "geog", TRUE), avg_date2, 
#          getDistCols(dist_avgs, "TP2", "temp", TRUE), avg_lat_2, avg_long_2, 
#          getDistCols(dist_avgs, "TP2", "geog", TRUE), first_tp1_flag, last_tp1_flag, 
#          first_tp2_flag, last_tp2_flag, tp1_cl_size, tp2_cl_size, 
#          actual_size_change, add_TP1, num_novs, actual_growth_rate, new_growth, type)
test_that("Part 10", {})

# Part 11
# step7 <- step6 %>% 
#   rename_with(., replaceDistName, getDistCols(colnames(.), "TP1", "temp")) %>% 
#   rename_with(., replaceDistName, getDistCols(colnames(.), "TP2", "temp")) %>% 
#   rename_with(., replaceDistName, getDistCols(colnames(.), "TP1", "geo")) %>% 
#   rename_with(., replaceDistName, getDistCols(colnames(.), "TP2", "geo"))
test_that("Part 11", {})

# Part 12
# step8 <- step7 %>% 
#   rename("TP1 cluster" = tp1_cl, "TP1 cluster size (1)" = TP1_T0_Size, 
#          "TP2 cluster" = tp2_cl, "TP2 cluster size (1)" = TP2_T0_Size, 
#          "TP1" = found_in_TP1, "TP2" = found_in_TP2, 
#          "Average TP1 date" = avg_date1, "Average TP1 latitude"	= avg_lat_1, 
#          "Average TP1 longitude" = avg_long_1, 
#          "Average TP2 date" = avg_date2, "Average TP2 latitude" = avg_lat_2, 
#          "Average TP2 longitude" = avg_long_2, 
#          "First time this cluster was seen in TP1" = first_tp1_flag, 
#          "Last time this cluster was seen in TP1" = last_tp1_flag, 
#          "First time this cluster was seen in TP2" = first_tp2_flag, 
#          "Last time this cluster was seen in TP2" = last_tp2_flag, 
#          "TP1 cluster size + 1 (2)"	= tp1_cl_size, "TP2 cluster size + 1 (2)"	= tp2_cl_size, 
#          "Actual cluster size (TP2 size - TP1 size)" = actual_size_change, 
#          "Number of additional TP1 strains in the TP2 match" = add_TP1, 
#          "Number of novels in the TP2 match" = num_novs, 
#          "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)" = actual_growth_rate, 
#          "Novel growth = (TP2 size) / (TP2 size - number of novels)" = new_growth, 
#          "Type" = type) %>% 
#   arrange(`TP2 cluster`, `TP1 cluster`, Strain)
test_that("Part 12", {})

# Part 13
# writeData(fp = "results/Merged_strain_results.tsv", df = step8)
test_that("Part 13", {})

# Part 14
# step8 %>% 
#   group_by(`TP2 cluster`) %>% slice(1) %>% 
#   select(-Strain) %>% ungroup() %>% 
#   writeData(fp = "results/Merged_cluster_results.tsv", df = .)
test_that("Part 14", {})


test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input combining results tests have passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}