libs <- c("optparse", "magrittr", "readr", "dplyr", "testit")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/formatprep.R")

test_results <- vector(length = 12) %>% 
  setNames(c("Part 1", "Part 2", "Part 3", "Part 4", "Part 5", "Part 6", 
             "Part 7", "Part 8", "Part 9", "Part 10", "Part 11", "Part 12", "Part 13"))


# Part 1
# params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>% 
#   set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
#               "th","nsTP1","nsTP2","temp_win","cnames"))
# a1 <- readData(arg$metadata, FALSE)
# a2 <- suppressWarnings(readData(arg$metadata, check_enc = TRUE))
# if (nrow(a1) > nrow(a2)) {strain_data <- a1}else {strain_data <- a2}
# time1 <- suppressWarnings(readData(arg$tp1, check_enc = TRUE))
# if (!exists("time1")) {time1 <- readData(arg$tp1, FALSE)}
# time2 <- suppressWarnings(readData(arg$tp2, check_enc = TRUE))
# if (!exists("time2")) {time2 <- readData(arg$tp2, FALSE)}
test_that("Part 1", {})

# Part 2
# initial_sizes <- tibble(type="initial", a=nrow(strain_data), b=nrow(time1), d=nrow(time2))
# assert("Has lineage info", as.logical(params$has_lin[2]))
# x <- updateStrains("lin_info", strain_data, time1, time2, initial_sizes)
# strain_data <- x$sd; time1 <- x$t1; time2 <- x$t2; initial_sizes <- x$sizes
test_that("Part 2", {})

# Part 3
# reqnames <- c("Strain", "Latitude", "Longitude", "Day", "Month", "Year")
# cnames <- params$cnames[2] %>% strsplit(split = ",") %>% unlist()
# if (!("none" %in% cnames)) {
#   fullcnames <- c(reqnames, cnames)
# }else {
#   fullcnames <- reqnames
# }
# strain_data <- strain_data %>% select(all_of(fullcnames)) %>% 
#   na.omit(Strain) %>% na.omit(Latitude) %>% na.omit(Longitude) %>% na.omit(Day) %>% 
#   na.omit(Month) %>% na.omit(Year)
test_that("Part 3", {})

# Part 4
# # add column to show which strains are found in TP1
# strain_data %<>% mutate(TP1 = ifelse(Strain %in% time1$Strain, 1, 0))
# # add column to show which strains are found in TP2
# strain_data %<>% mutate(TP2 = ifelse(Strain %in% time2$Strain, 1, 0))
test_that("Part 4", {})

# Part 5
# if (params$reg[2] != "All" & "Region" %in% colnames(strain_data)) {
#   strain_data <- strain_data %>% filter(Region %in% params$reg[2])
# }
test_that("Part 5", {})

# Part 6
# if (params$cou[2] != "All" & "Country" %in% colnames(strain_data)) {
#   strain_data <- strain_data %>% filter(Country %in% params$cou[2])
# }
test_that("Part 6", {})

# Part 7
# if (as.logical(params$has_date[2])) {
#   strain_data <- strain_data %>% 
#     filter(!is.na(Day)) %>% filter(!is.na(Month)) %>% filter(!is.na(Year)) %>% 
#     filter(Day != "") %>% filter(Month != "") %>% filter(Year != "")
# }
test_that("Part 7", {})

# Part 8
# if (nchar(params$temp_win[2]) > nchar("[,]")) {
#   tempwindow <- params$temp_win[2] %>% gsub("\\[|\\]", "", .) %>% 
#     strsplit(., ",") %>% unlist()
#   strain_data %<>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
#     filter(Date >= tempwindow[1] & Date <= tempwindow[2])
#   # Strains with metadata and defined lineage info at TP1
#   y <- updateStrains("temp_win", strain_data, time1, time2, initial_sizes)
#   strain_data <- y$sd; time1 <- y$t1; time2 <- y$t2; initial_sizes <- y$sizes
# }
test_that("Part 8", {})

# Part 9
# if (as.logical(params$has_prov[2])) {
#   strain_data <- strain_data %>% filter(!is.na(Province)) %>% filter(Province != "")
#   if (params$prov[2] != "All") {
#     strain_data <- strain_data %>% filter(Province %in% params$prov[2])
#   }
#   # Strains with metadata and defined lineage info at TP2
#   z <- updateStrains("aft_prov", strain_data, time1, time2, initial_sizes)
#   strain_data <- z$sd; time1 <- z$t1; time2 <- z$t2; initial_sizes <- z$sizes
# }
test_that("Part 9", {})

# Part 10
# processed_tp1 <- intClusters(time1)
# processed_tp2 <- intClusters(time2)
test_that("Part 10", {})

# Part 11
# th <- params$th[2]
# if (as.logical(params$nsTP1[2])) {
#   remove_strains <- strainsInSingletons(processed_tp1$new_cols, th)
#   processed_tp1$new_cols %<>% filter(!(Strain %in% remove_strains))
#   processed_tp2$new_cols %<>% filter(!(Strain %in% remove_strains))
#   strain_data %<>% filter(!(Strain %in% remove_strains))
#   initial_sizes %<>% add_row(tibble(type="tp1_ns", a=nrow(strain_data), b=nrow(time1), d=nrow(time2)))
# }
test_that("Part 11", {})

# Part 12
# if (as.logical(params$nsTP2[2])) {
#   remove_strains <- strainsInSingletons(processed_tp2$new_cols, th)
#   processed_tp1$new_cols %<>% filter(!(Strain %in% remove_strains))
#   processed_tp2$new_cols %<>% filter(!(Strain %in% remove_strains))
#   strain_data %<>% filter(!(Strain %in% remove_strains))
#   initial_sizes %<>% add_row(tibble(type="tp2_ns", a=nrow(strain_data), b=nrow(time1), d=nrow(time2)))
# }
test_that("Part 12", {})

# Part 13
# writeData(strain_data, file.path(arg$inputdir, "processed", "strain_info.txt"))
# writeData(processed_tp1$new_cols, file.path(arg$inputdir, "processed", "tp1_clusters.txt"))
# saveRDS(processed_tp1, file.path(arg$inputdir, "processed", "allTP1.Rds"))
# writeData(processed_tp2$new_cols, file.path(arg$inputdir, "processed", "tp2_clusters.txt"))
# saveRDS(processed_tp2, file.path(arg$inputdir, "processed", "allTP2.Rds"))
test_that("Part 13", {})

test_results %<>% unlist(test_results)
if (all(test_results)) {
  cat("\nAll valid-input prepare_input tests have passed.\n")
}else {
  cat(paste0("\nTests ", paste0(names(test_results)[which(test_results == FALSE)], 
                                collapse = ", "), " failed.\n"))
}