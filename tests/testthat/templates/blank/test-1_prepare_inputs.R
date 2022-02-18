libs <- c("optparse", "magrittr", "readr", "dplyr", "testit", "data.table", "tibble")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

# source("scripts/Misc/formatprep.R")

test_results <- vector(length = 12) %>% 
  setNames(c("Part 1", "Part 2", "Part 3", "Part 4", "Part 5", "Part 6", 
             "Part 7", "Part 8", "Part 9", "Part 10", "Part 11", "Part 12"))

# FUNCTIONS ----------------------------------------------------------------------------------------------------
checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

# function for reading raw strains and time point clusters
readData <- function(path, check_enc = TRUE) {
  if (check_enc) {
    enc <- checkEncoding(file.path(path))
  }else {
    enc <- ""
  }
  
  file.path(path) %>% 
    read.table(sep="\t", header=TRUE, fileEncoding=enc, fill=TRUE, quote="") %>% 
    as_tibble() %>% return()
}

# START TESTING --------------------------------------------------------------------------------------------
init_inputs <- c("form_inputs.txt", "tp2_clusters_init.txt", "strain_info.txt")
post_inputs <- c("allTP2.Rds", "strain_info.txt")

# Part 1
# required files exist and 1_prepare_inputs.rscript has been run successfully
test_that("Part 1", {
  inps_pre_script1 <- init_inputs %in% list.files("inputs")
  
  expect_true(all(inps_pre_script1))
  expect_true(dir.exists("inputs/processed"))
  
  inps_post_script1 <- post_inputs %in% list.files("inputs/processed/")

  expect_true(all(inps_post_script1))
})


test_that("Part 2", {
  a1 <- readData(file.path("inputs", init_inputs[3]))
  b1 <- c("Source", "Country", "Province", "City",
            "Day", "Month", "Year") %in% colnames(a1)
  b2 <- a1 %>% pull(Day) %>% unique() %>% sort() %>% as.integer()
  b3 <- a1 %>% pull(Month) %>% unique() %in% 1:12 %>% all()

  # required columns found
  expect_true(all(b1))
  # day columns are in [1,31]
  expect_true(all(b2 %in% 1:31))
  # month columns are in [1,12]
  expect_true(all(b3 %in% 1:12))
})

# check that updateStrains was done correctly, processed input files should have 
# the same sets of strains
strain_data <- read_tsv(file.path("inputs/processed", post_inputs[2]))
tp2_data <- readRDS(file.path("inputs/processed", post_inputs[1]))


# Check that dates were processed correctly
test_that("Part 3", {
  col1 <- c("Strain", "Latitude", "Longitude", "Date", "YearMonth", "YearWeek") %in% colnames(strain_data)
  expect_true(all(col1))
  
  x1 <- strain_data %>% pull(Date)
  x2 <- strain_data %>% pull(YearMonth)
  lapply(1:nrow(strain_data), function(i) {
    substr(x1[i], 1, 7) %>% identical(., x2[i])
    }) %>% unlist() %>% all() %>% expect_true()

  x3 <- strain_data %>% pull(YearWeek)
  weeks <- strain_data %>% pull(Date) %>% strftime(., format = "%Y-%V")
  expect_identical(x3, weeks)
})

# check that new_cols, lookup, and original are correct
test_that("Part 4a", {
  expect_identical(c("lookup_table", "new_cols", "original", "strain_pango", "pango_clusters"), names(tp2_data))
  
  lookup <- tp2_data$lookup_table %>% arrange(old_h, old_cl) %>% mutate(across(new_h, as.character))
  
  tp2_raw <- read_tsv(file.path("inputs", init_inputs[2]))
  tp2_raw %>% as.data.table() %>% melt.data.table(id.vars = "Strain") %>% select(-Strain) %>% unique()
  
  part_new <- tp2_data$new_cols %>% as.data.table() %>% melt.data.table(id.vars = "Strain") %>% 
    set_colnames(c("Strain", "new_h", "new_cl")) %>% mutate(across(new_h, as.character)) %>% 
    mutate(across(new_h, as.integer))
  
  part_old <- tp2_data$original %>% as.data.table() %>% melt.data.table(id.vars = "Strain") %>% 
    set_colnames(c("Strain", "old_h", "old_cl"))
  
  old_h_vals <- part_old$old_h %>% unique() %>% sort()
  df <- data.table(old_h = old_h_vals, new_h = 0:(length(old_h_vals)-1))
  
  new_lookup <- part_old %>% 
    merge.data.table(., df, by = "old_h") %>% 
    merge.data.table(., part_new, by = c("Strain", "new_h")) %>% 
    select(old_h, old_cl, new_h, new_cl) %>% unique() %>% 
    arrange(old_h, old_cl) %>% 
    mutate(across(old_h, as.character)) %>% 
    mutate(across(new_h, as.character))
    
  expect_equal(new_lookup, lookup)
})

pango_lineages <- read.csv("inputs/GISAID Lineages (770000 isolates).csv") %>% 
  select(Strain, T0.original, Pango_lineage) %>% 
  as.data.table()

test_that("Part 4d", {
  tp2_data$strain_pango %>% as.data.table() %>% 
    expect_equal(., pango_lineages)
})

test_that("Part 4e", {
  
  df <- tp2_data$new_cols %>% 
    melt.data.table(id.vars = "Strain") %>% 
    set_colnames(c("Strain", "new_h", "new_cl"))
  
  lookup <- tp2_data$lookup_table %>% arrange(old_h, old_cl) %>% mutate(across(new_h, as.character))
  
  actual <- tp2_data$pango_clusters %>% 
    mutate(across(new_h, as.integer))
  
  new_df <- tp2_data$strain_pango %>% 
    merge.data.table(., df, by = "Strain") %>% 
    merge.data.table(., lookup, by = c("new_h", "new_cl")) %>% 
    select(colnames(tp2_data$pango_clusters))
    
  expect_equal(new_df, actual)
    
  
  
})

# IMPORTANT!
# Still need to check column by column of the processed data

# strain_pango
# pango_clusters

# Part 2
test_that("Part 2", {})


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