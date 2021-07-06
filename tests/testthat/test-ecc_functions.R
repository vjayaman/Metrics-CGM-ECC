libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
library(testthat)
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/ECC/functions/ecc_functions.R")
source("scripts/ECC/functions/dist_functions.R")

fakeRaw <- function(ndrs, n, nc) {
  data.table(th = sample(1:nc, n, replace = TRUE) %>% sort(), 
             dr = rep(1:ndrs, ceiling(n/ndrs))[1:n]) %>% 
    rownames_to_column("Strain") %>% return()
}

fakeSecCluster <- function(ndrs, n, nc) {
  testcase <- fakeRaw(ndrs, n, nc)
  nunique <- length(unique(testcase$dr))
  
  assignments <- data.table(
    dr = testcase$dr %>% unique(), 
    Date = sample(seq(as.Date('2019/01/01'), as.Date('2021/01/01'), by="day"), nunique), 
    Latitude = sample(seq(25.00000, 65.00000, by = 0.00001), nunique), 
    Longitude = sample(seq(-100.00000, 125.00000, by = 0.00001), nunique))
  
  # fake sectionClusters() result object:
  parts <- formatForSectioning(testcase, 5) %>% sectionTypingData() %>% 
    list("drs" = testcase, "results" = .)
  return(list("a" = assignments, "b" = parts, "c" = testcase))
}

fakeStrains <- function(nrows, cnames = NULL) {
  sample_dates <- seq(as.Date('2019/01/01'), as.Date('2021/01/01'), by="day") %>% 
    sample(., nrows) %>% as.character() %>% strsplit(split = "-")
  
  basic_strains <- tibble(
    Strain = paste0("Strain.", 1:nrows), 
    Latitude = sample(seq(25.00000, 65.00000, by = 0.00001), size = nrows, replace = TRUE),
    Longitude = sample(seq(-100.00000, 125.00000, by = 0.00001), size = nrows, replace = TRUE), 
    Day = sapply(sample_dates, "[[", 3) %>% as.double(),
    Month = sapply(sample_dates, "[[", 2) %>% as.double(),
    Year = sapply(sample_dates, "[[", 1) %>% as.double()
  )
  
  if ("Country" %in% cnames) {
    basic_strains %<>% add_column(Country = paste0("Country.", 1:nrow(basic_strains)))
  }
  
  if ("Province" %in% cnames) {
    basic_strains %<>% add_column(Province = paste0("Province.", 1:nrow(basic_strains)))
  }
  
  if ("City" %in% cnames) {
    basic_strains %<>% add_column(City = paste0("City.", 1:nrow(basic_strains)))
  }
  
  return(basic_strains)
}

# (3) timeTaken() function
# - Indicates length of a process in hours, minutes, and seconds, when given a name of the process 
# ("pt") and a two-element named vector with Sys.time() values named "start_time" and "end_time"

test_that("timeTaken() - null inputs", {
  stopwatch <- list("start_time" = NULL, "end_time" = NULL)
  expect_identical(timeTaken("test", stopwatch), "Neither start nor end time were collected")
  
  stopwatch$end_time <- as.character.POSIXt(Sys.time())
  expect_identical(timeTaken("test", stopwatch), "Start time was not collected.")
  
  stopwatch$start_time <- as.character.POSIXt(Sys.time())
  stopwatch$end_time <- NULL
  expect_identical(timeTaken("test", stopwatch), "End time was not collected.")
})

test_that("timeTaken() - time inputs", {
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


# (4) outputMessages() function
test_that("outputMessages()", {
  expect_identical(outputMessages(""), cat(paste0("\n","")))
  expect_output(outputMessages("Test"), "\nTest")
  expect_silent(outputMessages())
})

# FAKES - used in the following tests
# In this fake object, there are 20 drs, 100 strains, and 15 clusters
fake_obj <- fakeSecCluster(sample.int(50,1), sample.int(200,1)*2, sample.int(20,1))
assignments <- fake_obj$a
results <- fake_obj$b$results
inds <- fake_obj$b$drs[["th"]] %in% pull(results[[1]], "th")
cluster_x <- fake_obj$b$drs[inds,-"Strain"]
cluster_asmts <- assignments[dr %in% pull(cluster_x, dr)]


test_that("transformData2() - temp", {
  dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
  
  returned <- transformData2(dm_temp, "temp", min(dm_temp), max(dm_temp))
  
  logtr <- dm_temp %>% add(10) %>% log10()
  logtr[logtr == -Inf] <- 0
  min_dm <- min(dm_temp) %>% add(10) %>% log10()
  max_dm <- max(dm_temp) %>% add(10) %>% log10()
  
  pairs <- expand.grid(cluster_asmts$dr, cluster_asmts$dr)
  
  for (i in 1:nrow(pairs)) {
    ri <- as.character(pairs$Var1[i])
    ci <- as.character(pairs$Var2[i])
    
    x1 <- (logtr[ri, ci] - min_dm) / (max_dm - min_dm)
    expect_equal(returned[ri, ci], x1)  
  }
})

test_that("transformData2() - geo", {
  dm_geo <- cluster_asmts %>% select(dr, Longitude, Latitude) %>% 
    distMatrix(., "geo", c("Longitude", "Latitude"))
  
  returned <- transformData2(dm_geo, "geo", min(dm_geo), max(dm_geo))
  
  logtr <- dm_geo %>% add(10) %>% log10()
  min_dm <- min(dm_geo) %>% add(10) %>% log10()
  max_dm <- max(dm_geo) %>% add(10) %>% log10()
  
  pairs <- expand.grid(cluster_asmts$dr, cluster_asmts$dr)
  
  if (max_dm == 1) {
    expect_true(all(unique(returned) == 0))
  }else {
    for (i in 1:nrow(pairs)) {
      ri <- as.character(pairs$Var1[i])
      ci <- as.character(pairs$Var2[i])
      
      x1 <- (logtr[ri, ci] - min_dm) / (max_dm - min_dm)
      expect_equal(returned[ri, ci], x1)  
    }  
  }
})

# checkEncoding - can check by inspection, don't necessarily need a test for this?
#   - unless the inputs are invalid, in which case will do this one in the next
#   iteration of test development (malicious/incorrect inputs)
# test_that("checkEncoding()", {})


# formatData - can check by inspection, don't necessarily need a test for this?
#   - unless the inputs are invalid, in which case will do this one in the next 
#   iteration of test development (malicious/incorrect inputs)
# test_that("formatData()", {formatData(dm, c("dr1","dr2","Temp.Dist"))})

test_that("processedStrains() - requirements test", {
  base_strains <- sample.int(200, 1) %>% fakeStrains()
  returned <- processedStrains(base_strains)
  
  strain_data <- base_strains %>% 
    mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
    select(colnames(returned$strain_data))
  
  assignments <- strain_data %>% select(Date, Latitude, Longitude) %>% 
    unique() %>% rownames_to_column("dr") %>% as.data.table()
  
  dr_matches <- left_join(strain_data, assignments, by = c("Latitude", "Longitude", "Date")) %>% 
    select(Strain, dr)
  
  expect_identical(strain_data, returned$strain_data)
  expect_identical(assignments, returned$assignments)
  expect_identical(dr_matches, returned$dr_matches)
})

test_that("processedStrains() - additional columns test", {
  additional_cols <- sample(c("Country","Province","City"), sample.int(3, 1)) %>% sort()
  base_strains <- sample.int(200, 1) %>% fakeStrains(., additional_cols)
  returned <- processedStrains(base_strains)
  
  strain_data <- base_strains %>% 
    mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")), 
           Location = do.call(paste, c(base_strains[additional_cols], sep = "_"))) %>% 
    select(colnames(returned$strain_data))
  
  assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
    unique() %>% rownames_to_column("dr") %>% as.data.table()
  
  dr_matches <- left_join(strain_data, assignments, by = c("Latitude", "Longitude", "Date")) %>% 
    select(Strain, dr)
  
  expect_identical(strain_data, returned$strain_data)
  expect_identical(assignments, returned$assignments)
  expect_identical(dr_matches, returned$dr_matches)
})

# # epi_cohesion_new(g_cuts, epi_melt)
# test_that("epi_cohesion_new()", {
#   fake_obj <- fakeSecCluster(12, 60, 10)
#   dr_matches <- fake_obj$c
#   dr_td1 <- dr_matches %>% select(-Strain)
#   g_cuts <- countDataReps(dr_td1)
# 
#   dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
#   dm_geo <- cluster_asmts %>% select(dr, Longitude, Latitude) %>%
#     distMatrix(., "geo", c("Longitude", "Latitude"))
#   
#   transformed_temp <- transformData2(dm_temp, "temp", min(dm_temp), max(dm_temp)) %>%
#     formatData(., c("dr1","dr2","Temp.Dist"))
#   transformed_geo <- transformData2(dm_geo, "geo", min(dm_geo), max(dm_geo)) %>%
#     formatData(., c("dr1","dr2","Geog.Dist"))
#   transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
# 
#   tau <- 1
#   gamma <- 0
#   epi_table <- transformed_dists %>%
#     mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>%
#     select(dr1, dr2, Total.Dist) %>% as.data.table()
# 
#   epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
#   epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)])
#   rownames(epi_matrix) <- colnames(epi_matrix)
# 
#   epi_melt <- as.matrix(1-epi_matrix) %>%
#     as.data.table(., keep.rownames = TRUE) %>%
#     melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>%
#     set_colnames(c("dr_1", "dr_2", "value")) %>% as.data.table()
# 
#   cx <- colnames(cluster_x)[1]
#   tallied_reps <- cluster_x %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
#   cnames <- intersect(colnames(tallied_reps), colnames(cluster_x))
#   g_cuts <- left_join(cluster_x, tallied_reps, by = cnames) %>%
#     unique() %>% mutate(across(dr, as.character))
# 
#   tpx <- 1
#   td_i <- epi_cohesion_new(g_cuts, epi_melt) %>%
#     set_colnames(c(paste0("TP", tpx, "_", colnames(.))))
#   colnames(td_i) %<>% gsub("ECC", paste0("ECC.0.", tau, ".", gamma), x = .)
# })
# 
# ndrs <- 10
# alldrs <- sample(seq(1,5500), ndrs)
# epi_melt <- data.table(
#   dr_1 = alldrs, 
#   dr_2 = sort(alldrs), 
#   value = sample(seq(0,0.95,0.0001), ndrs)
# )
# 
# nstrains <- ceiling(sample(25:75, 1))
# g_cuts <- data.table(
#   th = sample(1:10, nstrains, replace = TRUE), 
#   dr = sample(alldrs, nstrains, replace = TRUE)
# )
# -----------------------------------------------------------------------------------
# epi_cohesion_new ------------------------------------------------------------------
# -----------------------------------------------------------------------------------
# 
# epi_cohesion_new(g_cuts, epi_melt)
# # Example:
# # epi_melt: 
# #     dr_1 (chr)    dr_2 (chr)      value (dbl)
# #      1              1              0.8957695
# #      10             1              0.4255450
# #      100            1              0.4081672
# #      101            1              0.4081672
# #      105            1              0.4081672
# #      ...           ...                ...
# # min:  1             1              0.0005198079
# # max: 5499          5499            0.9019662
# 
# # g_cuts
# #       T0 (int)      dr (chr)          n (int)
# #         1              1                 1
# #         1              2                 1
# #         1              3                 5
# #         1              4                 1
# #         2              5                 1
# #         3              5                 2
# #         1              5                 3
# #         1              6                 1
# #         4              7                 1
# #         1              7                 1
# #        ...            ...               ...
# # min:    1              1                 1
# # max:    10            5499              154


# epiCollectionByCluster(strain_data, tau, gamma, transformed_dists, tpx, cluster_x)
# collectECCs(k, m, parts, extremes, c1, fpath1, fpath2)


# -----------------------------------------------------------------------------------
# epiCollectionByCluster ------------------------------------------------------------
# -----------------------------------------------------------------------------------

# # ### Incorporating the allele data with the epidemiological data 
# # epiCollectionByCluster(strain_data, tau, gamma, transformed_dists, tpx, cluster_x)
# k <- 1
# m <- read_tsv(params$strains) %>% processedStrains() #list(strain_data, assignments, dr_matches)
# parts <- sectionClusters(k, typing_data, m)
# extremes <- readRDS("results/dist_extremes.Rds")
# c1 <- "001" %>% strsplit(., split = "") %>% unlist() %>% as.numeric()
# fpath1 <- paste0("results/TP", k, "/dists/")
# fpath2 <- paste0("results/TP", k, "/eccs/", "001", "/")
# 
# df <- parts$drs
# cx <- setdiff(colnames(df), c("Strain", "dr"))
# results <- parts$results
# 
# # epiCollectionByCluster(m$strain_data, tau, gamma, transformed_dists, k, cluster_x)
# strain_data <- m$strain_data
# tau <- c1[2]
# gamma <- c1[3]
# 
# j <- 1
# cluster_x <- df[df[[cx]] %in% pull(results[[j]], cx),-"Strain"]
# 
# dms <- paste0(fpath1, "group", j, ".Rds") %>% readRDS()
# transformed_temp <- dms$temp %>% 
#   transformData2(., "temp", extremes$temp$min, extremes$temp$max) %>% 
#   formatData(., c("dr1","dr2","Temp.Dist"))
# transformed_geo <- dms$geo %>%
#   transformData2(., "geo", extremes$geo$min, extremes$geo$max) %>%
#   formatData(., c("dr1","dr2","Geog.Dist"))
# transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
# 
# tpx <- k
# 
# epiCollectionByCluster(strain_data, tau, gamma, transformed_dists, k, cluster_x)
# 
# ### Incorporating the allele data with the epidemiological data 
# epiCollectionByCluster <- function(strain_data, tau, gamma, transformed_dists, tpx, cluster_x) {
  # cat(paste0("\n   Collecting ECC values for temporal = ", tau, ", geo = ", gamma))
  # 
  # outputMessages(paste0("      Preparing table of distances: sqrt(", tau, "*(temp^2) + ", gamma, "*(geo^2))"))
  # epi_table <- transformed_dists %>%
  #   mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>%
  #   select(dr1, dr2, Total.Dist) %>% as.data.table()
  # invisible(gc())
  # 
  # outputMessages("      Generating epi distance matrix ...")
  # epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
  # epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)])
  # rownames(epi_matrix) <- colnames(epi_matrix)
  # rm(epi_table)
  # invisible(gc())
  # 
  # outputMessages("      Formatting into table of similarity values from epi distance matrix ...")
  # epi_melt <- as.matrix(1-epi_matrix) %>%
  #   as.data.table(., keep.rownames = TRUE) %>%
  #   melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>%
  #   set_colnames(c("dr_1", "dr_2", "value")) %>% as.data.table()
  # 
  # rm(epi_matrix)
  # invisible(gc())
  # 
  # outputMessages("      Incorporating the metadata with the clusters (typing data) ...")
  # # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
  # cx <- colnames(cluster_x)[1]
  # tallied_reps <- cluster_x %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
  # cnames <- intersect(colnames(tallied_reps), colnames(cluster_x))
  # g_cuts <- left_join(cluster_x, tallied_reps, by = cnames) %>%
  #   unique() %>% mutate(across(dr, as.character))
  # 
  # td_i <- epi_cohesion_new(g_cuts, epi_melt) %>%
  #   set_colnames(c(paste0("TP", tpx, "_", colnames(.))))
  # colnames(td_i) %<>% gsub("ECC", paste0("ECC.0.", tau, ".", gamma), x = .)
  # 
  # return(td_i)
# }

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------


