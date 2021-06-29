#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_distmatrices.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

files <- paste0("scripts/ECC") %>% 
  list.files(., full.names = TRUE) %>% 
  grep("dist_matrices.R", ., invert = TRUE, value = TRUE)
invisible(sapply(files, source))

# Title: "EpiQuant - Salmonella Enteritidis Project (2019-2020)"
# Authors of original work and initial modifications: Ben Hetman, Elissa Giang, Dillon Barker
# Responsible for changes made during and after merge with CGM process: Vasena Jayamanna

option_list <- list(
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ", 
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

m <- read_tsv(params$strains) %>% processedStrains()

# Note: dr stands for data representative
# in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
# cat(paste0("\n\nStep 1:"))
# cat(paste0("\n   Note that the source coefficent is always 0 in this version"))
# outputMessages("   Removing redundancy (comparing date, lat-long, etc. - not every pair of strains)")
# outputMessages("   Identifying which strains match which non-redundant 'data representatives'")

# COLLECT dist matrices using TP2 clusters ---------------------------------------------------------

# df has colnames Strain, T0, dr (in example):
sectionTypingData <- function(df, maxsize) {
  cx <- setdiff(colnames(df), c("Strain", "dr"))
  drsizes <- df %>% select(-Strain) %>% unique() %>% group_by(!!as.symbol(cx)) %>% summarise(n = n())
  
  m <- max(drsizes$n, maxsize)
  rows_x <- c()
  total <- 0
  results <- vector("list", nrow(drsizes))
  
  for (i in 1:nrow(drsizes)) {
    row_i <- drsizes[i,]
    updated <- row_i$n + total
    
    if (updated < m) {
      rows_x <- bind_rows(rows_x, row_i)
      total <- total + row_i$n
      
    }else if (updated == m) {
      rows_x <- bind_rows(rows_x, row_i)
      total <- total + row_i$n
      
      results[[i]] <- rows_x
      rows_x <- c()
      total <- 0
      
    }else { # updated > m
      results[[i]] <- rows_x
      rows_x <- row_i
      total <- row_i$n
    }
  }
  
  if (total != 0) {results[[i]] <- rows_x}
  
  results[sapply(results, is.null)] <- NULL
  
  return(results)
}

sectionClusters <- function(k, typing_data, m) {
  df <- typing_data[[k]] %>% rownames_to_column("Strain") %>%
    as.data.table() %>% left_join(., m$dr_matches, by = "Strain")
  
  gc()
  
  results <- sectionTypingData(df, 1000)
  assert("No clusters overlooked", length(setdiff(pull(df,2), pull(rbindlist(results),1))) == 0)
  
  return(list("drs" = df, "results" = results))
}

distMatrix <- function(input_data, dtype, cnames) {
  dm_names <- pull(input_data, 1) # dr, or Strain
  if (dtype == "temp") {
    dm <- input_data %>% select(all_of(cnames)) %>% pull() %>% 
      dist(diag = FALSE, upper = FALSE, method = "euclidean")
    dm %>% as.matrix(nrow = nrow(input_data), ncol = nrow(input_data)) %>% 
      set_rownames(dm_names) %>% 
      set_colnames(dm_names) %>% return()
    
  }else if (dtype == "geo") {
    # consider using geosphere::distm() for this
    dm <- input_data %>% select(all_of(cnames)) %>% as.data.frame() %>% earth.dist(dist = TRUE)
    dm %>% as.matrix() %>% 
      set_rownames(dm_names) %>% 
      set_colnames(dm_names) %>% return()
  }
}

collectDistances <- function(k, m, parts) {
  
  df <- parts$drs
  cx <- setdiff(colnames(df), c("Strain", "dr"))
  results <- parts$results
  p <- length(results)
  
  min_geo <- min_temp <- Inf
  max_geo <- max_temp <- -Inf
  
  for (j in 1:p) {
    outputMessages(paste0("Working through group ", j, " / ", p))
    cluster_x <- df[df[[cx]] %in% pull(results[[j]], cx),-"Strain"]
    cluster_asmts <- m$assignments[dr %in% pull(cluster_x, dr)]
    
    outputMessages("   Generating all possible date pair distances ...")
    dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
    
    outputMessages("   Generating all possible lat-long pair distances ...")
    dm_geo <- cluster_asmts %>% select(dr, Latitude, Longitude) %>% 
      distMatrix(., "geo", c("Latitude", "Longitude"))
    
    min_temp <- min(min_temp, min(dm_temp))
    max_temp <- max(max_temp, max(dm_temp))
    
    min_geo <- min(min_geo, min(dm_geo + 10))
    max_geo <- max(max_geo, max(dm_geo))
    
    list(temp = dm_temp, geo = dm_geo) %>% 
      saveRDS(., paste0("results/tmp/TP", k, "-dists/group", j, ".Rds"))
    
    rm(dm_temp)
    rm(dm_geo)
    gc()
  }
  
  if (k == 2) {
    tibble(temp = list("max" = max_temp, "min" = min_temp), 
           geo = list("max" = max_geo, "min" = min_geo)) %>% 
      saveRDS(., paste0("results/tmp/TP", k, "-dists/extremes.Rds"))
  }
}

for (k in 1:2) {
  outputMessages(paste0("\nCollecting and saving distances for groups of clusters at TP", k))
  dir.create(paste0("results/dists/TP", k, "-dists"), showWarnings = FALSE)
  sectionClusters(k, typing_data, m) %>% collectDistances(k, m, .)
}

outputMessages("\nFinished saving distance matrices.")
