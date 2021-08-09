#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

# msg <- file("logs/logfile_datacollection.txt", open="wt")
# sink(msg, type="message")

libs <- c("R6", "tibble", "optparse", "magrittr", "dplyr", "reshape2", "progress", 
          "testit", "data.table", "readr")
y <- lapply(libs, require, character.only = TRUE)

files <- paste0("scripts/CGM") %>% list.files(., full.names = TRUE)
invisible(sapply(files, source))

# READING IN THE INPUTS ----------------------------------------------------------------------------------------
# Change the default values to read in your own files, or feed through terminal arguments
option_list <- list(
  make_option(c("-m", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain metadata file"), 
  make_option(c("-a", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "Time point 1 file name (TP1)"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "Time point 2 file name (TP2)"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = paste0("A string of comma-delimited numbers, e.g. '50,75,100' to ", 
                            "use as heights for metric generation (strain table outputs)")))

arg <- parse_args(OptionParser(option_list=option_list))

# BASIC STARTUP MESSAGES ---------------------------------------------------------------------------------------
outputDetails(paste0("\n||", paste0(rep("-", 32), collapse = ""), " Cluster metric generation ", 
                     paste0(rep("-", 32), collapse = ""), "||\nStarted process at: ", Sys.time()))
cat(paste0("\nIf at any point the process cuts off with no success message, please see the log file.\n"))
outputDetails("\nStep 1 OF 3: Data processing ", newcat = TRUE)

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

# TP DATA PREPARATION ------------------------------------------------------------------------------------------
metadata <- suppressMessages(read_tsv(arg$strains)) %>% 
  mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
  mutate(Week = strftime(Date, format = "%V")) %>% select(-TP1, -TP2) %>% 
  mutate(num_week = as.integer(Week)) %>% 
  arrange(Week) %>% as.data.table()

f2 <- readBaseData(arg$tp2, 2, reader::get.delim(arg$tp2)) %>% as.data.table()

week_list <- sort(unique(metadata$num_week))

for (n2 in week_list[length(week_list):2]) {
  n1 <- n2 - 1
  
  tpx2 <- f2 %>% filter(Strain %in% metadata[num_week <= n2, Strain])
  tpx1 <- tpx2 %>% filter(Strain %in% metadata[num_week <= n1, Strain])
  
  print(paste0(n2, ": ", nrow(tpx2)))
  print(paste0(n1, ": ", nrow(tpx1)))
  
  colnames(tpx1)[1] <- colnames(tpx2)[1] <- "isolate"
  heights <- strsplit(as.character(arg$heights), split = ",") %>% unlist()
  ph <- max(nchar(colnames(tpx1)[-1]), nchar(colnames(tpx2)[-1]))
  pc <- tpx2 %>% select(-isolate) %>% max(., tpx2 %>% select(-isolate)) %>% nchar()
  
  tplist <- tpDataSetup(tpx1, tpx2, ph, pc)
  tp1 <- tplist[["tp1"]]
  tp2 <- tplist[["tp2"]]
  novels <- tplist[["novs"]]
 
  # BASE CASE (FIRST HEIGHT) -------------------------------------------------------------------------------------
  outputDetails("\nStep 2 OF 3: Tracking and flagging clusters for base case ", newcat = TRUE)
  outputDetails(paste0("  Collecting height data for base case, height ", heights[1], "..."), newcat = TRUE)
  
  hx <- Heightdata$new(starter = heights[1], t1_comps = tp1$comps, hvals = heights)$
    clust_tracking(tp2$comps, tp2$cnames, tp1$coded, tp2$coded, TRUE)$
    update_iteration()
  
  outputDetails("  Identifying and counting 'additional TP1 strains'.\n", newcat = FALSE)
  
  clusters_just_tp1 <- lapply(heights, function(h) {
    hx$results[[h]] %>% left_join(., tp1$flagged) %>% 
      left_join(., tp2$flagged) %>% 
      arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
      findingSneakers(novels, tp1$status, tp2$status, .) %>% return()
  }) %>% bind_rows()
  
  outputDetails("  Handling novel tracking, adding to dataset.\n", newcat = FALSE)
  
  isolates_file <- novelHandling(tp1, tp2, clusters_just_tp1)
  
  isolates_file %<>% 
    mutate(novel = ifelse(isolate %in% setdiff(tp2$raw$isolate, tp1$raw$isolate), 1, 0)) %>% 
    rename(Strain = isolate)
  
  isolates_file[,c("tp1_h", "tp2_h")] %<>% apply(., 2, padCol, padval = ph, padchr = "h")
  isolates_file[,c("tp1_cl", "tp2_cl")] %<>% apply(., 2, padCol, padval = pc, padchr = "c")
  
  outputDetails("  Incrementing all cluster sizes by 1, then calculating growth columns.\n", newcat = FALSE)
  outputDetails("  Also adding 'type' column to CGM results table.\n", newcat = FALSE)
  
  isolates_file %<>% 
    mutate(tp1_cl_size = tp1_cl_size + 1, tp2_cl_size = tp2_cl_size + 1) %>% 
    oneHeight() %>% 
    addingType(.)
  
  outputDetails("  Saving the data in a file with cluster identifiers.\n", newcat = FALSE)
  
  isolates_file %>% 
    select(tp1_id, tp1_cl_size, first_tp1_flag, last_tp1_flag, 
           first_tp2_flag, tp2_cl_size, last_tp2_flag, add_TP1, novel, 
           num_novs, actual_size_change, actual_growth_rate, new_growth, type) %>% 
    unique() %>% as.data.table() %>% 
    addInterval(., n1, n2) %>% 
    saveRDS(., paste0("intermediate_data/cgms/interval-", 
                      unique(metadata[num_week==n1,Week]), "-", 
                      unique(metadata[num_week==n2,Week]), ".Rds"))
}

# WRAPPING THINGS UP -------------------------------------------------------------------------------------------
stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())

outputDetails(paste0("\nSuccessfully collected data for all heights."), newcat = TRUE)
timeTaken(pt = "CGM data collection", stopwatch) %>% outputDetails(., newcat = TRUE)
outputDetails(paste0("||", paste0(rep("-", 28), collapse = ""), " End of cluster metric generation ", 
                     paste0(rep("-", 29), collapse = ""), "||\n"))
closeAllConnections()
