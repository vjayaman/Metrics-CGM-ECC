#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/cgm_collection.txt", open="wt")
sink(msg, type="message")

libs <- c("R6", "tibble", "optparse", "magrittr", "dplyr", "reshape2", "progress", 
          "testit", "data.table", "readr")
y <- lapply(libs, require, character.only = TRUE); rm(libs); rm(y)

# READING IN THE INPUTS ----------------------------------------------------------------------------------------
# Change the default values to read in your own files, or feed through terminal arguments
option_list <- list(
  make_option(c("-f", "--intervalfile"), metavar = "file", default = "inputs/processed/clustersets.Rds"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"), 
  make_option(c("-n", "--tpn"), metavar = "file", 
              default = "inputs/processed/allTP2.Rds", help = "TPN data"))

arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)

files <- paste0("scripts/CGM") %>% list.files(., full.names = TRUE)
invisible(sapply(files, source)); rm(files)

# BASIC STARTUP MESSAGES ---------------------------------------------------------------------------------------
outputDetails(paste0("\n||", paste0(rep("-", 32), collapse = ""), " Cluster metric generation ", 
                     paste0(rep("-", 32), collapse = ""), "||\nStarted process at: ", Sys.time()))
cat(paste0("\nIf at any point the process cuts off with no success message, please see the log file.\n\n"))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

# TP DATA PREPARATION ------------------------------------------------------------------------------------------
params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

heights <- strsplit(as.character(params$th[2]), split = ",") %>% unlist()
clustersets <- readRDS(arg$intervalfile)
interval_list <- names(clustersets)

if (params$int_type[2] == "multiset") {
  interval <- "Multiset"
}else if (params$int_type[2] == "monthly") {
  interval <- "YearMonth"
}else if (params$int_type[2] == "weekly") {
  interval <- "Week"
}

save_to <- file.path("intermediate_data", tolower(params$int_type[2]), "cgms")
tpn <- readRDS(arg$tpn)$new_cols

# rowx <- readRDS("results/rowx.Rds")
for (i in 1:(length(interval_list)-1)) {
  
  n1 <- as.character(interval_list[i])
  tpx1a <- clustersets[[n1]]$sofar %>% select(-ivl,-Date) %>% set_colnames(c("isolate", heights))
  tpx1 <- tpn %>% rename("isolate" = "Strain") %>% 
    left_join(tpx1a, ., by = intersect(colnames(tpx1a), colnames(.)))
  
  n2 <- as.character(interval_list[i+1])
  tpx2a <- clustersets[[n2]]$sofar %>% select(-ivl,-Date) %>% set_colnames(c("isolate", heights))
  tpx2 <- tpn %>% rename("isolate" = "Strain") %>% 
    left_join(tpx2a, ., by = intersect(colnames(tpx2a), colnames(.)))
  
  # if (i > 1) {
  #   fullset <- clustersets[[n1]]$sofar
  #   ivl_i <- clustersets[[n1]]$ivl
  #   unchanged_clusters <- setdiff(fullset, ivl_i) %>% pull(heightx) %>% unique()
  #   strains <- fullset[heightx %in% unchanged_clusters] %>% pull(isolate)
  #   rm(fullset); rm(unchanged_clusters); rm(ivl_i)
  #   unchanged_data <- tmp %>% filter(Strain %in% strains)
  # }
  
  ph <- max(nchar(colnames(tpx1)[-1]), nchar(colnames(tpx2)[-1]))
  pc <- tpx2 %>% select(-isolate) %>% max(., tpx2 %>% select(-isolate)) %>% nchar()
  
  msgtexts <- c(
    paste0("  Constructing ", n1, " data object, ", nrow(tpx1), " (", i, " / ", length(interval_list), "):\n"), 
    paste0("  Constructing ", n2, " data object, ", nrow(tpx2), " (", i+1, " / ", length(interval_list), "):\n")
  )
  
  tplist <- tpDataSetup(tpx1, tpx2, ph, pc, FALSE, msgtexts)#; rm(tpx1); rm(tpx2)
  tp1 <- tplist[["tp1"]]
  tp2 <- tplist[["tp2"]]
  novels <- tplist[["novs"]]
  rm(tplist)
  
  # BASE CASE (FIRST HEIGHT) -------------------------------------------------------------------------------------
  outputDetails(paste0("  Tracking clusters from ", n1, " to ", n2, ", from height ", heights[1], " ..."), newcat = TRUE)
  
  hx <- Heightdata$new(starter = heights[1], t1_comps = tp1$comps, hvals = heights)$
    clust_tracking(tp2$comps, tp2$cnames, tp1$coded, tp2$coded, TRUE)$
    update_iteration()
  
  # outputDetails("  Identifying and counting 'additional TP1 strains'.\n", newcat = FALSE)
  clusters_just_tp1 <- lapply(heights, function(h) {
    df1 <- hx$results[[h]]
    df2 <- left_join(df1, tp1$flagged, by = intersect(colnames(df1), colnames(tp1$flagged)))
    
    left_join(df2, tp2$flagged, by = intersect(colnames(df2), colnames(tp2$flagged))) %>% 
      arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
      findingSneakers(novels, tp1$status, tp2$status, .) %>% return()
  }) %>% bind_rows()
  
  # outputDetails("  Handling novel tracking, adding to dataset.\n", newcat = FALSE)
  isolates_file <- novelHandling(tp1, tp2, clusters_just_tp1, heights)
  
  isolates_file %<>% 
    mutate(novel = ifelse(isolate %in% setdiff(tp2$raw$isolate, tp1$raw$isolate), 1, 0)) %>% 
    rename(Strain = isolate)
  
  isolates_file[,c("tp1_h", "tp2_h")] %<>% apply(., 2, padCol, padval = ph, padchr = "h")
  isolates_file[,c("tp1_cl", "tp2_cl")] %<>% apply(., 2, padCol, padval = pc, padchr = "c")
  
  # outputDetails("  Incrementing all cluster sizes by 1, then calculating growth columns.\n", newcat = FALSE)
  # outputDetails("  Also adding 'type' column to CGM results table.\n", newcat = FALSE)
  isolates_file %<>% 
    mutate(tp1_cl_size = tp1_cl_size + 1, tp2_cl_size = tp2_cl_size + 1) %>% 
    oneHeight()
  
  # if (i > 1) {isolates_file <- unchanged_data %>% bind_rows(isolates_file, .)}
  
  # tmp <- isolates_file
  isolates_file %<>% addingType(.)
  
  # outputDetails("  Saving the data in a file with cluster identifiers.\n", newcat = FALSE)
  # strains removed
  isolates_file <- isolates_file %>% 
    select(tp1_id, tp1_cl_size, first_tp1_flag, last_tp1_flag, 
           first_tp2_flag, tp2_cl_size, last_tp2_flag, add_TP1, novel, 
           num_novs, actual_size_change, actual_growth_rate, new_growth, type) %>% 
    unique() %>% as.data.table()
  
  correcthere <- isolates_file
  
  indices <- which(is.na(isolates_file$tp1_id) & isolates_file$tp1_cl_size == 1)
  isolates_file[indices]$tp1_id <- isolates_file$first_tp2_flag[indices] %>% 
    sapply(., strsplit, "_") %>% sapply(., '[[', 3) %>% 
    paste0("AbsentAtTP1-", "TP2_", .)
  
  cgm_results <- isolates_file %>% 
    mutate(across(colnames(isolates_file), as.character)) %>% 
    melt.data.table(id.vars = "tp1_id") %>% 
    add_column(Interval = paste0(n1, "-", n2), .after = 1) %>% 
    set_colnames(c("Cluster", paste0(arg$intervaltype, "Interval"), "Field", "Value"))
  
  cgm_results <- isolates_file %>% 
    add_column(interval = paste0(n1, "-", n2), .before = 1)# %>% add_column(TP = n2, .before = 1)
  
  saveRDS(cgm_results, file.path(save_to, paste0("TP", n2, ".Rds")))
}

if (params$int_type[2] == "multiset") {
  res_file <- gsub("-", "", params$divs[2]) %>% gsub(",", "-", .) %>% 
    paste0("results/CGM-", ., "-midpoints.Rds")  
}else {
  res_file <- paste0("results/CGM-", params$int_type[2], "-intervals.Rds")
}

cgmfiles <- list.files(save_to, full.names = TRUE)
lapply(cgmfiles, function(fi) {readRDS(fi)}) %>% bind_rows() %>% saveRDS(., res_file)

# WRAPPING THINGS UP -------------------------------------------------------------------------------------------
stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())

outputDetails(paste0("\nSuccessfully collected data for all heights."), newcat = TRUE)
timeTaken(pt = "CGM data collection", stopwatch) %>% outputDetails(., newcat = TRUE)
outputDetails(paste0("||", paste0(rep("-", 28), collapse = ""), " End of cluster metric generation ", 
                     paste0(rep("-", 29), collapse = ""), "||\n"))
closeAllConnections()

# for 264052 strains, and 53 weeks, took 32 min and 10 sec to run the whole thing