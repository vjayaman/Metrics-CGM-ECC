#! /usr/bin/env Rscript

msg <- file("logs/logfile_epiquant1.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source))

assert("Distances were collected and saved", file.exists("intermediate_data/dist_extremes.Rds"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ",
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))

option_list <- list(
  make_option(c("-m", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-t", "--trio"), metavar = "character", default = "0-1-0",
              help = "source, temporal, geographic coefficients"), 
  make_option(c("-i", "--intervaltype"), metavar = "char", default = "weekly", 
              help = "Type of intervals, choices are: weekly, monthly, multiset. If multiset, provide a time to split the dataset at."))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$trio %>% strsplit(., ",") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))
m <- read_tsv(params$strains) %>% processedStrains()
metadata <- m$strain_data %>% 
  mutate(YearMonth = format(Date, "%Y-%m")) %>% 
  mutate(wk = strftime(Date, format = "%V")) %>% 
  select(-TP1, -TP2) %>% arrange(wk) %>% as.data.table()

tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
f2 <- tp2$filedata %>% rownames_to_column("isolate") %>% as.data.table() %>% 
  select(isolate, all_of(params$heights)) %>% arrange(isolate)

if (params$intervaltype == "weekly") {
  print("weekly")
  interval <- "Week"
  interval_list <- metadata$wk %>% unique() %>% sort()
  
  interval_clusters <- metadata %>% select(c(Strain, wk)) %>% 
    inner_join(., f2, by = c("Strain" = "isolate")) %>% 
    set_colnames(c("isolate", "ivl", "heightx"))
  
}else if (params$intervaltype == "monthly") {
  print("monthly")
  interval <- "Month"
  interval_list <- metadata$YearMonth %>% unique() %>% sort()
  
  interval_clusters <- metadata %>% select(c(Strain, YearMonth)) %>% 
    inner_join(., f2, by = c("Strain" = "isolate")) %>% 
    set_colnames(c("isolate", "ivl", "heightx"))
  

}else if (params$intervaltype == "multiset") {
  if (interactive()) {
    divider <- paste0("Enter date(s) between ", min(metadata$Date), " and ", 
                      max(metadata$Date), " to split set into two timepoints ", 
                      "\n(YYYY-MM-DD format, separated by commas): ") %>% 
      readline() %>% as.character()
  }else {
    cat(paste0("Enter date(s) between ", min(metadata$Date), " and ", 
               max(metadata$Date), " to split set into two timepoints ", 
               "\n(YYYY-MM-DD format, separated by commas): "))
    divider <- readLines("stdin", n = 1) %>% as.character()
  }
  
  setdivider <- strsplit(divider, split = ",") %>% unlist() %>% as.Date(., format = "%Y-%m-%d")
  
  assertion1 <- lapply(setdivider, function(x) nchar(as.character(x)) == 10) %>% unlist()
  assert("Provided input(s) has/have 10 characters", all(assertion1))
  
  assert("Provided input(s) is/are date type, correct format", !any(is.na(setdivider)))
  
  assertion3 <- lapply(setdivider, function(x) x >= min(metadata$Date) & x <= max(metadata$Date)) %>% unlist()
  assert("Provided input(s) is/are between min date and max date", all(assertion3))
  
  sdiv <- c(min(metadata$Date), setdivider, max(metadata$Date))
  fdivs <- tibble(start = sdiv[1:(length(sdiv)-1)], end = sdiv[2:length(sdiv)]) %>% 
    mutate(ivl = paste0("set", 1:nrow(.)))

  metadata <- metadata %>% add_column(Tpt = "")
  for (i in 1:nrow(fdivs)) {
    metadata[metadata$Date >= fdivs$start[i] & metadata$Date < fdivs$end[i]]$Tpt <- fdivs$ivl[i]
  }
  metadata[metadata$Date == fdivs$end[i]]$Tpt <- fdivs$ivl[i]
  
  interval <- "Multiset"
  interval_list <- metadata$Tpt %>% unique() %>% sort()

  interval_clusters <- metadata %>% select(c(Strain, Tpt)) %>%
    inner_join(., f2, by = c("Strain" = "isolate")) %>%
    set_colnames(c("isolate", "ivl", "heightx"))
}

clusters <- vector(mode = "list", length = length(interval_list)) %>% set_names(interval_list)

for (xj in interval_list) {
  # cluster assignments for clusters that changed when interval i strains were added
  int_j <- interval_clusters[heightx %in% interval_clusters[ivl == xj]$heightx]
  sofar <- interval_clusters[heightx %in% interval_clusters[ivl <= xj]$heightx]
  clusters[[xj]] <- list(int_j, sofar) %>% set_names(c("ivl", "sofar"))
}

for (i in 1:length(interval_list)) {
  n1 <- as.character(interval_list[i])
  
  tp1strains <- metadata[wk <= n1]$Strain

  dfz <- tp2$filedata %>% rownames_to_column("isolate") %>% 
    select(isolate, all_of(params$heights)) %>% 
    filter(isolate %in% tp1strains) %>% column_to_rownames("isolate")
  
  typing_data <- list(dfz[,hx$h[1],drop=FALSE] %>% set_colnames(hx$th[i])) %>% set_names(n1)

  extremes <- readRDS("intermediate_data/dist_extremes.Rds")

  dfx <- expand.grid(x = combos, k = n1) %>% as.data.frame()
  
  for (k in unique(dfx$k)) {
    dir.create(paste0("intermediate_data/TP", k), showWarnings = FALSE)
    dir.create(paste0("intermediate_data/TP", k, "/ecc_groups/"), showWarnings = FALSE)
  }

  stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)
  
  for (i in 1:nrow(dfx)) {
    cat(paste0("\n\nStep ", i, " / ", nrow(dfx) + 2, ":"))
    c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "-") %>% unlist() %>% as.numeric()
  
    paste0("intermediate_data/TP", dfx$k[i], "/ecc_groups/", as.character(dfx$x[i])) %>% 
      dir.create(., showWarnings = FALSE)
  
    outputMessages(paste0("\nCollecting and saving ECCs for groups of clusters at TP", 
                          dfx$k[i], ", for ECC coefficient triple ", as.character(dfx$x[i])))
    
    k <- dfx$k[i]
    
    parts <- m$dr_matches %>% filter(Strain %in% rownames(typing_data$`2`[1])) %>% 
      sectionClusters(k, typing_data, .)
    read_from <- paste0("intermediate_data/TP", 2, "/dists/")
    save_to <- paste0("intermediate_data/TP", k, "/ecc_groups/", as.character(dfx$x[i]), "/")
  
    collectECCs(k, m, parts, extremes, c1, read_from, save_to)
  
    c1 %>% paste0(., collapse = ", ") %>% 
      paste0("\nMerging ECC files at TP", dfx$k[i], ", for ECC parameters ", .) %>% outputMessages()
  }

  cat(paste0("\n\nStep ", nrow(dfx) + 1, " / ", nrow(dfx) + 2, ":"))
  for (i in 1:nrow(dfx)) {
    fnames <- list.files(paste0("intermediate_data/TP", dfx$k[i], 
                                "/ecc_groups/", dfx$x[i], "/"), full.names = TRUE)
    
    tpk <- lapply(fnames, function(f) {readRDS(f)}) %>% bind_rows()
    
    saveRDS(tpk, paste0("intermediate_data/TP", dfx$k[i], "/", dfx$x[i], "-eccs.Rds"))
  }
}

# typing_data <- tp2$height_list %>% set_names("2")


# Generating ECC results file ----------------------------------------------------
# source("scripts/5_merging_eccs.R")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ",
           paste0(rep("-", 31), collapse = ""), "||\n"))
