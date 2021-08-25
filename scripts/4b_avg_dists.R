#! /usr/bin/env Rscript

msg <- file("logs/avg_dists.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "Rcpp", "RcppArmadillo")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y))); rm(libs); rm(y)

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)
sourceCpp("scripts/epicohversions.cpp")

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source)); rm(files)

groupXAvgDists <- function(dm, groupX, assignments) {
  clustersX <- pull(groupX, 1)
  
  lapply(clustersX, function(x) {
    cx <- assignments[Hx == x]
    if (length(unique(cx$dr)) == 1) {return(0)
    }else {
      dr_counts <- cx %>% group_by(dr) %>% summarise(n_drs = n()) %>% as.data.table()
      dmx <- dm[rownames(dm) %in% dr_counts$dr, colnames(dm) %in% dr_counts$dr]
      dr_counts <- dr_counts[match(rownames(dmx), dr_counts$dr)]
      dmxcols <- sweep(dmx, MARGIN = 2, dr_counts$n_drs, `*`)
      dmxrest <- sweep(dmxcols, MARGIN = 1, dr_counts$n_drs, `*`)
      return(mean(dmxrest))
    }
  }) %>% unlist() %>% data.table(Hx = clustersX, AvgDist = .) %>% return()
}

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ",
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", 
              help = "Time point 2 file name (TP2)"), 
  make_option(c("f", "--intervalfile"), metavar = "file", default = "inputs/processed/clustersets.Rds"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))

arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)

params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

hx <- strsplit(as.character(params$th[2]), split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))
tp2 <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)
m <- read_tsv(arg$metadata) %>% processedStrains()
metadata <- m$strain_data %>% as.data.table()

if (params$int_type[2] == "multiset") {
  interval <- "Multiset"
}else if (params$int_type[2] == "monthly") {
  interval <- "YearMonth"
}else if (params$int_type[2] == "weekly") {
  interval <- "Week"
}

clustersets <- readRDS(arg$intervalfile)
interval_list <- names(clustersets)
rm(clustersets)

basedir <- file.path("intermediate_data", params$int_type[2], "avgdists")
dir.create(basedir, showWarnings = FALSE, recursive = TRUE)

parts <- readRDS("intermediate_data/TPN/parts.Rds")
groups <- names(parts$results)

typing_data <- lapply(1:length(interval_list), function(i) {
  n1 <- as.character(interval_list[i])
  tpkstrains <- metadata[get(interval) <= n1]$Strain
  dfz <- tp2$filedata %>% rownames_to_column("isolate") %>%
    select(isolate, all_of(hx$h)) %>%
    filter(isolate %in% tpkstrains) %>% column_to_rownames("isolate")
  dfz[,hx$h[1],drop=FALSE] %>% set_colnames(hx$th[1])
}) %>% set_names(as.character(interval_list))


for (tdx in names(typing_data)) {
  print(tdx)
  td <- typing_data[[tdx]] %>% rownames_to_column("Strain") %>% as.data.table()
  
  assignments <- inner_join(m$dr_matches, td, by = "Strain") %>% 
    set_colnames(c("Strain", "dr", "Hx")) %>% as.data.table()
  
  sizes <- td %>% set_colnames(c("Strain", "Hx")) %>% group_by(Hx) %>% summarise(Size = n())
  
  temp_avg_dists <- lapply(groups, function(x_i) {
    dms <- paste0("group", x_i, ".Rds") %>% file.path("intermediate_data/TPN/dists", .) %>% readRDS()
    groupX <- parts$results[[x_i]] %>% as.data.table() %>% set_colnames(c("Hx", "n"))
    groupX <- groupX[Hx %in% assignments[Strain %in% td$Strain]$Hx]
    groupXAvgDists(dms[["temp"]], groupX, assignments)
  }) %>% bind_rows() %>% inner_join(., sizes, by = "Hx") %>% rename(Temp.Avg.Dist = AvgDist)

  geo_avg_dists <- lapply(groups, function(x_i) {
    dms <- paste0("group", x_i, ".Rds") %>%
      file.path("intermediate_data/TPN/dists", .) %>% readRDS()
    groupX <- parts$results[[x_i]] %>% as.data.table() %>% set_colnames(c("Hx", "n"))
    groupX <- groupX[Hx %in% assignments[Strain %in% td$Strain]$Hx]
    groupXAvgDists(dms[["geo"]], groupX, assignments)
  }) %>% bind_rows() %>% inner_join(., sizes, by = "Hx") %>% rename(Geo.Avg.Dist = AvgDist)
  
  assert("No NA avg geo distances", nrow(geo_avg_dists[is.na(Geo.Avg.Dist)]) == 0)
  assert("No NA avg temp distances", nrow(temp_avg_dists[is.na(Temp.Avg.Dist)]) == 0)
  
  # Check that clusters with size > 1 where the average pairwise distance is 0 have only one
  # unique date / coordinate pair, as is required to get a distance result of 0
  renamed_tping <- td %>% set_colnames(c("Strain", "Hx"))
  
  temp_zero_clusters <- temp_avg_dists[Temp.Avg.Dist == 0 & Size > 1]$Hx
  for (x in temp_zero_clusters) {
    y <- renamed_tping[Hx == x] %>% pull(Strain)
    unidate <- m$strain_data %>% filter(Strain %in% y) %>% pull(Date) %>% unique() %>% length()
    assert(paste0("Cluster ", x, " has only one unique date"), unidate == 1)
  }
  
  geo_zero_clusters <- geo_avg_dists[Geo.Avg.Dist == 0 & Size > 1]$Hx
  for (x in geo_zero_clusters) {
    y <- renamed_tping[Hx == x] %>% pull(Strain)
    uniloc <- m$strain_data %>% filter(Strain %in% y) %>% select(Latitude, Longitude) %>% unique() %>% nrow()
    assert(paste0("Cluster ", x, " has only one unique coordinate pair"), uniloc == 1)
  }
  
  inner_join(temp_avg_dists, geo_avg_dists, by = c("Hx", "Size")) %>% 
    add_column(TP = tdx, .before = 1) %>% 
    saveRDS(., file.path(basedir, paste0("TP", tdx, ".Rds")))
}

avg_dists <- lapply(names(typing_data), function(tdx) {
  readRDS(file.path(basedir, paste0("TP", tdx, ".Rds")))
}) %>% set_names(names(typing_data)) %>% bind_rows()

if (params$int_type[2] == "multiset") {
  res_file <- gsub("-", "", params$divs[2]) %>% gsub(",", "-", .) %>% 
    paste0("results/AVGS-", ., "-midpoints.Rds")  
}else {
  res_file <- paste0("results/AVGS-", params$int_type[2], "-intervals.Rds")
}

saveRDS(avg_dists, res_file)
