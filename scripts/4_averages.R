#! /usr/bin/env Rscript

msg <- file("logs/avg_dists.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "Rcpp", "RcppArmadillo")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y))); rm(libs); rm(y)

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " Average distances collection ",
           paste0(rep("-", 31), collapse = ""), "||\nStarted process at: ", Sys.time(), "\n"))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source)); rm(files)

outputDetails("Sourcing required functions, reading in input data", newcat = TRUE)

sourceCpp("scripts/epicohversions.cpp")

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", 
              help = "Time point 2 file name (TP2)"), 
  make_option(c("-f", "--intervalfile"), metavar = "file", default = "inputs/processed/clustersets.Rds"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))

arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)

# dm <- dms[["temp"]]; cnames <- "Date"
groupXAvgDists <- function(dm, groupX, assignments, strain_data, cnames, metadata) {
  clustersX <- pull(groupX, 1)

  avgdists <- lapply(clustersX, function(x) {
    cx <- assignments[Hx == x]

    num_unique_cases <- strain_data[Strain %in% cx$Strain] %>%
      select(all_of(cnames)) %>% unique() %>% nrow()
    if (num_unique_cases == 1) {
      return(0)
    }else {
      dr_counts <- cx %>% group_by(dr) %>% summarise(n_drs = n()) %>% as.data.table()
      dmx <- dm[rownames(dm) %in% dr_counts$dr, colnames(dm) %in% dr_counts$dr,drop=FALSE]
      
      dr_counts <- dr_counts[match(rownames(dmx), dr_counts$dr)]
      dmxcols <- sweep(dmx, MARGIN = 2, dr_counts$n_drs, `*`)
      dmxrest <- sweep(dmxcols, MARGIN = 1, dr_counts$n_drs, `*`)
      dmxrest <- dmxrest[upper.tri(dmxrest)]
      
      allstrains <- matrix(nrow = sum(dr_counts$n_drs), ncol = sum(dr_counts$n_drs))
      avg_dist <- sum(dmxrest)/length(allstrains[upper.tri(allstrains)])
      return(avg_dist)
    }
    # if (is.na(avg_dist)) {print(x)}
  }) %>% unlist() %>% data.table(Hx = clustersX, AvgDist = .)
  assert("No NA avg dists", nrow(avgdists[is.na(AvgDist)]) == 0)
  return(avgdists)
}

outputDetails("Processing input data, formatting ...", newcat = TRUE)

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
interval_list <- names(clustersets); rm(clustersets)
basedir <- file.path("intermediate_data", params$int_type[2], "avgdists")
parts <- readRDS("intermediate_data/TPN/parts.Rds")
groups <- names(parts$results)

tpnfiles <- lapply(names(parts$results), function(x) parts$results[[x]] %>% add_column(f = x) %>% 
                     select(-n)) %>% bind_rows() %>% as.data.table() %>% set_colnames(c("Hx", "f"))

strain_data <- as.data.table(m$strain_data)

typing_data <- lapply(1:length(interval_list), function(i) {
  n1 <- as.character(interval_list[i])
  tpkstrains <- metadata[get(interval) <= n1]$Strain
  dfz <- tp2$filedata %>% rownames_to_column("isolate") %>%
    select(isolate, all_of(hx$h)) %>%
    filter(isolate %in% tpkstrains) %>% column_to_rownames("isolate")
  dfz[,hx$h[1],drop=FALSE] %>% set_colnames(hx$th[1])
}) %>% set_names(as.character(interval_list))

outputDetails(paste0("\nCollecting average distances (for each cluster) at ..."), newcat = TRUE)

for (tdx in names(typing_data)) {
  outputDetails(paste0("TP", tdx, ":"), newcat = TRUE)
  
  td <- typing_data[[tdx]] %>% rownames_to_column("Strain") %>% as.data.table()
  pd1 <- parts$drs[Strain %in% td$Strain] %>% set_colnames(c("Strain", "Hx", "dr")) %>% 
    inner_join(., tpnfiles, by = "Hx") %>% arrange(f)
  pd2 <- pd1 %>% group_by(Hx) %>% summarise(n = n()) %>% left_join(pd1, ., by = "Hx") %>% 
    select(-Strain) %>% unique()
  
  assignments <- inner_join(m$dr_matches, td, by = "Strain") %>% 
    set_colnames(c("Strain", "dr", "Hx")) %>% as.data.table()
  
  sizes <- td %>% set_colnames(c("Strain", "Hx")) %>% group_by(Hx) %>% summarise(Size = n())
  
  pb <- txtProgressBar(min = 0, max = length(groups), initial = 0, style = 3)
  avg_dists <- lapply(1:length(groups), function(i) {
    setTxtProgressBar(pb, i)
    # print(i)
    x_i <- groups[i]
    dms <- paste0("group", x_i, ".Rds") %>% file.path("intermediate_data/TPN/dists", .) %>% readRDS()
    groupX <- pd2[f == x_i] %>% select(Hx, n) %>% unique()
    
    if (nrow(groupX) > 0) {
      temp_dists <- dms[["temp"]] %>% 
        groupXAvgDists(., groupX, assignments, strain_data, "Date", metadata) %>% 
        rename(Temp.Avg.Dist = AvgDist)
      
      geo_dists <- dms[["geo"]] %>% 
        groupXAvgDists(., groupX, assignments, strain_data, c("Longitude", "Latitude"), metadata) %>% 
        rename(Geo.Avg.Dist = AvgDist)  
      
      return(merge.data.table(temp_dists, geo_dists))
    }
  }) %>% bind_rows() %>% inner_join(., sizes, by = "Hx")
  close(pb)
  
  outputDetails(paste0("   Checking that only clusters with just one unique info pair have average distance 0"), newcat = TRUE)
  assert("No NA avg geo distances", !any(is.na(avg_dists$Geo.Avg.Dist)))
  assert("No NA avg temp distances", !any(is.na(avg_dists$Temp.Avg.Dist)))
  
  # Check that clusters with the average pairwise distance being 0 have only one
  # unique date / coordinate pair, as is required to get a distance result of 0
  renamed_tping <- td %>% set_colnames(c("Strain", "Hx"))
  
  temp_zero_clusters <- avg_dists[Temp.Avg.Dist == 0]$Hx
  for (x in temp_zero_clusters) {
    y <- renamed_tping[Hx == x] %>% pull(Strain)
    unidate <- strain_data[Strain %in% y, Date] %>% unique() %>% length()
    assert(paste0("Cluster ", x, " has only one unique date"), unidate == 1)
  }
  
  geo_zero_clusters <- avg_dists[Geo.Avg.Dist == 0]$Hx
  for (x in geo_zero_clusters) {
    y <- renamed_tping[Hx == x] %>% pull(Strain)
    uniloc <- strain_data[Strain %in% y] %>% select(Latitude, Longitude) %>% unique() %>% nrow()
    assert(paste0("Cluster ", x, " has only one unique coordinate pair"), uniloc == 1)
  }

  colnames(avg_dists)[which(colnames(avg_dists) == "Hx")] <- hx$th
  
  # Collects raw cluster averages as well
  curr_data <- left_join(td, strain_data, by = "Strain") %>% 
    set_colnames(gsub(hx$th, "heightx", colnames(.)))
  clusters <- curr_data %>% pull(heightx) %>% unique()
  
  reg_avgs <- lapply(clusters, function(c_i) {
    dfx <- curr_data[heightx == c_i]
    data.table(c_i, mean(dfx$Date), mean(dfx$Longitude), mean(dfx$Latitude))
  }) %>% bind_rows() %>% 
    set_colnames(c(hx$th, "Avg.Date", "Avg.Longitude", "Avg.Latitude"))
  
  merge.data.table(avg_dists, reg_avgs) %>% 
    set_colnames(gsub("Size", paste0(hx$th, "_Size"), colnames(.))) %>% 
    add_column(TP = tdx, .before = 1) %>% 
    saveRDS(., file.path(basedir, paste0("TP", tdx, ".Rds")))
}

all_avg_dists <- lapply(names(typing_data), function(tdx) {
  readRDS(file.path(basedir, paste0("TP", tdx, ".Rds")))
}) %>% set_names(names(typing_data)) %>% bind_rows()


if (params$int_type[2] == "multiset") {
  res_file <- gsub("-", "", params$divs[2]) %>% gsub(",", "-", .) %>% 
    paste0("results/AVGS-", ., "-midpoints.Rds")  
}else {
  res_file <- paste0("results/AVGS-", params$int_type[2], "-intervals.Rds")
}

outputDetails(paste0("\nMerging average distance results and saving to '", res_file, "'"), newcat = TRUE)

saveRDS(all_avg_dists, res_file)

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())

timeTaken(pt = "average distances collection", stopwatch) %>% outputDetails(., newcat = TRUE)
cat(paste0("||", paste0(rep("-", 27), collapse = ""), " End of average distances collection ",
           paste0(rep("-", 27), collapse = ""), "||\n"))