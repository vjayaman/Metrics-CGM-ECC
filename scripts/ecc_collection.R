#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_epiquant.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

files <- paste0("scripts/ECC") %>% list.files(., full.names = TRUE)
invisible(sapply(files, source))

# Title: "EpiQuant - Salmonella Enteritidis Project (2019-2020)"
# Authors of original work and initial modifications: Ben Hetman, Elissa Giang, Dillon Barker
# Responsible for changes made during and after merge with CGM process: Vasena Jayamanna

checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

mergeECCs <- function(eccs, tpx, typing_data) {
  tbl1 <- as.data.table(typing_data)
  sapply(eccs, "[", tpx) %>% Reduce(function(...) merge(...), .) %>% as.data.table() %>% 
    merge.data.table(tbl1, ., by = intersect(colnames(tbl1), colnames(.))) %>% return()
}

option_list <- list(
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"),
  make_option(c("-t", "--trio"), metavar = "character", default = "010-001",
              help = "temporal, geographic coefficients"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ", 
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$trio %>% strsplit(., "-") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)

base_strains <- read_tsv(params$strains)
loc_cols <- length(intersect(c("Country", "Province", "City"), colnames(base_strains)))
if (loc_cols == 3) {
  strain_data <- base_strains %>% 
    mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")), 
           Location = paste(Country, Province, City, sep = "_")
    )
}else {
  strain_data <- base_strains %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")))
}

typing_data <- tp1$height_list %>% append(tp2$height_list)

cat(paste0("\n\nStep 1:"))
cat(paste0("\n   Note that the source coefficent is always 0 in this version"))

# Note: dr stands for data representative
# in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
outputMessages("   Removing redundancy (comparing date, lat-long, etc. - not every pair of strains)")
if (loc_cols == 3) {
  assignments <- strain_data %>% select(Date, Latitude, Longitude, Location)  
}else {
  assignments <- strain_data %>% select(Date, Latitude, Longitude)
}
assignments %<>% unique() %>% rownames_to_column("dr")

outputMessages("   Generating all possible date pair distances ...")
uni_temp <- assignments %>% select(Date) %>% unique() %>% rownames_to_column("drt") %>% 
  left_join(., assignments, by = "Date")
dm_temp <- uni_temp %>% select(drt, Date) %>% unique() %>% rename(dr = drt) %>% 
  distMatrix(., "temp", "Date")
transformed_temp <- dm_temp %>% pairwiseDists(., "temp", c("drt1", "drt2", "Temp.Dist"))

outputMessages("   Generating all possible lat-long pair distances ...")
uni_geo <- assignments %>% select(Latitude, Longitude) %>% unique() %>% 
  rownames_to_column("drg") %>% 
  left_join(., assignments, by = c("Latitude", "Longitude"))
dm_geo <- uni_geo %>% select(drg, Latitude, Longitude) %>% unique() %>% 
  rename(dr = drg) %>% distMatrix(., "geo", c("Latitude", "Longitude"))
transformed_geo <- dm_geo %>% pairwiseDists(., "geo", c("drg1", "drg2", "Geog.Dist"))

dr_tempgeo <- uni_geo %>% 
  full_join(., uni_temp, by = c("dr", "Date", "Latitude", "Longitude", "Location")) %>% 
  select(dr, drt, drg)

# transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
# rm(transformed_geo)
# rm(transformed_temp)


# outputMessages("   Generating all possible date pair distances ...")
# dm_temp <- assignments %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
# transformed_temp <- dm_temp %>% pairwiseDists(., "temp", c("dr1", "dr2", "Temp.Dist"))
# 
# outputMessages("   Generating all possible lat-long pair distances ...")
# dm_geo <- assignments %>% select(dr, Latitude, Longitude) %>% 
#   distMatrix(., "geo", c("Latitude", "Longitude"))
# transformed_geo <- dm_geo %>% pairwiseDists(., "geo", c("dr1", "dr2", "Geog.Dist"))
# 
# transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
# rm(transformed_geo)
# rm(transformed_temp)

outputMessages("   Identifying which strains match which non-redundant 'data representatives'")
if (loc_cols == 3) {
  match_names <- c("Latitude", "Longitude", "Date", "Location")
}else {
  match_names <- c("Latitude", "Longitude", "Date")
}

dr_matches <- left_join(strain_data, assignments, by = match_names) %>% 
  select(Strain, dr) %>% 
  left_join(., dr_tempgeo, by = "dr")

t2clusters <- typing_data[[2]] %>% set_colnames("Threshold") %>% 
  rownames_to_column("Strain") %>% as_tibble()

# avgdistvals <- lapply(1:length(typing_data), function(i) {
#   dr_td1 <- typing_data[[i]] %>% rownames_to_column("Strain") %>% as_tibble() %>%
#     left_join(., dr_matches, by = "Strain") %>%
#     mutate(across(dr, as.character)) %>% select(-Strain)
#   
#   # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
#   cx <- colnames(dr_td1)[1]
#   tallied_reps <- dr_td1 %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
#   g_cuts <- left_join(dr_td1, tallied_reps, by = intersect(colnames(tallied_reps), colnames(dr_td1))) %>%
#     unique() %>% mutate(across(dr, as.character))
#   
#   outputMessages(paste0("      Calculating average (not transformed) distances for timepoint ", i))
#   a2 <- avgDists(g_cuts, dm_temp, "Temp.Dist", paste0("TP", i, "_", colnames(g_cuts)[1]))
#   b2 <- avgDists(g_cuts, dm_geo, "Geog.Dist", paste0("TP", i, "_", colnames(g_cuts)[1]))
#   return(list(temp = a2, geo = b2))
# })


collected_eccs <- lapply(1:length(combos), function(j) {
  c1 <- unlist(strsplit(combos[j], split = "")) %>% as.numeric()
  cat(paste0("\n\nStep ", j + 1, ":"))
  tau <- c1[2]
  gamma <- c1[3]
  transformed_temp %<>% mutate(Temp.New = (Temp.Dist^2)*tau)
  transformed_geo %<>% mutate(Geog.New = (Geog.Dist^2)*gamma)
  epiCollection(strain_data, tau, gamma, typing_data, 
                transformed_dists, dm_temp, dm_geo, dr_matches, avgdistvals, j)
})

cat(paste0("\n\nStep ", length(combos) + 2, ":"))
outputMessages("   Merging collected ECCs ...\n")
full_set <- mergeECCs(collected_eccs, 1, tp1$proc) %>%
  merge.data.table(., mergeECCs(collected_eccs, 2, tp2$proc), by = "Strain", all.y = TRUE) %>%
  mutate(TP1 = ifelse(is.na(TP1), 0, TP1))

write.table(full_set, file = "results/ECCs.tsv", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))
