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
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"),
  make_option(c("-t", "--duo"), metavar = "character", default = "01-10",
              help = "temporal, geographic coefficients"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ", 
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$duo %>% strsplit(., "-") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)

cat(paste0("\n\nPart 1:"))
# outputMessages("   Removing Type 4 cases before ECC generation (will handle separately) ...")
# # Type IV modifications: TP1 < 3, TP2 < 3
# tp1data <- tp1$height_list[[1]] %>% rownames_to_column() %>% 
#   as_tibble() %>% set_colnames(c("Strain", "tp1_cl"))
# 
# tp2data <- tp2$height_list[[1]] %>% rownames_to_column() %>% 
#   as_tibble() %>% set_colnames(c("Strain", "tp2_cl"))
# 
# tp1clusters <- tp1data %>% group_by(tp1_cl) %>% summarise(n = n()) %>% filter(n < 3) %>% left_join(., tp1data)
# tp2clusters <- tp2data %>% group_by(tp2_cl) %>% summarise(n = n()) %>% filter(n < 3) %>% left_join(., tp2data)
# type4_strains <- intersect(tp1clusters$Strain, tp2clusters$Strain)
# 
# tp1singletons <- tp1clusters %>% filter(n == 1) %>% pull(Strain)
# tp2singletons <- tp2clusters %>% filter(n == 1) %>% pull(Strain)
# 
# tp1$height_list[[1]] %<>% filter(!(rownames(.) %in% c(tp1singletons, type4_strains)))
# tp2$height_list[[1]] %<>% filter(!(rownames(.) %in% c(tp2singletons, type4_strains)))

typing_data <- tp1$height_list %>% append(tp2$height_list)

strain_data <- read_tsv(params$strains) %>% 
  mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")), 
         Location = paste(Country, Province, City, sep = "_")) %>% 
  filter(Strain %in% rownames(typing_data[[2]]))

cat(paste0("\n\nPart 2:"))
cat(paste0("\nNote that source = 0"))

# Note: dr stands for data representative
# in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
outputMessages("   Removing redundancy (comparing date, lat-long, etc. - not every pair of strains)")
assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
  unique() %>% rownames_to_column("dr")

outputMessages("   Generating all possible date pair distances ...")
dm_temp <- assignments %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
transformed_temp <- dm_temp %>% pairwiseDists(., "temp", c("dr1", "dr2", "Temp.Dist"))

outputMessages("   Generating all possible lat-long pair distances ...")
dm_geo <- assignments %>% select(dr, Latitude, Longitude) %>% distMatrix(., "geo", c("Latitude", "Longitude"))
transformed_geo <- dm_geo %>% pairwiseDists(., "geo", c("dr1", "dr2", "Geog.Dist"))

transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
rm(transformed_geo)
rm(transformed_temp)

outputMessages("   Identifying which strains match which non-redundant 'data representatives'")
dr_matches <- strain_data %>% 
  left_join(., assignments, by = c("Latitude", "Longitude", "Date", "Location")) %>% 
  select(Strain, dr)

avgdistvals <- lapply(1:length(typing_data), function(i) {
  dr_td1 <- typing_data[[i]] %>% rownames_to_column("Strain") %>% as_tibble() %>%
    left_join(., dr_matches, by = "Strain") %>%
    mutate(across(dr, as.character)) %>% select(-Strain)
  
  # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
  cx <- colnames(dr_td1)[1]
  tallied_reps <- dr_td1 %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
  g_cuts <- left_join(dr_td1, tallied_reps, by = intersect(colnames(tallied_reps), colnames(dr_td1))) %>%
    unique() %>% mutate(across(dr, as.character))
  
  outputMessages(paste0("      Calculating average (not transformed) distances for timepoint ", i))
  a2 <- avgDists(g_cuts, dm_temp, "Temp.Dist", paste0("TP", i, "_", colnames(g_cuts)[1]))
  b2 <- avgDists(g_cuts, dm_geo, "Geog.Dist", paste0("TP", i, "_", colnames(g_cuts)[1]))
  return(list(temp = a2, geo = b2))
})

collected_eccs <- lapply(1:length(combos), function(j) {
  c1 <- unlist(strsplit(combos[j], split = "")) %>% as.numeric()
  tau <- c1[1]
  gamma <- c1[2]
  cat(paste0("\n\nPart ", j + 2, ":"))
  epiCollection(strain_data, tau, gamma, typing_data, transformed_dists, 
                dm_temp, dm_geo, dr_matches, avgdistvals)
})

cat(paste0("\n\nPart ", length(combos) + 3, ":"))
outputMessages("Merging collected ECCs ...")
full_set <- mergeECCs(collected_eccs, 1, tp1$proc) %>% 
  merge.data.table(., mergeECCs(collected_eccs, 2, tp2$proc), by = "Strain", all.y = TRUE) %>% 
  mutate(TP1 = ifelse(is.na(TP1), 0, TP1))

# outputMessages("Handling Type 4 cases, adding into table before saving ...")

write.table(full_set, file = "results/ECCs.tsv", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))

