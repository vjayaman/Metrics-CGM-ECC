#! /usr/bin/env Rscript

libs <- c("optparse","magrittr","tibble", "dplyr", "readr", "testit", "data.table", "R6")
y <- suppressMessages(lapply(libs, require, character.only = TRUE))

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-b", "--tp2"), metavar = "file", 
              default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"), 
  make_option(c("f", "--intervalfile"), metavar = "file", 
              default = "inputs/processed/clustersets.Rds"))

arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)

source("scripts/Misc/type_handling.R"); source("scripts/ECC/classes_ecc.R")

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), " Merging CGM and ECC results ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))

# INPUT PREP -------------------------------------------------------------------------------------------------
params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs","coeffs", "numcl"))

hx <- strsplit(as.character(params$th[2]), split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

if (params$int_type[2] == "multiset") {
  interval <- "Multiset"
}else if (params$int_type[2] == "monthly") {
  interval <- "YearMonth"
}else if (params$int_type[2] == "weekly") {
  interval <- "Week"
}

clustersets <- readRDS(arg$intervalfile)
interval_list <- names(clustersets) %>% sort(); rm(clustersets)

# ------------------------------------------------------------------------------------------------------------
# NOW SAVING OUTPUTS AND MERGING ECCS WITH CGM DATA ----------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# For padding height and cluster columns with h0..0.., and c0..0.., respectively
padCol <- function(cvals, padval, padchr) {
  ifelse(!is.na(cvals), formatC(cvals, width = padval, format = "d", flag = "0") %>% 
           paste0(padchr, .), NA) %>% return()
}

# read in CGM results, and use "first_tp2_flag" get the padding values for height and cluster
cgms <- readRDS("results/CGM-monthly-intervals.Rds") %>% select(-TP)

a1 <- unlist(strsplit(cgms$first_tp2_flag[1], split = "_"))[2:3] %>% 
  nchar() %>% `-`(1) %>% set_names(c("ph", "pc"))
h_id <- padCol(as.double(hx$h), a1[['ph']], "h")

# read in AVERAGE results and use the TP column to create a TP1_ID column and an identical TP2_ID column: 
# since each TP set is both a TP1 and a TP2, depending on which two timepoints we are comparing
avgs <- readRDS("results/AVGS-monthly-intervals.Rds")
colnames(avgs)[grep("Avg", colnames(avgs))] %<>% paste0(hx$th, "_", .)
avgs <- avgs %>% 
  mutate(tp1_id = paste0("TP1_", h_id, "_", padCol(!!as.symbol(hx$th), a1[['pc']], "c"))) %>% 
  mutate(first_tp2_flag = tp1_id %>% gsub("TP1", "TP2", .))

# read in ECC results and repeat the process
eccs <- readRDS("results/ECC-monthly-intervals.Rds") %>% 
  mutate(tp1_id = paste0("TP1_", h_id, "_", padCol(!!as.symbol(hx$th), a1[['pc']], "c"))) %>% 
  mutate(first_tp2_flag = tp1_id %>% gsub("TP1", "TP2", .))

assert("ECC results have data for all intervals", identical(unique(eccs$TP), interval_list))
assert("AVG results have data for all intervals", identical(unique(avgs$TP), interval_list))

# cgms[,hx$th] <-
#   lapply(cgms$first_tp2_flag, function(idx) {
#     strsplit(idx, split = "c") %>% unlist() %>% extract2(2) %>% as.double()
#   }) %>% unlist()

step0 <- merge.data.table(eccs, avgs, by = intersect(colnames(eccs), colnames(avgs)))
stepcols <- ncol(step0)

step1 <- lapply(1:(length(interval_list)-1), function(i) {
  ivl1 <- interval_list[i]
  ivl2 <- interval_list[i+1]

  # selecting columns tp1_id, Size, ECC.0.1.0, ECC.0.0.1, Temp.Avg.Dist, Geo.Avg.Dist
  partA <- step0[TP == ivl1] %>% select(4,3,6:all_of(stepcols)) %>% 
    set_colnames(gsub(hx$th, "TP1", colnames(.)))
  
  # selecting columns first_tp2_flag, Size, ECC.0.1.0, ECC.0.0.1, Temp.Avg.Dist, Geo.Avg.Dist
  partB <- step0[TP == ivl2] %>% select(5,3,6:all_of(stepcols)) %>% 
    set_colnames(gsub(hx$th, "TP2", colnames(.)))
  
  partC <- cgms[interval == paste0(c(ivl1, ivl2), collapse = "-")] %>% 
    merge.data.table(., partA, all.x = TRUE) %>% 
    merge.data.table(., partB, by = "first_tp2_flag")
  
  partC %>% select(
    interval, tp1_id, tp1_cl_size, TP1_Size, TP1_ECC.0.1.0, TP1_ECC.0.0.1, 
    first_tp2_flag, tp2_cl_size, TP2_Size, TP2_ECC.0.1.0, TP2_ECC.0.0.1, 
    # delta_ECC_0.1.0, delta_ECC_0.0.1, 
    TP1_Avg.Date, TP1_Temp.Avg.Dist, TP1_Avg.Longitude, 
    TP1_Avg.Latitude, TP1_Geo.Avg.Dist, TP2_Avg.Date, TP2_Temp.Avg.Dist, 
    TP2_Avg.Longitude, TP2_Avg.Latitude, TP2_Geo.Avg.Dist, first_tp1_flag, 
    last_tp1_flag, first_tp2_flag, last_tp2_flag, tp1_cl_size, tp2_cl_size, 
    actual_size_change, add_TP1, novel, num_novs, actual_growth_rate, new_growth, type    )
}) %>% bind_rows()

tpn <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)

cnames <- params$cnames[2] %>% strsplit(split = ",") %>% unlist()
metadata <- read_tsv(arg$metadata) %>% as.data.table() %>%
  select(Strain, Latitude, Longitude, Day, Month, Year, all_of(cnames),
         Date, YearMonth, Week)

typing_data <- lapply(1:length(interval_list), function(i) {
  n1 <- as.character(interval_list[i])
  tpkstrains <- metadata[get(interval) <= n1]$Strain
  dfz <- tpn$filedata %>% rownames_to_column("isolate") %>%
    select(isolate, all_of(hx$h)) %>%
    filter(isolate %in% tpkstrains) %>% column_to_rownames("isolate")
  dfz[,hx$h[1],drop=FALSE] %>% set_colnames(hx$th[1])
}) %>% set_names(as.character(interval_list))

# tp1eccs <- grep("TP1", ecccols, value = TRUE)
# tp2eccs <- grep("TP2", ecccols, value = TRUE)

assert("No clusters with unassigned type", checkTypes(step1))

# TYPE MODIFICATIONS FOR ECCs --------------------------------------------------------------------

# Type I modifications: TP1 > 2, TP2 > 2, TP1 = TP2
#   - None required, we can use the ECC stats for TP1 & TP2, we can use the cluster averages for TP1 & TP2


# Type II modifications: TP1 > 2, TP2 > 2, TP2 > TP1
#   -	Main problem is that the novel strains in TP2 do not have TP1 data
#   -	TP1, no modification required
#   -	TP2: 
#     - no change for strains also in TP1; 
#     - for novel strains in TP2: 
#       - in TP1, needs to have the cluster size and ECC stats from the TP1 strains they cluster with in TP2, 
#       - need to have the TP1 cluster number

cases2a <- step1[type == "Type2"]
if (nrow(cases2a) > 0) {
  cases2 <- cases2a %>% type2Inheritance(.)
  step1 <- step1[ type != "Type2" ] %>% bind_rows(cases2)
}

# Type III modifications: TP1 < 3, TP2 > 2
# -	Main problem is that TP1 cluster doesn't have ECC stats, 
#   - impacts the change vector calculation; 
#   - also, if TP1 = 0 then cluster size for bubble plot & the cluster growth have no data 
#     - (no bubble for TP1 & -Inf growth rate)
# - TP1 needs to have a size of 1 
#   - (+ 1 adjustment for every cluster) so that the denominator is not 0 for cluster growth
#   - (initially wanted ECC bubbles of size at least 1, NOW not doing cluster size increment for the ECCs)

cases3a <- step1 %>% filter(type == "Type3")
if (nrow(cases3a) > 0) {
  cases3b <- type3Inheritance(cases3a)
  step1 <- step1[ type != "Type3" ] %>% bind_rows(cases3b)
}


# Type IV modifications: TP1 < 3, TP2 < 3
# -	Main problem is that TP1 and TP2 are both small and do not have ECC stats 
#   since they are singletons or non-existent
# -	Force TP1 and TP2 ECC stats to blanks
# -	Filter these strains prior to analysis & give ECC blanks
# -	Eventually, include in analysis but maybe do not include them in EpiMatrix calculation

index_eccs <- grep("ECC", colnames(step1))

# should already be NA
assert("All TP1 and TP2 ECCs for these should be blank (NA, for now)", 
       all(is.na(step1[ type == "Type4", ..index_eccs])))

# step1 %<>% mutate(across(all_of(index_eccs), as.character))
cases4 <- step1 %>% filter(type == "Type4")

assert("All singletons or nonexistent (at both TP1 and TP2)", 
       all(c(cases4$tp1_cl_size, cases4$tp2_cl_size) < 3))
# for (index_j in index_eccs) {set(cases4, j = index_j, value = "")}
# step1 <- step1[ type != "Type4"] %>% bind_rows(cases4)

# check that the only NA ECCs are for Type4 cases now:
ecccols <- grep("ECC", colnames(step1), value = TRUE) %>% sort(decreasing = TRUE)
assert("Only blank ECCs are for Type4 cases", identical(unique(step1[is.na(get(ecccols)),type]), "Type4"))

# END OF TYPE MODIFICATIONS FOR ECCs --------------------------------------------------------------

# adding basic delta ECC columns
step2 <- step1 %>% as_tibble()
for (i in 1:(length(ecccols)/2)) {
  a <- tp1eccs[i] %>% strsplit(., split = "ECC") %>% unlist() %>% extract2(2) %>% 
    substr(., 2, nchar(.)) %>% paste0("delta_ECC_", .)
  z <- pull(step2, tp2eccs[i]) - pull(step2, tp1eccs[i])
  step2[a] <- z
}

# adding average lat and long columns
step3 <- step2 %>% 
  left_join(., getAverage(step2, tp1_cl, Latitude, "avg_lat_1"), by = "tp1_cl") %>% 
  left_join(., getAverage(step2, tp1_cl, Longitude, "avg_long_1"), by = "tp1_cl") %>% 
  left_join(., getAverage(step2, tp2_cl, Latitude, "avg_lat_2"), by = "tp2_cl") %>% 
  left_join(., getAverage(step2, tp2_cl, Longitude, "avg_long_2"), by = "tp2_cl")

# adding average date columns
step4 <- step3 %>% 
  left_join(., getAverage(step3, tp1_cl, Date, "avg_date1"), by = "tp1_cl") %>% 
  left_join(., getAverage(step3, tp2_cl, Date, "avg_date2"), by = "tp2_cl")

dist_avgs <- grep("avg", colnames(step4), value = TRUE) %>% grep("dist", ., value = TRUE) %>% sort()

time1 <- readRDS(arg$tp1)$lookup_table %>% 
  select(1, as.double(filtering_ags$th[2])+2) %>% 
  set_colnames(c("Strain", "Actual TP1 cluster"))

time2 <- readRDS(arg$tp2)$lookup_table %>% 
  select(1, as.double(filtering_ags$th[2])+2) %>% 
  set_colnames(c("Strain", "Actual TP2 cluster"))

step5 <- step4 %>% 
  left_join(., time1, by = intersect(colnames(.), colnames(time1))) %>% 
  left_join(., time2, by = intersect(colnames(.), colnames(time2)))

step6 <- step5 %>% 
  select(Strain, Country, Province, City, Latitude, Longitude, Day, Month, Year, found_in_TP1, 
         `Actual TP1 cluster`, tp1_cl, TP1_T0_Size, all_of(tp1eccs), found_in_TP2, 
         `Actual TP2 cluster`, tp2_cl, TP2_T0_Size, all_of(tp2eccs), 
         grep("delta_ECC", colnames(step5), value = TRUE), avg_date1, 
         getDistCols(dist_avgs, "TP1", "temp", TRUE), avg_lat_1, avg_long_1, 
         getDistCols(dist_avgs, "TP1", "geog", TRUE), avg_date2, 
         getDistCols(dist_avgs, "TP2", "temp", TRUE), avg_lat_2, avg_long_2, 
         getDistCols(dist_avgs, "TP2", "geog", TRUE), first_tp1_flag, last_tp1_flag, 
         first_tp2_flag, last_tp2_flag, tp1_cl_size, tp2_cl_size, 
         actual_size_change, add_TP1, num_novs, actual_growth_rate, new_growth, type)

step7 <- step6 %>% 
  rename_with(., replaceDistName, getDistCols(colnames(.), "TP1", "temp")) %>% 
  rename_with(., replaceDistName, getDistCols(colnames(.), "TP2", "temp")) %>% 
  rename_with(., replaceDistName, getDistCols(colnames(.), "TP1", "geo")) %>% 
  rename_with(., replaceDistName, getDistCols(colnames(.), "TP2", "geo"))

step8 <- step7 %>% 
  rename("TP1 cluster" = tp1_cl, "TP1 cluster size (1)" = TP1_T0_Size, 
         "TP2 cluster" = tp2_cl, "TP2 cluster size (1)" = TP2_T0_Size, 
         "TP1" = found_in_TP1, "TP2" = found_in_TP2, 
         "Average TP1 date" = avg_date1, "Average TP1 latitude"	= avg_lat_1, 
         "Average TP1 longitude" = avg_long_1, 
         "Average TP2 date" = avg_date2, "Average TP2 latitude" = avg_lat_2, 
         "Average TP2 longitude" = avg_long_2, 
         "First time this cluster was seen in TP1" = first_tp1_flag, 
         "Last time this cluster was seen in TP1" = last_tp1_flag, 
         "First time this cluster was seen in TP2" = first_tp2_flag, 
         "Last time this cluster was seen in TP2" = last_tp2_flag, 
         "TP1 cluster size + 1 (2)"	= tp1_cl_size, "TP2 cluster size + 1 (2)"	= tp2_cl_size, 
         "Actual cluster size (TP2 size - TP1 size)" = actual_size_change, 
         "Number of additional TP1 strains in the TP2 match" = add_TP1, 
         "Number of novels in the TP2 match" = num_novs, 
         "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)" = actual_growth_rate, 
         "Novel growth = (TP2 size) / (TP2 size - number of novels)" = new_growth, 
         "Type" = type) %>% 
  arrange(`TP2 cluster`, `TP1 cluster`, Strain)

writeData(fp = "results/Merged_strain_results.tsv", df = step8)

step8 %>% 
  group_by(`TP2 cluster`) %>% slice(1) %>% 
  select(-Strain) %>% ungroup() %>% 
  writeData(fp = "results/Merged_cluster_results.tsv", df = .)

cat(paste0("See 'results' folder for cluster-specific and strain-specific files.\n"))
cat(paste0("\n||", paste0(rep("-", 35), collapse = ""), " End of merging step ", 
           paste0(rep("-", 35), collapse = ""), "||\n"))
