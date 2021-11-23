#! /usr/bin/env Rscript

libs <- c("optparse","magrittr","tibble", "dplyr", "readr", "testit", "data.table", "R6")
y <- suppressMessages(lapply(libs, require, character.only = TRUE))

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", 
              default = "inputs/processed/strain_info.txt", help = "Metadata file"),
  make_option(c("-b", "--tpn"), metavar = "file", 
              default = "inputs/processed/allTP2.Rds", help = "TP2 data"), 
  make_option(c("-f", "--intervalfile"), metavar = "file", 
              default = "inputs/processed/clustersets.Rds"), 
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))

arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)

source("scripts/Misc/type_handling.R"); source("scripts/ECC/classes_ecc.R")

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), " Merging CGM and ECC results ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))

dir.create("results/Merged_strain_results", showWarnings = FALSE)

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
interval_list <- names(clustersets) %>% sort()

# ------------------------------------------------------------------------------------------------------------
# NOW SAVING OUTPUTS AND MERGING ECCS WITH CGM DATA ----------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# For padding height and cluster columns with h0..0.., and c0..0.., respectively
padCol <- function(cvals, padval, padchr) {
  ifelse(!is.na(cvals), formatC(cvals, width = padval, format = "d", flag = "0") %>% 
           paste0(padchr, .), NA) %>% return()
}

# read in CGM results, and use "first_tp2_flag" get the padding values for height and cluster
results_files <- list.files("results/", full.names = TRUE) %>% grep(params$int_type[2], ., value = TRUE)
cgms <- grep("CGM", results_files, value = TRUE) %>% readRDS()

first_ivl <- cgms$interval %>% unique() %>% sort() %>% first()

tmp <- cgms[interval == first_ivl] %>% 
  select(interval, grep("tp1_", colnames(cgms), value = TRUE)) %>% 
  filter(!is.na(tp1_id)) %>% 
  set_colnames(c("interval", "tp2_id", "tp2_cl_size", "first_tp2_flag", "last_tp2_flag"))
tmp2 <- tmp[grep("Absent", tmp$tp2_id, invert = TRUE),] %>% 
  select(-tp2_id) %>% 
  mutate(across(tp2_cl_size, as.double), 
         interval = paste0("set0-", interval_list[1]), 
         first_tp2_flag = gsub("TP1", "TP2", first_tp2_flag), 
         last_tp2_flag = gsub("TP1", "TP2", last_tp2_flag))
cgms <- cgms %>% add_row(tmp2, .before = 1)





a1 <- unlist(strsplit(cgms$first_tp2_flag[1], split = "_"))[2:3] %>% 
  nchar() %>% `-`(1) %>% set_names(c("ph", "pc"))
h_id <- padCol(as.double(hx$h), a1[['ph']], "h")

# read in AVERAGE results and use the TP column to create a TP1_ID column and an identical TP2_ID column: 
# since each TP set is both a TP1 and a TP2, depending on which two timepoints we are comparing
avgs <- grep("AVGS", results_files, value = TRUE) %>% readRDS()
colnames(avgs)[grep("Avg", colnames(avgs))] %<>% paste0(hx$th, "_", .)
avgs <- avgs %>% 
  mutate(tp1_id = paste0("TP1_", h_id, "_", padCol(!!as.symbol(hx$th), a1[['pc']], "c"))) %>% 
  mutate(first_tp2_flag = tp1_id %>% gsub("TP1", "TP2", .))

# read in ECC results and repeat the process
eccs <- grep("ECC", results_files, value = TRUE) %>% readRDS() %>% 
  mutate(tp1_id = paste0("TP1_", h_id, "_", padCol(!!as.symbol(hx$th), a1[['pc']], "c"))) %>% 
  mutate(first_tp2_flag = tp1_id %>% gsub("TP1", "TP2", .))

assert("ECC results have data for all intervals", identical(unique(eccs$TP), interval_list))
assert("AVG results have data for all intervals", identical(unique(avgs$TP), interval_list))

step0 <- merge.data.table(eccs, avgs, by = intersect(colnames(eccs), colnames(avgs)))
stepcols <- ncol(step0)

all_intervals <- c("set0", interval_list)

step1 <- lapply(1:length(all_intervals), function(i) {
  ivl1 <- all_intervals[i]
  ivl2 <- all_intervals[i+1]
  
  partC <- cgms[which(cgms$interval == paste0(c(ivl1, ivl2), collapse = "-"))]
  
  # selecting columns tp1_id, Size, ECC.0.1.0, ECC.0.0.1, Temp.Avg.Dist, Geo.Avg.Dist
  # Avg.Date, Avg.Longitude, Avg.Latitude
  if (ivl1 == "set0") {
    partC$type <- "Type0"
    partA_cnames <- step0 %>% select(4,3,6:all_of(stepcols)) %>% colnames()
    partA <- matrix(ncol = length(partA_cnames)) %>% set_colnames(partA_cnames) %>% 
      set_colnames(gsub(hx$th, "TP1", colnames(.)))
  }else {
    partA <- step0[TP == ivl1] %>% select(4,3,6:all_of(stepcols)) %>% 
      set_colnames(gsub(hx$th, "TP1", colnames(.)))
  }
  
  # selecting columns first_tp2_flag, Size, ECC.0.1.0, ECC.0.0.1, Temp.Avg.Dist, Geo.Avg.Dist, 
  # Avg.Date, Avg.Longitude, Avg.Latitude
  partB <- step0[TP == ivl2] %>% select(5,3,6:all_of(stepcols)) %>% 
    set_colnames(gsub(hx$th, "TP2", colnames(.)))

  partD <- partC %>% 
    merge.data.table(., partA, all.x = TRUE) %>% 
    merge.data.table(., partB, by = "first_tp2_flag")
  
  ecc_cols <- grep("ECC", colnames(partD), value = TRUE)
  partD %>% select(c(interval, tp1_id, tp1_cl_size, TP1_Size, grep("TP1", ecc_cols, value = TRUE),
                     first_tp2_flag, tp2_cl_size, TP2_Size, grep("TP2", ecc_cols, value = TRUE),
                     TP1_Avg.Date, TP1_Temp.Avg.Dist, TP1_Avg.Longitude,TP1_Avg.Latitude, 
                     TP1_Geo.Avg.Dist, TP2_Avg.Date, TP2_Temp.Avg.Dist, TP2_Avg.Longitude, 
                     TP2_Avg.Latitude, TP2_Geo.Avg.Dist, first_tp1_flag, last_tp1_flag, 
                     first_tp2_flag, last_tp2_flag, tp1_cl_size, tp2_cl_size, actual_size_change, 
                     add_TP1, novel, num_novs, actual_growth_rate, new_growth, type))
}) %>% bind_rows()

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

# Type III modifications: TP1 < 3, TP2 > 2
# -	Main problem is that TP1 cluster doesn't have ECC stats, 
#   - impacts the change vector calculation; 
#   - also, if TP1 = 0 then cluster size for bubble plot & the cluster growth have no data 
#     - (no bubble for TP1 & -Inf growth rate)
# - TP1 needs to have a size of 1 
#   - (+ 1 adjustment for every cluster) so that the denominator is not 0 for cluster growth
#   - (initially wanted ECC bubbles of size at least 1, NOW not doing cluster size increment for the ECCs)

# Type IV modifications: TP1 < 3, TP2 < 3
# -	Main problem is that TP1 and TP2 are both small and do not have ECC stats 
#   since they are singletons or non-existent
# -	Force TP1 and TP2 ECC stats to blanks
# -	Filter these strains prior to analysis & give ECC blanks
# -	Eventually, include in analysis but maybe do not include them in EpiMatrix calculation

cat(paste0("Checking type modifications: \n"))
step3 <- data.table()
for (i in 1:length(interval_list)) {
# for (i in 1:(length(interval_list)-1)) {
  cat(paste0("Interval ", all_intervals[i], " - ", all_intervals[i+1], "\n"))
  int_i <- all_intervals[i:(i+1)] %>% paste0(., collapse = "-")
  step2 <- step1[interval == int_i]
  
  cases2a <- step2[type == "Type2"]
  if (nrow(cases2a) > 0) {
    cases2 <- cases2a %>% type2Inheritance(.)
    step2 <- step2[type != "Type2"] %>% bind_rows(cases2)
  }
  
  # the TP1 ECCs should be equal to 1
  cases3a <- step2[type == "Type3"]
  # typeX <- cases3a
  if (nrow(cases3a) > 0) {
    cases3b <- type3Inheritance(cases3a)
    assert("Type 3 inheritances was dealt with during ECC collection", !is.null(cases3b))
    step2 <- step2[type != "Type3"] %>% bind_rows(cases3b)
  }
  
  index_eccs <- grep("ECC", colnames(step1))
  
  # should already be NA
  assert("All TP1 and TP2 ECCs for these should be blank (NA, for now)", 
         all(is.na(step2[ type == "Type4", ..index_eccs])))
  cases4 <- step2 %>% filter(type == "Type4")
  
  assert("All singletons or nonexistent (at both TP1 and TP2)", 
         all(c(cases4$tp1_cl_size, cases4$tp2_cl_size) < 3))
  
  # check that the only NA ECCs are for Type4 cases now:
  ecccols <- grep("ECC", colnames(step2), value = TRUE) %>% sort(decreasing = TRUE)
  
  x1 <- which(is.na(step2[,..ecccols]), arr.ind = TRUE) %>% as.data.frame() %>% pull(row)
  assert("Only blank ECCs are for Type0 or Type4 cases", 
         unique(step2[x1]$type) %in% c("Type0", "Type4"))
  
  step3 <- step3 %>% bind_rows(step2)
}

# END OF TYPE MODIFICATIONS FOR ECCs --------------------------------------------------------------
tp1eccs <- grep("TP1", ecccols, value = TRUE)
tp2eccs <- grep("TP2", ecccols, value = TRUE)

# adding basic delta ECC columns
step4 <- step3 %>% as_tibble()
for (i in 1:(length(ecccols)/2)) {
  a <- tp1eccs[i] %>% strsplit(., split = "ECC") %>% unlist() %>% extract2(2) %>% 
    substr(., 2, nchar(.)) %>% paste0("delta_ECC_", .)
  z <- pull(step4, tp2eccs[i]) - pull(step4, tp1eccs[i])
  step4[a] <- z
}

dist_avgs <- grep("Avg", colnames(step4), value = TRUE) %>% grep("Dist", ., value = TRUE) %>% sort()

step5 <- step4 %>% add_column(tp1_cl = NA, tp2_cl = NA) %>% as.data.table(); rm(step4)
cl_vals <- grep("Absent", step5[!is.na(tp1_id)]$tp1_id, invert = TRUE, value = TRUE)
cl_indices <- which(step5$tp1_id %in% cl_vals)

step5$tp1_cl[cl_indices] <- lapply(step5$tp1_id[cl_indices], function(p) {
  strsplit(p, split = "_c") %>% unlist() %>% extract2(., 2) %>% as.double()
}) %>% unlist()

step5$tp2_cl <- lapply(step5$first_tp2_flag, function(p) {
  strsplit(p, split = "_c") %>% unlist() %>% extract2(., 2) %>% as.double()
}) %>% unlist()

timeN <- readRDS(arg$tpn)
matched_clusters <- timeN$lookup_table[new_h == params$th[2]] %>% select(old_cl, new_cl) %>%
  unique() %>% set_colnames(c("ActCol", "IntCol"))

step6 <- step5 %>% 
  left_join(., matched_clusters, by = c("tp1_cl" = "IntCol")) %>% rename(actual_tp1_cl = ActCol) %>% 
  left_join(., matched_clusters, by = c("tp2_cl" = "IntCol")) %>% rename(actual_tp2_cl = ActCol)

step7 <- step6 %>% add_column(novel_only_ECCs = 0) %>% 
  select(interval, actual_tp1_cl, tp1_cl, TP1_Size, all_of(tp1eccs), 
         actual_tp2_cl, tp2_cl, TP2_Size, all_of(tp2eccs),
         grep("delta_ECC", colnames(step6), value = TRUE), novel_only_ECCs,
         TP1_Avg.Date, TP1_Temp.Avg.Dist, TP1_Avg.Latitude, TP1_Avg.Longitude, TP1_Geo.Avg.Dist, 
         TP2_Avg.Date, TP2_Temp.Avg.Dist, TP2_Avg.Latitude, TP2_Avg.Longitude, TP2_Geo.Avg.Dist,
         first_tp1_flag, last_tp1_flag, first_tp2_flag, last_tp2_flag, tp1_cl_size,
         tp2_cl_size, actual_size_change, add_TP1, num_novs, actual_growth_rate, new_growth, type)

step7[TP2_Size == num_novs]$novel_only_ECCs <- 1
assert("Novel only ECCs indicated with binary value", all(is.na(step7[novel_only_ECCs == 1]$TP1_Size)))

date_order <- lapply(2:length(all_intervals), function(i) {
  paste0(all_intervals[i-1], "-", all_intervals[i])
}) %>% unlist()

step8 <- step7 %>% 
  rename("Actual TP1 cluster" = actual_tp1_cl, 
         "Actual TP2 cluster" = actual_tp2_cl, 
         "TP1 cluster" = tp1_cl, "TP1 cluster size (1)" = TP1_Size, 
         "TP2 cluster" = tp2_cl, "TP2 cluster size (1)" = TP2_Size, 
         "Average TP1 date" = TP1_Avg.Date, "Average TP1 latitude"	= TP1_Avg.Latitude, 
         "Average TP1 longitude" = TP1_Avg.Longitude, 
         "Average TP2 date" = TP2_Avg.Date, "Average TP2 latitude" = TP2_Avg.Latitude, 
         "Average TP2 longitude" = TP2_Avg.Longitude, 
         "Novel only ECCs" = novel_only_ECCs,
         "First time this cluster was seen in TP1" = first_tp1_flag, 
         "Last time this cluster was seen in TP1" = last_tp1_flag, 
         "First time this cluster was seen in TP2" = first_tp2_flag, 
         "Last time this cluster was seen in TP2" = last_tp2_flag, 
         "TP1 cluster size + 1 (2)"	= tp1_cl_size, "TP2 cluster size + 1 (2)"	= tp2_cl_size, 
         "Actual cluster growth (TP2 size - TP1 size)" = actual_size_change, 
         "Number of additional TP1 strains in the TP2 match" = add_TP1, 
         "Number of novels in the TP2 match" = num_novs, 
         "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)" = actual_growth_rate, 
         "Novel growth = (TP2 size) / (TP2 size - number of novels)" = new_growth, 
         "Type" = type) %>% 
  arrange(`TP2 cluster`, `TP1 cluster`) %>% 
  unique() %>% 
  arrange(factor(interval, levels = date_order))

cat(paste0("Saving cluster files to results/ ... \n"))
writeData(fp = "results/Wide_merged_cluster_results.tsv", df = step8)

step9 <- suppressWarnings(melt.data.table(step8, id.vars = c("interval", "Actual TP1 cluster")))
step9 %>% arrange(interval, `Actual TP1 cluster`, variable) %>% 
  writeData(fp = "results/Long_merged_cluster_results.tsv", df = .)

metadata <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()
cnames <- strsplit(params$cnames[2], split = ",") %>% unlist()
strains <- metadata$strain_data %>% 
  select(Strain, all_of(cnames), Latitude, Longitude, Day, Month, Year) %>% 
  mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")))

# ivls <- names(clustersets)
# for (i in 2:length(ivls)) {
#   ivl_i <- paste0(ivls[i-1], "-", ivls[i])
#   step10 <- clustersets[[i]]$ivl %>% select(isolate, heightx) %>% rename(Strain = isolate) %>% 
#     add_column(interval = ivl_i) %>% 
#     inner_join(step8, ., by = c("TP2 cluster" = "heightx", "interval")) %>% 
#     inner_join(., strains, by = "Strain") %>% 
#     select(interval, colnames(strains), colnames(step8))
#   writeData(fp = paste0("results/Merged_strain_results/", ivl_i, ".tsv"), df = step10)
# }

cat(paste0("Saving strain files to results/Merged_strain_results/ ... \n"))
ivls <- all_intervals
for (i in 2:length(ivls)) {
  
  ivl_i <- paste0(ivls[i-1], "-", ivls[i])
  
  step10 <- clustersets[[ivls[i]]]$sofar %>% select(isolate, heightx) %>% 
    rename(Strain = isolate) %>% add_column(interval = ivl_i) %>% 
    inner_join(step8, ., by = c("TP2 cluster" = "heightx", "interval")) %>% 
    inner_join(., strains, by = "Strain") %>% 
    select(interval, colnames(strains), colnames(step8))
  
  writeData(fp = paste0("results/Merged_strain_results/", ivl_i, ".tsv"), df = step10)
  rm(step10)
}

cat(paste0("See 'results' folder for cluster-specific and strain-specific files.\n"))
cat(paste0("\n||", paste0(rep("-", 35), collapse = ""), " End of merging step ", 
           paste0(rep("-", 35), collapse = ""), "||\n"))
