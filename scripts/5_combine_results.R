#! /usr/bin/env Rscript

libs <- c("optparse","magrittr","tibble", "dplyr", "readr", "testit", "data.table")
y <- suppressMessages(lapply(libs, require, character.only = TRUE))

option_list <- list(
  make_option(c("-m", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain metadata file"), 
  make_option(c("-a", "--tp1"), metavar = "file", 
              default = "inputs/processed/allTP1.Rds", help = "TP1 cluster assignments"), 
  make_option(c("-b", "--tp2"), metavar = "file", 
              default = "inputs/processed/allTP2.Rds", help = "TP2 cluster assignments"), 
  make_option(c("-c", "--CGMs"), metavar = "file", default = "results/CGM_strain_results.tsv", help = "CGM result file"),
  make_option(c("-d", "--details"), metavar = "file", 
              default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"), 
  make_option(c("-e", "--ECCs"), metavar = "file", default = "results/ECCs.tsv", help = "ECC result file")
  )

arg <- parse_args(OptionParser(option_list=option_list))

source("scripts/Misc/type_handling.R")

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), " Merging CGM and ECC results ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))

# ------------------------------------------------------------------------------------------------------------
# NOW SAVING OUTPUTS AND MERGING ECCS WITH CGM DATA ----------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
eccs <- readData(arg$ECCs)
cgms <- readData(arg$CGMs)

filtering_params <- readLines(arg$details, warn = FALSE) %>% 
  strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP1","nsTP2","temp_win","cnames"))
cnames <- filtering_params$cnames[2] %>% strsplit(split = ",") %>% unlist()

# actually assigned a cluster at TP2, not NA (just in case we don't have typing data for some strains)
strain_data <- suppressMessages(read_tsv(arg$strains)) %>% 
  mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
  select(Strain, Latitude, Longitude, Day, Month, Year, Date, all_of(cnames), TP1, TP2) %>% 
  as.data.table() %>% 
  filter(TP2 == 1)

step1 <- merge.data.table(cgms, eccs, by = "Strain") %>% 
  merge.data.table(., strain_data, by = "Strain") %>% 
  rename(found_in_TP1 = TP1, found_in_TP2 = TP2)

ecccols <- grep("ECC", colnames(step1), value = TRUE) %>% sort(decreasing = TRUE)
tp1eccs <- grep("TP1", ecccols, value = TRUE)
tp2eccs <- grep("TP2", ecccols, value = TRUE)

assert("No clusters with unassigned type", checkTypes(step1))

# TYPE MODIFICATIONS FOR ECCs --------------------------------------------------------------------

# Type I modifications: TP1 > 2, TP2 >2, TP1 = TP2
#   - None required, we can use the ECC stats for TP1 & TP2, we can use the cluster averages for TP1 & TP2


# Type II modifications: TP1 > 2, TP2 > 2, TP2 > TP1
#   -	Main problem is that the novel strains in TP2 don’t have TP1 data
#   -	TP1, no modification required
#   -	TP2: 
#     - no change for strains also in TP1; 
#     - for novel strains in TP2: 
#       - in TP1, needs to have the cluster size and ECC stats from the TP1 strains they cluster with in TP2, 
#       - need to have the TP1 cluster number

cases2a <- step1 %>% filter(type == "Type2")
if (nrow(cases2a) > 0) {
  cases2 <- cases2a %>% type2Inheritance(.)
  step1 <- step1[ type != "Type2" ] %>% bind_rows(cases2)
}

# Type III modifications: TP1 < 3, TP2 > 2
# -	Main problem is that TP1 cluster doesn't have ECC stats, 
#   - impacts the change vector calculation; 
#   - also, if TP1 = 0 then cluster size for bubble plot & the cluster growth have no data 
#     - (no bubble for TP1 & “Inf” growth rate)
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
  select(1, as.double(filtering_params$th[2])+2) %>% 
  set_colnames(c("Strain", "Actual TP1 cluster"))

time2 <- readRDS(arg$tp2)$lookup_table %>% 
  select(1, as.double(filtering_params$th[2])+2) %>% 
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
