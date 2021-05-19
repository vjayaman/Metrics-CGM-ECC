#! /usr/bin/env Rscript

libs <- c("optparse","magrittr","tibble", "dplyr", "readr", "testit", "data.table")
y <- suppressMessages(lapply(libs, require, character.only = TRUE))

option_list <- list(
  make_option(c("-e", "--ECCs"), metavar = "file", default = "results/ECCs.tsv", help = "ECC result file"),
  make_option(c("-c", "--CGMs"), metavar = "file", default = "results/CGM_strain_results.tsv", help = "CGM result file"),
  make_option(c("-s", "--strains"), metavar = "file", default = "inputs/strain_info.txt", help = "Strain metadata file"))

arg <- parse_args(OptionParser(option_list=option_list))

# repNA <- function(vec_x, i) {ifelse(is.na(vec_x), i, vec_x)}

writeData <- function(fp, df) {
  write.table(df, fp, row.names = FALSE, quote = FALSE, sep = "\t")
}

checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

readData <- function(fp) {
  read.table(fp, stringsAsFactors = FALSE, header = TRUE, fileEncoding = checkEncoding(fp)) %>% 
    as.data.table() %>% return()
}

# inheritCol <- function(vec) {
#   x <- vec %>% pull() %>% na.omit() %>% unique()
#   assert("There are more than one ECC assigned to this cluster - incorrect merge!", length(x) == 1)
#   vec <- x
#   return(vec)
# }

checkTypes <- function(df) {
  # verify types
  x <- df %>% select(tp1_cl_size, tp2_cl_size, type) %>% as_tibble() %>% 
    rownames_to_column("id") %>% mutate(across(id, as.integer))
  
  x1 <- x %>% filter(tp1_cl_size > 2, tp2_cl_size > 2, tp1_cl_size == tp2_cl_size)
  x2 <- x %>% filter(tp1_cl_size > 2, tp2_cl_size > 2, tp2_cl_size > tp1_cl_size)
  x3 <- x %>% filter(tp1_cl_size < 3, tp2_cl_size > 2)
  x4 <- x %>% filter(tp1_cl_size < 3, tp2_cl_size < 3)
  
  assert("Assigned types are correct", all(
    unique(x1$type) == "Type1", 
    unique(x2$type) == "Type2", 
    unique(x3$type) == "Type3", 
    unique(x4$type) == "Type4"))
  
  return(identical(sort(c(x1$id, x2$id, x3$id, x4$id)), x$id))
}

getAverage <- function(df, grp, avg, newname) {
  df %>% group_by({{grp}}) %>% summarise(avg_col = mean({{avg}})) %>% 
    rename_at(2, ~newname)
}

getDistCols <- function(cnames, x, y, returnVal = TRUE) {
  grep(paste0("(?=.*", x, ")(?=.*", y, ")"), cnames, value = returnVal, perl = TRUE) %>% return()
}

replaceDistName <- function(x) {
  a1 <- ifelse(grepl("TP1", x), "TP1", "TP2")
  if (grepl("temp", x)) {
    paste0(a1, " temp average cluster distance (days)") %>% return()
  }else {
    paste0(a1, " geo average cluster distance (km)") %>% return()
  }
}

colsToInherit <- function(df) {
  grep("tp1|TP1", colnames(df), value = TRUE) %>% 
    grep("found_in", ., value = TRUE, invert = TRUE) %>% return()
}

inheritTP1Data <- function(typeX, clustersX) {
  cnames <- colsToInherit(typeX)
  
  for (k in 1:length(clustersX)) {
    to_inherit <- typeX[first_tp2_flag == clustersX[k] & !is.na(tp1_id), ..cnames] %>% unique()
    
    if (nrow(to_inherit) == 1) {
      rowsY <- typeX[   first_tp2_flag == clustersX[k] & is.na(tp1_id), -(..cnames)] %>% cbind(to_inherit)
      typeX <- typeX[ !(first_tp2_flag == clustersX[k] & is.na(tp1_id)) ] %>% bind_rows(., rowsY)
      
    }else {
      stop("Number of TP1 clusters to inherit from is 0 or > 1")
    }  
  }
  return(typeX)
}

type2Inheritance <- function(typeX) {
  clustersX <- typeX$first_tp2_flag %>% unique()
  inheritTP1Data(typeX, clustersX) %>% return()
}

# for Type 3, 
#   - in some cases we have TP1 singletons to inherit cluster info from (but no ECC), TP1 ECC := 1
#   - and in other cases no TP1 cluster to inherit from --> in this case, TP1 ECC := 1
type3Inheritance <- function(typeX) {
  
  x1 <- getDistCols(colnames(typeX), "TP1", "ECC", FALSE)
  assert("All TP1 clusters have ECCs of NA", all(is.na(typeX[,..x1])))
  
  set(typeX, j = x1, value = 1)
  
  # for those that have NA tp1_id, we know they do not exist at TP1 for Type 3, so we just set 
  # the TP1 ECC to 1, there is no TP1 cluster info to inherit 
  clustersX <- typeX[!is.na(tp1_id)] %>% pull(first_tp2_flag) %>% unique()
  
  # if 0, cluster k is entirely novel, so no TP1 data to inherit, we just force ECCs to be 1
  inheritTP1Data(typeX, clustersX) %>% return()
}

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), " Merging CGM and ECC results ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))

# ------------------------------------------------------------------------------------------------------------
# NOW SAVING OUTPUTS AND MERGING ECCS WITH CGM DATA ----------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
eccs <- readData(arg$ECCs)
cgms <- readData(arg$CGMs)

# actually assigned a cluster at TP2, not NA (just in case we don't have typing data for some strains)
strain_data <- suppressMessages(read_tsv(arg$strains)) %>% 
  filter(TP2 == 1) %>% 
  mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
  select(Strain, Source, City, Province, Country, Latitude, Longitude, Date, Day, Month, Year, TP1, TP2) %>% 
  as.data.table()

step1 <- merge.data.table(cgms, eccs, by = "Strain") %>% select(-TP1, -TP2) %>% 
  merge.data.table(., strain_data, by = "Strain") %>% 
  rename(found_in_TP1 = TP1, found_in_TP2 = TP2)

tp1cnames <- grep("TP1", colnames(step1), value = TRUE)
tp2cnames <- grep("TP2", colnames(step1), value = TRUE)

ecccols <- grep("ECC", colnames(step1), value = TRUE) %>% sort(decreasing = TRUE)
tp1eccs <- grep("TP1", ecccols, value = TRUE)
tp2eccs <- grep("TP2", ecccols, value = TRUE)

assert("No clusters with unassigned type", checkTypes(step1))


# TYPE MODIFICATIONS FOR ECCs --------------------------------------------------------------------

# dim(step1): 35,627 by 41

# Type I modifications: TP1 > 2, TP2 >2, TP1 = TP2
# ○	None required
# ○	we can use the ECC stats for TP1 & TP2, we can use the cluster averages for TP1 & TP2

# cases1 <- step1 %>% filter(type == "Type1")
# nrow(cases1) == 8

# Type II modifications: TP1 > 2, TP2 > 2, TP2 > TP1
# ○	Main problem is that the novel strains in TP2 don’t have TP1 data
# ○	TP1, no modification required
# ○	TP2: 
#   - no change for strains also in TP1; 
#   - for novel strains in TP2: 
#     - in TP1, needs to have the cluster size and ECC stats from the TP1 strains they cluster with in TP2, 
#     - need to have the TP1 cluster number

cases2 <- step1 %>% filter(type == "Type2") %>% type2Inheritance(.)
step1 <- step1[ type != "Type2" ] %>% bind_rows(cases2)


# Type III modifications: TP1 < 3, TP2 > 2
# ○	Main problem is that TP1 cluster doesn’t have ECC stats, 
#   - impacts the change vector calculation; 
#   - also, if TP1 = 0 then cluster size for bubble plot & the cluster growth have no data 
#     - (no bubble for TP1 & “Inf” growth rate)
# ○	TP1 needs to have a size of 1 
#   - (+ 1 adjustment for every cluster) so that the denominator is not 0 for cluster growth
#   - (initially wanted ECC bubbles of size at least 1, NOW not doing cluster size increment for the ECCs)


cases3 <- step1 %>% filter(type == "Type3") %>% type3Inheritance()
step1 <- step1[ type != "Type3" ] %>% bind_rows(cases3)


# Type IV modifications: TP1 < 3, TP2 < 3
# ○	Main problem is that TP1 and TP2 are both small and do not have ECC stats 
#   since they are singletons or non-existent
# ○	Force TP1 and TP2 ECC stats to blanks
# ○	Filter these strains prior to analysis & give ECC blanks
# ○	Eventually, include in analysis but maybe do not include them in EpiMatrix calculation

cases4 <- step1 %>% filter(type == "Type4")
assert("All singletons or nonexistent (at both TP1 and TP2)", 
       all(c(cases4$tp1_cl_size, cases4$tp2_cl_size) < 3))

for (index_j in grep("ECC", colnames(cases4))) {
  set(cases4, j = index_j, value = 0)  
}
step1 <- step1[ type != "Type4"] %>% bind_rows(cases4)

# END OFTYPE MODIFICATIONS FOR ECCs --------------------------------------------------------------



# adding basic delta ECC columns
step2a <- step1 %>% as_tibble()
for (i in 1:(length(ecccols)/2)) {
  a <- tp1eccs[i] %>% strsplit(., split = "ECC") %>% unlist() %>% extract2(2) %>% 
    substr(., 2, nchar(.)) %>% paste0("delta_ECC_", .)
  z <- pull(step2a, tp2eccs[i]) - pull(step2a, tp1eccs[i])
  step2a[a] <- z
}

# adding average lat and long columns
step3 <- step2a %>% 
  left_join(., getAverage(step2a, tp1_cl, Latitude, "avg_lat_1"), by = "tp1_cl") %>% 
  left_join(., getAverage(step2a, tp1_cl, Longitude, "avg_long_1"), by = "tp1_cl") %>% 
  left_join(., getAverage(step2a, tp2_cl, Latitude, "avg_lat_2"), by = "tp2_cl") %>% 
  left_join(., getAverage(step2a, tp2_cl, Longitude, "avg_long_2"), by = "tp2_cl")

# adding average date columns
step4 <- step3 %>% 
  left_join(., getAverage(step3, tp1_cl, Date, "avg_date1"), by = "tp1_cl") %>% 
  left_join(., getAverage(step3, tp2_cl, Date, "avg_date_2"), by = "tp2_cl")

dist_avgs <- grep("avg", colnames(step4), value = TRUE) %>% grep("dist", ., value = TRUE) %>% sort()

step6 <- step4 %>% 
  select(Strain, Country, Province, City, Latitude, Longitude, Day, Month, Year, 
         found_in_TP1, tp1_cl, TP1_T0_Size, 
         all_of(tp1eccs), 
         found_in_TP2, tp2_cl, TP2_T0_Size, 
         all_of(tp2eccs), 
         grep("delta_ECC", colnames(step4), value = TRUE),
         avg_date1, 
         getDistCols(dist_avgs, "TP1", "temp", TRUE),
         avg_lat_1, avg_long_1, 
         getDistCols(dist_avgs, "TP1", "geog", TRUE),
         avg_date_2, 
         getDistCols(dist_avgs, "TP2", "temp", TRUE),
         avg_lat_2, avg_long_2, 
         getDistCols(dist_avgs, "TP2", "geog", TRUE),
         first_tp1_flag, last_tp1_flag, first_tp2_flag, last_tp2_flag, tp1_cl_size, tp2_cl_size, 
         actual_size_change, add_TP1, num_novs, actual_growth_rate, new_growth, type) %>% 
  rename_with(., replaceDistName, getDistCols(colnames(.), "TP1", "temp")) %>% 
  rename_with(., replaceDistName, getDistCols(colnames(.), "TP2", "temp")) %>% 
  rename_with(., replaceDistName, getDistCols(colnames(.), "TP1", "geo")) %>% 
  rename_with(., replaceDistName, getDistCols(colnames(.), "TP2", "geo"))

step7 <- step6 %>% 
  rename("TP1 cluster" = tp1_cl, 
         "TP1 cluster size (1)" = TP1_T0_Size, 
         "TP2 cluster" = tp2_cl, 
         "TP2 cluster size (1)" = TP2_T0_Size, 
         "Average TP1 date" = avg_date1, 
         "Average TP1 latitude"	= avg_lat_1, 
         "Average TP1 longitude" = avg_long_1, 
         "Average TP2 date" = avg_date_2, 
         "Average TP2 latitude" = avg_lat_2, 
         "Average TP2 longitude" = avg_long_2, 
         "First time this cluster was seen in TP1" = first_tp1_flag, 
         "Last time this cluster was seen in TP1" = last_tp1_flag, 
         "First time this cluster was seen in TP2" = first_tp2_flag, 
         "Last time this cluster was seen in TP2" = last_tp2_flag, 
         "TP1 cluster size (2)"	= tp1_cl_size, 
         "TP2 cluster size (2)"	= tp2_cl_size, 
         "Actual cluster size (TP2 size – TP1 size)" = actual_size_change, 
         "Number of additional TP1 strains in the TP2 match" = add_TP1, 
         "Number of novels in the TP2 match" = num_novs, 
         "Actual growth rate = (TP2 size – TP1 size) / (TP1 size)" = actual_growth_rate, 
         "Novel growth = (TP2 size) / (TP2 size – number of novels)" = new_growth) %>% 
    arrange(`TP2 cluster`, `TP1 cluster`, Strain)



# # Clusters that are completely novel at TP2 should have a value of 1
# completely_novel <- which(step1$novel == step1$tp2_cl_size)
# if (length(completely_novel) > 0) {
#   step1[completely_novel, ecccols] <- 1  
# }
#  
# # Mixed TP2 clusters that contain novels and were not singletons at TP1
# stx <- step1 %>% 
#   filter(num_novs > 0 & num_novs != tp2_cl_size & tp1_cl_size != 1) %>% 
#   pull(first_tp2_flag) %>% unique()
# 
# # The novels in these clusters inherit the TP1 ECCs of the originals in their TP2 cluster
# for (x in stx) {
#   inds <- which(step1$first_tp2_flag == x)
#   step1[inds, tp1eccs[1]] %<>% inheritCol()
#   step1[inds, tp1eccs[2]] %<>% inheritCol()
# }
#  
# # Replacing NAs in the ECC columns with 1 - TP1 ECCs for TP1 singletons
# step1[which(step1$tp1_cl_size <= 1), tp1eccs] %<>% apply(., 2, repNA, i = 1.000) %>% as_tibble()
# # Replacing NAs in the ECC columns with 1 - TP2 ECCs for TP2 singletons
# step1[which(step1$tp2_cl_size <= 1), tp2eccs] %<>% apply(., 2, repNA, i = 1.000) %>% as_tibble()
# # If this fails, then there are NA ECCs with a reason I did not consider
# assert("Not all NA ECCs have been accounted for", !any(is.na(step1[,ecccols])))
# 
# # Replacing NAs in the TP1 id column with a blank character
# step1[,grep("tp1_id", colnames(step1))] %<>% apply(., 2, repNA, i = "") %>% as_tibble()



writeData(fp = "results/Merged_strain_results.tsv", df = step7)

step7 %>% 
  group_by(`TP2 cluster`) %>% slice(1) %>% 
  select(-Strain) %>% ungroup() %>% 
  writeData(fp = "results/Merged_cluster_results.tsv", df = .)

cat(paste0("See 'results' folder for cluster-specific and strain-specific files.\n"))
cat(paste0("\n||", paste0(rep("-", 35), collapse = ""), " End of merging step ", 
           paste0(rep("-", 35), collapse = ""), "||\n"))