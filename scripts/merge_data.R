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

# distCols <- function(x, y, cnames) {
#   grep(paste0("(?=.*", x, ")(?=.*", y, ")"), cnames, perl = TRUE, value = TRUE) %>% return()
# }

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
  merge.data.table(., strain_data, by = "Strain")

tp1cnames <- grep("TP1", colnames(step1), value = TRUE)
tp2cnames <- grep("TP2", colnames(step1), value = TRUE)

ecccols <- grep("ECC", colnames(step1), value = TRUE) %>% sort(decreasing = TRUE)
tp1eccs <- grep("TP1", ecccols, value = TRUE)
tp2eccs <- grep("TP2", ecccols, value = TRUE)

assert("No clusters with unassigned type", checkTypes(step1))

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


# # type modifications for ECCs
# 
# # type 1 - no modifications required
# # type 2
# step5a <- step4 %>% filter(type == "Type2")# %>% filter(novel == 1)
# # type 3
# step3a$TP1_T0_Size <- step3a$tp1_cl_size - 1
# 
# # type 4
# 
# 
# 
dist_avgs <- grep("avg", colnames(step4), value = TRUE) %>% grep("dist", ., value = TRUE) %>% sort()

step6 <- step4 %>% 
  select(Strain, Country, Province, City, Latitude, Longitude, Day, Month, Year, 
         TP1, tp1_cl, TP1_T0_Size, 
         all_of(tp1eccs), 
         TP2, tp2_cl, TP2_T0_Size, 
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