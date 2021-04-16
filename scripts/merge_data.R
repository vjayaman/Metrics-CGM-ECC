#! /usr/bin/env Rscript

libs <- c("optparse","magrittr","tibble", "dplyr", "readr", "testit")
y <- suppressMessages(lapply(libs, require, character.only = TRUE))

option_list <- list(
  make_option(c("-e", "--ECCs"), metavar = "file", default = NULL, help = "ECC result file"),
  make_option(c("-c", "--CGMs"), metavar = "file", default = NULL, help = "CGM result file"),
  make_option(c("-s", "--strains"), metavar = "file", default = NULL, help = "Strain metadata file"))

arg <- parse_args(OptionParser(option_list=option_list))

repNA <- function(vec_x, i) {
  ifelse(is.na(vec_x), i, vec_x)
}

writeData <- function(fp, df) {
  write.table(df, fp, row.names = FALSE, quote = FALSE, sep = "\t")
}

checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

readData <- function(fp) {
  read.table(fp, stringsAsFactors = FALSE, header = TRUE, fileEncoding = checkEncoding(fp)) %>% as_tibble() %>% return()
}

inheritCol <- function(vec) {
  x <- vec %>% pull() %>% na.omit() %>% unique()
  assert("There are more than one ECC assigned to this cluster - incorrect merge!", length(x) == 1)
  vec <- x
  return(vec)
}

cat(paste0("\n||", paste0(rep("-", 31), collapse = ""), " Merging CGM and ECC results ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))

# ------------------------------------------------------------------------------------------------------------
# NOW SAVING OUTPUTS AND MERGING ECCS WITH CGM DATA ----------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
eccs <- readData(arg$ECCs)
cgms <- readData(arg$CGMs)

# actually assigned a cluster at TP2, not NA (185 such cases)
strain_data <- suppressMessages(read_tsv(arg$strains)) %>% 
  filter(TP2 == 1) %>% 
  select(Strain, Source, City, Province, Country, Latitude, Longitude, Day, Month, Year, TP1, TP2)

step1 <- left_join(cgms, eccs, by = "Strain") %>% select(-TP1, -TP2) %>% 
  left_join(., strain_data, by = "Strain")

tp1cnames <- grep("TP1", colnames(step1), value = TRUE)
tp2cnames <- grep("TP2", colnames(step1), value = TRUE)

ecccols <- grep("ECC", colnames(step1), value = TRUE) %>% sort(decreasing = TRUE)
tp1eccs <- grep("TP1", ecccols, value = TRUE)
tp2eccs <- grep("TP2", ecccols, value = TRUE)


# Clusters that are completely novel at TP2 should have a value of 1
step1[which(step1$novel == step1$tp2_cl_size), ecccols] <- 1

# Mixed TP2 clusters that contain novels and were not singletons at TP1
stx <- step1 %>% filter(novel == 1 & novel != tp2_cl_size & tp1_cl_size != 1) %>% 
  pull(first_tp2_flag) %>% unique()

# The novels in these clusters inherit the TP1 ECCs of the originals in their TP2 cluster
for (x in stx) {
  inds <- which(step1$first_tp2_flag == x)
  step1[inds, tp1eccs[1]] %<>% inheritCol()
  step1[inds, tp1eccs[2]] %<>% inheritCol()
}

# Replacing NAs in the ECC columns with 1 - TP1 ECCs for TP1 singletons
step1[which(step1$tp1_cl_size <= 1), tp1eccs] %<>% apply(., 2, repNA, i = 1.000) %>% as_tibble()
# Replacing NAs in the ECC columns with 1 - TP2 ECCs for TP2 singletons
step1[which(step1$tp2_cl_size <= 1), tp2eccs] %<>% apply(., 2, repNA, i = 1.000) %>% as_tibble()
# If this fails, then there are NA ECCs with a reason I did not consider
assert("Not all NA ECCs have been accounted for", !any(is.na(step1[,ecccols])))

# Replacing NAs in the TP1 id column with a blank character
step1[,grep("tp1_id", colnames(step1))] %<>% apply(., 2, repNA, i = "") %>% as_tibble()

for (coeff in unique(substr(ecccols, 12, 16))) {
  ecc_coeff <- grep(coeff, ecccols, value = TRUE)
  col2 <- grep("TP2", ecc_coeff, value = TRUE)
  col1 <- grep("TP1", ecc_coeff, value = TRUE)
  
  step1 <- step1 %>% mutate(delta = step1 %>% pull(col2) - step1 %>% pull(col1))
  colnames(step1)[which(colnames(step1) == "delta")] <- paste0("delta_ECC_", coeff)
}
  
step2 <- step1 %>% rename("TP1 cluster" = tp1_id) %>% 
  mutate("TP2 cluster" = first_tp2_flag, "TP1 cluster size" = tp1_cl_size, "TP2 cluster size" = tp2_cl_size) %>% 
  
  select(Strain, Country, Province, City, Latitude, Longitude, Day, Month, Year, 
         
         TP1, `TP1 cluster`, tp1_cl_size, all_of(tp1eccs), 
         TP2, `TP2 cluster`, tp2_cl_size, all_of(tp2eccs), 
         
         grep("delta", colnames(step1), value = TRUE), 
         grep("TP1_avg", colnames(step1), value = TRUE), 
         grep("TP2_avg", colnames(step1), value = TRUE), 
         
         first_tp1_flag, last_tp1_flag, first_tp2_flag, last_tp2_flag, `TP1 cluster size`, 
         `TP2 cluster size`, actual_size_change, add_TP1, num_novs, actual_growth_rate, new_growth) %>% 
  
  rename("TP1 cluster size (2)" = tp1_cl_size, 
         "TP2 cluster size (2)" = tp2_cl_size, 
         "TP1 temp average cluster distance (days)" = TP1_avg_temp_dist, 
         "TP1 geo average cluster distance (km)" = TP1_avg_geo_dist, 
         "TP2 temp average cluster distance (days)" = TP2_avg_temp_dist, 
         "TP2 geo average cluster distance (km)" = TP2_avg_geo_dist, 
         "Average TP1 date" = grep("avg_date", tp1cnames, value = TRUE), 
         "Average TP2 date" = grep("avg_date", tp2cnames, value = TRUE), 
         "Average TP1 longitude" = grep("avg_long", tp1cnames, value = TRUE), 
         "Average TP2 longitude" = grep("avg_long", tp2cnames, value = TRUE), 
         "Average TP1 latitude" = grep("avg_lat", tp1cnames, value = TRUE), 
         "Average TP2 latitude" = grep("avg_lat", tp2cnames, value = TRUE), 
         "First time this cluster was seen in TP1" = first_tp1_flag, 
         "Last time this cluster was seen in TP1" = last_tp1_flag, 
         "First time this cluster was seen in TP2" = first_tp2_flag, 
         "Last time this cluster was seen in TP2" = last_tp2_flag, 
         "Actual cluster size change (TP2 size - TP1 size)" = actual_size_change, 
         "Number of additional TP1 strains in the TP2 match" = add_TP1, 
         "Number of novels in the TP2 match" = num_novs, 
         "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)" = actual_growth_rate, 
         "Novel growth = (TP2 size) / (TP2 size - number of novels)" = new_growth) %>% 
  arrange(`TP2 cluster`, `TP1 cluster`, Strain)

writeData(fp = "results/Merged_strain_results.tsv", df = step2)

step2 %>% 
  group_by(`TP2 cluster`) %>% slice(1) %>% 
  select(-Strain) %>% ungroup() %>% 
  writeData(fp = "results/Merged_cluster_results.tsv", df = .)

cat(paste0("See 'results' folder for cluster-specific and strain-specific files.\n"))
cat(paste0("\n||", paste0(rep("-", 35), collapse = ""), " End of merging step ", 
           paste0(rep("-", 35), collapse = ""), "||\n"))