x <- c("magrittr", "tibble", "dplyr", "readxl", "purrr", "testit", "readr", "optparse")
y <- lapply(x, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

option_list <- list(
  make_option(c("-e", "--ECCs"), metavar = "file", default = "results/ECCs.tsv", help = "ECC result file"),
  make_option(c("-c", "--CGMs"), metavar = "file", default = "results/CGM_strain_results.txt", help = "CGM result file"),
  make_option(c("-s", "--strains"), metavar = "file", default = "inputs/strain_data.tsv", help = "Strain metadata file"))

arg <- parse_args(OptionParser(option_list=option_list))

repNA <- function(vec_x, i) {
  ifelse(is.na(vec_x), i, vec_x)
}

writeData <- function(fp, df) {
  write.table(df, fp, row.names = FALSE, quote = FALSE, sep = "\t")
}

readData <- function(fp) {
  read.table(fp, stringsAsFactors = FALSE, header = TRUE) %>% as_tibble() %>% return()
}
# ------------------------------------------------------------------------------------------------------------
# NOW SAVING OUTPUTS AND MERGING ECCS WITH CGM DATA ----------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
eccs <- readData(arg$ECCs)
cgms <- readData(arg$CGMs) %>% rename(Strain = strain)

strain_data <- read_tsv(arg$strains) %>% 
  filter(Province != "England") %>% filter(Country != "United Kingdom") %>% 
  filter(TP2 == 1) # actually assigned a cluster at TP2, not NA (185 such cases)
# ------------------------------------------------------------------------------------------------------------

step1 <- left_join(strain_data, cgms, by = "Strain") %>% 
  left_join(., eccs, by = c("Strain", "TP1", "TP2"))
step1[,grep("ECC", colnames(step1))] %<>% apply(., 2, repNA, i = 1.000) %>% as_tibble()
step1[,grep("tp1_id", colnames(step1))] %<>% apply(., 2, repNA, i = "") %>% as_tibble()

step2 <- step1 %>% rename("TP1 cluster" = tp1_id) %>% 
  mutate("TP2 cluster" = first_tp2_flag, "TP1 cluster size" = tp1_cl_size, 
         "TP2 cluster size" = tp2_cl_size) %>% 
  select(Strain, Country, Province, City, Latitude, Longitude, Day, Month, Year, 
         TP1, `TP1 cluster`, tp1_cl_size, TP1_T0_ECC_0.0.1, TP1_T0_ECC_0.1.0, 
         TP2, `TP2 cluster`, tp2_cl_size, TP2_T0_ECC_0.0.1, TP2_T0_ECC_0.1.0,
         first_tp1_flag, last_tp1_flag, first_tp2_flag, last_tp2_flag, `TP1 cluster size`, 
         `TP2 cluster size`, actual_size_change, add_TP1, num_novs, actual_growth_rate, new_growth) %>% 
  rename("TP1 cluster size (2)" = tp1_cl_size, 
         "TP2 cluster size (2)" = tp2_cl_size, 
         "First time this cluster was seen in TP1" = first_tp1_flag, 
         "Last time this cluster was seen in TP1" = last_tp1_flag, 
         "First time this cluster was seen in TP2" = first_tp2_flag, 
         "Last time this cluster was seen in TP2" = last_tp2_flag, 
         "Actual cluster size change (TP2 size - TP1 size)" = actual_size_change, 
         "Number of additional TP1 strains in the TP2 match" = add_TP1, 
         "Number of novels in the TP2 match" = num_novs, 
         "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)" = actual_growth_rate, 
         "Novel growth = (TP2 size) / (TP2 size - number of novels)" = new_growth)

writeData(fp = "results/Strain_results.tsv", df = step2)


step2 %>% arrange(`TP2 cluster`) %>% group_by(`TP2 cluster`) %>% slice(1) %>% 
  select(-Strain) %>% ungroup() %>% 
  writeData(fp = "results/Cluster_results.tsv", df = .)


# We realized that the easiest way to get actual numbers for the temporal window and the geographic radius 
# was to use the equations themselves. 
# And it turns out that the first steps in the calculations are to generate temporal and geo matrices already. 
# And then the data in these matrices is used to generate the ECC matrix
# So we figure we may as well also report the cluster averages for geo and temporal
# So we don't have to go back to a calibration curve to figure it out later
