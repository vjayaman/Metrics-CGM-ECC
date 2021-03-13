x <- c("magrittr", "tibble", "dplyr", "readxl", "purrr", "testit", "readr")
y <- lapply(x, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

repNA <- function(vec_x) {
  ifelse(is.na(vec_x), 1.000, vec_x)
}

# ------------------------------------------------------------------------------------------------------------
# NOW SAVING OUTPUTS AND MERGING WITH CGM DATA ---------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
eccs <- read.table("results/ECCs.tsv", stringsAsFactors = FALSE, header = TRUE) %>% as_tibble()
cgms <- read.table("results/CGM_strain_results.txt", sep = "\t", 
                   stringsAsFactors = FALSE, header = TRUE) %>% as_tibble() %>% rename(Strain = strain)

strain_data <- read_tsv("inputs/strain_data.tsv") %>% 
  filter(Province != "England") %>% filter(Country != "United Kingdom") %>% 
  filter(TP2 == 1) # actually assigned a cluster at TP2, not NA (185 such cases)
# ------------------------------------------------------------------------------------------------------------

# strain_results <- left_join(step1, cem_st_results, by = "Strain") %>%
  # select(Strain, Country, Province, City, Latitude, Longitude, Day, Month, Year,
  #        TP1, first_tp1_flag, tp1_cl_size, TP1_ECC_0.0.1, TP1_ECC_0.1.0, 
  #        TP2, first_tp2_flag, tp2_cl_size, TP2_ECC_0.0.1,
  #        TP2_ECC_0.1.0, last_tp1_flag, last_tp2_flag, actual_size_change, add_tp1, num_novs,
  #        actual_growth, novel_growth) %>%
  # rename(TP1_cluster = first_tp1_flag, TP2_cluster = first_tp2_flag) %>%
  # mutate(first_tp1_flag = TP1_cluster, first_tp2_flag = TP2_cluster,
  #   tp1_cluster_size = tp1_cl_size, tp2_cluster_size = tp2_cl_size)

# strain_results[,c(13,14,18,19)] %<>% apply(., 2, repNA) %>% as_tibble()
alldata <- left_join(strain_data, cgms) %>% left_join(., eccs)
step1 <- alldata %>% 
  select(Strain, Country, Province, City, Latitude, Longitude, Day, Month, Year, 
         TP1, tp1_id, tp1_cl_size, TP1_T0_ECC_0.0.1, TP1_T0_ECC_0.1.0, TP2, first_tp2_flag, tp2_cl_size, 
         TP2_T0_ECC_0.0.1, TP2_T0_ECC_0.1.0, first_tp1_flag, last_tp1_flag, last_tp2_flag, tp1_cl_size, 
         tp2_cl_size, actual_size_change, add_TP1, num_novs, actual_growth_rate, new_growth)

# strain_results %>% set_colnames(c(
#   colnames(strain_results)[1:10], "TP1 cluster", "TP1 cluster size", colnames(strain_results)[13:15],
#   "TP2 cluster", "TP2 cluster size", colnames(strain_results)[18:19], "First time this cluster was seen in TP1",
#   "Last time this cluster was seen in TP1", "First time this cluster was seen in TP2",
#   "Last time this cluster was seen in TP2", "TP1 cluster size", "TP2 cluster size",
#   "Actual cluster size change (TP2 size - TP1 size)", "Number of additional TP1 strains in TP2 match",
#   "Number of novels in TP2 match", "Actual growth = (TP2 size - TP1 size) / (TP1 size)",
#   "Novel growth = (TP2 size) / (TP2 size - number of novels)"
# )) %>%
#   write.table(., "processed_results/strain_data.txt", row.names = FALSE, quote = FALSE, sep = "\t")

