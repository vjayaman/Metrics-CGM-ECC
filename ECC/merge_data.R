x <- c("magrittr", "tibble", "dplyr", "readxl", "purrr")
lapply(x, require, character.only = TRUE)

repNA <- function(vec_x) {
  ifelse(is.na(vec_x), 1.000, vec_x)
}

# ------------------------------------------------------------------------------------------------------------
# NOW SAVING OUTPUTS AND MERGING WITH CGM DATA ---------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
strain_file <- "Europe TP2 metadata.xlsx"
fdir <- "results/adjusted/"
ip <- "input_data/"
cem_res <- "../../Project-CGM/ClusterGrowthMetrics/outputs/CGM_strain_results.txt"

ecc_files <- list.files(fdir) %>% 
  lapply(., function(f) {
    f_parta <- strsplit(f, "_") %>% unlist()
    tptype <- f_parta[1]
    f_partb <- unlist(strsplit(f_parta[2], "\\("))
    ecc_thresh <- f_partb[1]
    ecc_params <- f_partb[2] %>% strsplit(., "\\)") %>% unlist() %>% extract2(1)
    
    read.table(file.path(fdir,f), header = TRUE, stringsAsFactors = FALSE) %>% 
      as_tibble() %>% select(-cut, -W_ECC) %>% 
      set_colnames(c(paste0(tptype, "_", ecc_thresh), paste0(tptype, "_Size_", ecc_thresh), 
                     paste0(tptype, "_ECC_", ecc_params)))
  }) %>% set_names(list.files(fdir))

#  merge ECC and strain data
strain_data <- read_excel(strain_file, sheet = 4, col_names = TRUE, .name_repair = "minimal") %>% 
  set_colnames(c("Strain", "Source", "Country", "Province", "City", "Latitude", "Longitude", "Day", 
                 "Month", "Year", "TP1", "TP1_T0", "TP1_T4", "TP1_T5", "TP1_T10", "TP1_T21", "TP1_T187", 
                 "TP1_Order", "TP1_Size_T0", "TP1_Size_T4", "TP1_Size_T5", "TP1_Size_T10", "TP1_Size_11", 
                 "TP1_Size_T187", "TP2", "TP2_T0", "TP2_9", "TP2_10", "TP2_14", "TP2_19", "TP2_336", 
                 "TP2_Size_T1", "TP2_Size_T9", "TP2_Size_T10", "TP2_Size_T14", "TP2_Size_T19", 
                 "TP2_Size_T336", "TP2_Order")) %>% 
  select(c("Strain", "Source", "Country", "Province", "City", "Latitude", "Longitude", "Day", 
           "Month", "Year", "TP1", "TP1_T0", "TP2", "TP2_T0")) %>% 
  filter(Province != "England") %>% filter(Country != "United Kingdom") %>% 
  filter(TP2 == 1) # actually assigned a cluster at TP2, not NA (185 such cases)
# ------------------------------------------------------------------------------------------------------------

typing_data_files <- list.files(ip, full.names = TRUE)
names(typing_data_files) <- typing_data_files %>% basename() %>% tools::file_path_sans_ext()

read_and_convert <- function(path, cnames) {
  read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, 
             quote = "", stringsAsFactors = FALSE) %>% 
    as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% 
    set_colnames(cnames)
}

typing_data <- lapply(1:length(typing_data_files), function(i) {
  read_and_convert(typing_data_files[i], c("Strain", names(typing_data_files[i])))
}) %>% set_names(names(typing_data_files))

# ------------------------------------------------------------------------------------------------------------
# Note that typing_data[[1]]$TP1_T0 %>% unique() %>% length() is the same as the number of clusters in geo_tp1
# i.e. some of the strains are in the same cluster, as expected
tp1_collected <- Reduce(function(y, z) {full_join(y, z)}, ecc_files[grep("tp1|TP1", names(ecc_files))]) %>% 
  left_join(typing_data[[1]], ., by = "TP1_T0")

# as with TP1, typing_data[[2]]$TP2_T0 %>% unique() %>% length() is the same as the number of clusters in geo_tp2
tp2_collected <- Reduce(function(y, z) {full_join(y, z)}, ecc_files[grep("tp2|TP2", names(ecc_files))]) %>% 
  left_join(typing_data[[2]], ., by = "TP2_T0")

# ------------------------------------------------------------------------------------------------------------

step1 <- right_join(tp1_collected, tp2_collected, by = "Strain") %>% 
  mutate(TP1 = ifelse(is.na(TP1_T0), 0, 1), TP2 = ifelse(is.na(TP2_T0), 0, 1)) %>% 
  left_join(strain_data, ., by = c("Strain", "TP1", "TP1_T0", "TP2", "TP2_T0")) %>% 
  select(1:10, TP1, TP1_T0, TP1_Size_T0, TP1_ECC_0.0.1, TP1_ECC_0.1.0, 
               TP2, TP2_T0, TP2_Size_T0, TP2_ECC_0.0.1, TP2_ECC_0.1.0)

cem_st_results <- read.table(cem_res, sep = "\t", stringsAsFactors = FALSE, header = TRUE) %>% as_tibble() %>% 
  set_colnames(c("Strain", "novel", "first_tp2_flag", "tp2_h", "tp2_cl", "tp2_cl_size", "last_tp2_flag", "tp1_id", 
                 "tp1_h", "tp1_cl", "tp1_cl_size", "first_tp1_flag", "last_tp1_flag", "add_tp1", "num_novs", 
                 "actual_size_change", "actual_growth", "novel_growth"))

strain_results <- left_join(step1, cem_st_results, by = "Strain") %>% 
  select(Strain, Country, Province, City, Latitude, Longitude, Day, Month, Year, TP1, first_tp1_flag, 
         tp1_cl_size, TP1_ECC_0.0.1, TP1_ECC_0.1.0, TP2, first_tp2_flag, tp2_cl_size, TP2_ECC_0.0.1, 
         TP2_ECC_0.1.0, last_tp1_flag, last_tp2_flag, actual_size_change, add_tp1, num_novs, 
         actual_growth, novel_growth) %>% 
  rename(TP1_cluster = first_tp1_flag, TP2_cluster = first_tp2_flag) %>% 
  mutate(
    first_tp1_flag = TP1_cluster, first_tp2_flag = TP2_cluster, 
    tp1_cluster_size = tp1_cl_size, tp2_cluster_size = tp2_cl_size)

strain_results[,c(13,14,18,19)] %<>% apply(., 2, repNA) %>% as_tibble()

strain_results %>% set_colnames(c(
  colnames(strain_results)[1:10], "TP1 cluster", "TP1 cluster size", colnames(strain_results)[13:15],
  "TP2 cluster", "TP2 cluster size", colnames(strain_results)[18:19], "First time this cluster was seen in TP1",
  "Last time this cluster was seen in TP1", "First time this cluster was seen in TP2",
  "Last time this cluster was seen in TP2", "TP1 cluster size", "TP2 cluster size",
  "Actual cluster size change (TP2 size - TP1 size)", "Number of additional TP1 strains in TP2 match",
  "Number of novels in TP2 match", "Actual growth = (TP2 size - TP1 size) / (TP1 size)",
  "Novel growth = (TP2 size) / (TP2 size - number of novels)"
)) %>%
  write.table(., "processed_results/strain_data.txt", row.names = FALSE, quote = FALSE, sep = "\t")

