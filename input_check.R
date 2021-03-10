
x <- c("reshape2", "dplyr", "tidyr", "readr", "stringr", "magrittr", "tibble", "purrr", "readxl")
lapply(x, require, character.only = TRUE)

# Reading and processing the data in the full metadata excel file ----------------------------------------------
# 6000 rows
file1 <- "Europe TP2 metadata.xlsx"
strain_data <- read_excel(file1, sheet = 4, col_names = TRUE, .name_repair = "minimal") %>% 
  set_colnames(c("Strain", "Source", "Country", "Province", "City", "Latitude", "Longitude", "Day", 
                 "Month", "Year", "TP1", "TP1_T0", "TP1_T4", "TP1_T5", "TP1_T10", "TP1_T21", "TP1_T187", 
                 "TP1_Order", "TP1_Size_T0", "TP1_Size_T4", "TP1_Size_T5", "TP1_Size_T10", "TP1_Size_11", 
                 "TP1_Size_T187", "TP2", "TP2_T0", "TP2_9", "TP2_10", "TP2_14", "TP2_19", "TP2_336", 
                 "TP2_Size_T1", "TP2_Size_T9", "TP2_Size_T10", "TP2_Size_T14", "TP2_Size_T19", 
                 "TP2_Size_T336", "TP2_Order")) %>% 
  select(c("Strain", "Source", "Country", "Province", "City", "Latitude", "Longitude", "Day", 
           "Month", "Year", "TP1", "TP1_T0", "TP1_Size_T0", "TP2", "TP2_T0")) %>% 
  filter(Province != "England") %>% filter(Country != "United Kingdom")

# 5628 rows with TP2 == 1, not 0 or NA
# 3565 rows
tp2_strain_data <- strain_data %>% filter(TP2 == 1)

# 1808 rows
tp1_strain_data <- strain_data %>% filter(TP2 == 1) %>% filter(TP1 == 1)

# Reading and processing the cluster data for TP1 and TP2, making sure they match the metadata file ------------
# TP1 DATA -----------------------------------------------------------------------------------------------------
timepoint1 <- read.csv(file = "european-t1_clusters.csv", stringsAsFactors = FALSE, 
                       numerals = "no.loss", check.names = FALSE, sep = ",") %>% as_tibble() %>% 
  set_colnames(c("Strain", 0:ncol(.)))

timepoint1$Strain <- timepoint1$Strain %>% gsub("hCoV-19/", "", .) %>% sub("\\|.*", "", .)
timepoint1$Place <- timepoint1$Strain %>% gsub("\\/.*", "", .)

# 1808 rows
timepoint1 <- timepoint1 %>% filter(!(Place %in% c("England", "Scotland", "Wales"))) %>% select(-Place)

# check that all strains are present in both
if (identical(sort(tp1_strain_data$Strain), sort(timepoint1$Strain))) {
  write.table(timepoint1, "t1_clusters_processed.csv", row.names = FALSE, quote = FALSE, sep = "\t")
}


# TP1 DATA -----------------------------------------------------------------------------------------------------
timepoint2 <- read.csv(file = "european-t2_clusters.csv", stringsAsFactors = FALSE,
                       numerals = "no.loss", check.names = FALSE, sep = ",") %>% as_tibble() %>% 
  set_colnames(c("Strain", 0:ncol(.)))

timepoint2$Strain <- timepoint2$Strain %>% gsub("hCoV-19/", "", .) %>% sub("\\|.*", "", .)
timepoint2$Place <- timepoint2$Strain %>% gsub("\\/.*", "", .)

# 3565 rows
timepoint2 <- timepoint2 %>% filter(!(Place %in% c("England", "Scotland", "Wales"))) %>% select(-Place)

# check that ll strains are present in both
if (identical(sort(tp2_strain_data$Strain), sort(timepoint2$Strain))) {
  write.table(timepoint2, "t2_clusters_processed.csv", row.names = FALSE, quote = FALSE, sep = "\t")
}

# Extracting clusters for T0 from both TP1 and TP2, to run EpiQuant on -----------------------------------------
tp1h0 <- timepoint1 %>% select(Strain, '0') %>% arrange(Strain) %>% rename(T0 = '0')
write.table(tp1h0, "input_data/TP1_T0.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

tp2h0 <- timepoint2 %>% select(Strain, '0') %>% arrange(Strain) %>% rename(T0 = '0')
write.table(tp2h0, "input_data/TP2_T0.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)





