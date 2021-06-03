libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr","reshape2","fossil","tidyr","purrr")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

source("../../Documents/tmp/Metrics-CGM-ECC/ECC/functions/classes_ecc.R")

option_list <- list(
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"), 
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  make_option(c("-x", "--heights"), metavar = "character", default = "0", 
              help = "Comma-delimited string of heights to collect ECCs for"), 
  make_option(c("-t", "--trio"), metavar = "character", default = "010-001", 
              help = "source, temporal, geographic coefficients"))

params <- parse_args(OptionParser(option_list=option_list))

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)

strain_data <- 
  read_tsv(params$strains) %>% 
  mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
         Location = paste(Country, Province, City, sep = "_"))
strain_data <- strain_data[1:9000,]

x1 <- tp1$height_list[[1]]
tp1$height_list[[1]] <- x1[rownames(x1) %in% strain_data$Strain,,drop=FALSE]
tp1$proc %<>% filter(Strain %in% strain_data$Strain)

x2 <- tp2$height_list[[1]]
tp2$height_list[[1]] <- x2[rownames(x2) %in% strain_data$Strain,,drop=FALSE]
tp2$proc %<>% filter(Strain %in% strain_data$Strain)

td <- tp1$height_list %>% append(tp2$height_list)

sd1 <- strain_data
td1 <- td
