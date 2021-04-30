libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr","reshape2","fossil","tidyr","purrr")
library(data.table)
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/ECC/classes_ecc.R")

source("scripts/ECC/epi-helper-modular.R")
source("scripts/ECC/collecting_eccs.R")
source("scripts/ECC/epi-helper.R")
source("scripts/ECC/ECC-sep_singletons.R")

option_list <- list(
  make_option(c("-a", "--source"), metavar = "file", default = "inputs/processed/source_data.tsv", help = "Source data"),
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/Strain_data.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"), 
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  make_option(c("-x", "--heights"), metavar = "character", default = "0", 
              help = "Comma-delimited string of heights to collect ECCs for"), 
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"),
  make_option(c("-t", "--trio"), metavar = "character", default = "010", 
              help = "source, temporal, geographic coefficients"))

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$trio %>% strsplit(., "-") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)

x <- combos[1]
c1 <- strsplit(x, split = "") %>% unlist() %>% as.numeric() %>% as.list() %>% set_names(c("sigma", "tau", "gamma"))

sigma <- c1$sigma
tau <- c1$tau
gamma <- c1$gamma
cpus <- params$cpus
source_file <- params$source

td <- td1 <- tp1$height_list %>% append(tp2$height_list)

strain_data <- read_tsv(params$strains) %>%
  mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
         Location = paste(Country, Province, City, sep = "_"))

incl_clusters <- td[[2]] %>% group_by(T0) %>% count() %>% ungroup() %>% 
  set_colnames(c("cluster", "cluster_size")) %>% 
  filter(cluster_size < 100) %>% pull(cluster)

td[[2]] <- td[[2]] %>% filter(T0 %in% incl_clusters)
typing_data <- td

strain_data %<>% filter(Strain %in% rownames(td[[2]]))


