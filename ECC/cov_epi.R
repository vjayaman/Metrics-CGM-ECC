#! /usr/bin/env Rscript
x <- c("optparse", "tibble", "magrittr", "readr", "dplyr", "tidyr", "purrr", "beepr", "stringr", 
       "reshape2", "gplots", "fossil")
lapply(x, require, character.only = TRUE)

source("epi-helper.R")
source("ECC-sep_singletons.R")  # source("helpers/ECC-helper.R")

read_typing <- function(path) {
  read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, quote = "",
             stringsAsFactors = FALSE)
}

# title: "EpiQuant - Salmonella Enteritidis Project (2019-2020)"
# author: "Elissa Giang, National Microbiology Laboratory (Guelph), elissagiang6@gmail.com"

option_list <- list(
  make_option(c("-a", "--source"), metavar = "file", default = "source.tsv", help = "Source data"),
  make_option(c("-b", "--strains"), metavar = "file", default = "strains.tsv", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "../input_data/tp1_clusters.txt", 
              help = "TP1 cluster assignments"), 
  make_option(c("-d", "--tp2"), metavar = "file", default = "../input_data/tp2_clusters.txt", 
              help = "TP2 cluster assignments"), 
  # make_option(c("-c", "--clusters"), metavar = "file", default = "../input_data/", help = "Cluster input file"),
  make_option(c("-h", "--heights"), metavar = "character", default = "0,5", 
              help = "Comma-delimited string of heights to collect ECCs for"), 
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"),
  make_option(c("-t", "--trio"), metavar = "character", default = "010-001", 
              help = "source, temporal, geographic coefficients"))

params <- parse_args(OptionParser(option_list=option_list))

dir.create("results", showWarnings = FALSE)

### Section 3: Incorporating the allele data with the epidemiological data 
oneCombo <- function(strains, source_file, sigma, tau, gamma, clusters, cpus, typing_data) {
  ### Section 2: Generating the EpiMatrix
  strain_data <- 
    read_tsv(strains) %>% 
    mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
           Location = paste(Country, Province, City, sep = "_")
    )
  
  source_pw <- read_tsv(source_file) %>% 
    mutate(Source.Dist = 1- value) %>% 
    select(-value)
  
  temp_pw <- temp_calc(strain_data)
  geog_pw <- geog_calc(strain_data)
  geog_temp <- left_join(geog_pw, temp_pw)

  ## This generates a table of your comparisons based on your epidemiological data (source, time, geographical) 
  ## with the assigned weights of s, t and g and then computes the similarity/distance and generates a matrix
  epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, geog_temp)
  epi.matrix <- EpiMatrix(epi.table)

  # ### Section 3: Incorporating the allele data with the epidemiological data 
  # typing_data_files <- list.files(clusters, full.names = TRUE)
  # names(typing_data_files) <- typing_data_files %>% basename() %>% tools::file_path_sans_ext()
  # 
  # read_typing <- function(path) {
  #   read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, quote = "",
  #              stringsAsFactors = FALSE)
  # }
  # 
  # typing_data <- lapply(typing_data_files, read_typing)
  
  # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  Sys.time()
  eccs <- lapply(typing_data, function(typing_datum) {
    
    g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
    epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
  })
  Sys.time()
  
  newnames <- sapply(strsplit(names(eccs), "_"), `[`, 1)
  
  lapply(1:length(eccs), function(i) {
    eccs[[i]] %>% select(-W_ECC) %>% 
      set_colnames(c(newnames[i], paste0(newnames[i], "_cluster"), paste0(newnames[i], "_cl_size"), 
                     paste0(newnames[i], "_ECC_", sigma, ".", tau, ".", gamma))) %>% ungroup()  
  }) %>% set_names(newnames) %>% return()  
}


hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))
combos <- params$trio %>% strsplit(., "-") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)


a1 <- read_typing(params$tp1)
x1 <- lapply(1:nrow(hx), function(i) a1[,hx$h[i],drop=FALSE] %>% set_colnames(hx$th[i])) %>% 
  set_names(paste0("TP1_", hx$th))
b1 <- read_typing(params$tp2)
td <- lapply(1:nrow(hx), function(i) b1[,hx$h[i],drop=FALSE] %>% set_colnames(hx$th[i])) %>% 
  set_names(paste0("TP2_", hx$th)) %>% append(x1, .)


for (x in combos) {
  c1 <- strsplit(x, split = "") %>% unlist() %>% as.numeric() %>% as.list() %>% set_names(c("sigma", "tau", "gamma"))
  z[[x]] <- oneCombo(params$strains, params$source, c1$sigma, c1$tau, c1$gamma, params$clusters, params$cpus, td)
}

# bind TP1 columns
Reduce(function(...) merge(..., by = c("TP1", "TP1_cluster", "TP1_cl_size")), sapply(z, `[`, 'TP1')) %>% as_tibble() %>% 
  write.table(., file = "results/TP1_ECCs.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE)
# bind TP2 columns
Reduce(function(...) merge(..., by = c("TP2", "TP2_cluster", "TP2_cl_size")), sapply(z, `[`, 'TP2')) %>% as_tibble() %>% 
  write.table(., file = "results/TP2_ECCs.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE)


