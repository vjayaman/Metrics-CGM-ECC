
libs <- c("R6","optparse","magrittr","tibble", "dplyr", "readr", "testit", "data.table", 
          "fossil", "tidyr", "purrr")
y <- suppressMessages(lapply(libs, require, character.only = TRUE))

option_list <- list(
  make_option(c("-e", "--ECCs"), metavar = "file", default = "results/ECCs.tsv", help = "ECC result file"),
  make_option(c("-a", "--tp1"), metavar = "file",
              default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-b", "--tp2"), metavar = "file",
              default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-p", "--timepoint"), metavar = "numeric", default = 1, 
              help = "Test ECCs of this timepoint manually"), 
  make_option(c("-s", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain metadata file"), 
  make_option(c("-t", "--smalldataset"), metavar = "logical", default = FALSE))

arg <- parse_args(OptionParser(option_list=option_list))

source("scripts/ECC/ecc_functions.R")
source("scripts/ECC/classes_ecc.R")
source("tests/testthat/check_ecc_functions.R")
source("scripts/ECC/dist_functions.R")
# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

gc()

outputMessages("Manually checking the ECC results file ...")
outputMessages("The default is to randomly select 10 clusters at TP1 and 10 at TP2 to verify ECCs")
outputMessages("To check all clusters (only recommended if your dataset is relatively small, < 5000 strains),")
outputMessages("  use argument --smalldataset TRUE\n")

tp1 <- Timepoint$new(arg$tp1, "tp1")$filedata[,"0",drop=FALSE]
tp2 <- Timepoint$new(arg$tp2, "tp2")$filedata[,"0",drop=FALSE]
typing_data <- list("tp1" = tp1, "tp2" = tp2)

strain_data <- suppressMessages(read_tsv(arg$strains))
# strain_data <- readData(arg$strains, check_enc = TRUE)

returned_eccs <- readData(arg$ECCs) %>% as.data.table()

triples <- colnames(returned_eccs) %>% grep("ECC", ., value = TRUE) %>% 
  strsplit(., split = "ECC.") %>% sapply(., '[[', 2) %>% unique()

combos <- expand_grid(timepoint = c("TP1", "TP2"), triples) %>% 
  separate(., triples, sep = "[.]", into = c("phi", "tau", "gamma")) %>% select(-phi)

if (arg$smalldataset) {
  outputMessages("Collecting the pairwise distances (all strain pairs - i.e. manually)")
  transformed_dists <- strain_data %>% testDists()
  
  results <- lapply(1:nrow(combos), function(i) {
    tpcname <- combos$timepoint[i]
    tpx <- tpcname %>% gsub("TP", "", .) %>% as.integer()
    g_cuts <- typing_data[[tpx]] %>% set_colnames(paste0("T", colnames(.))) %>%
      rownames_to_column("genome")
    
    tau <- combos$tau[i] %>% as.integer()
    gamma <- combos$gamma[i] %>% as.integer()
    
    outputMessages(paste0("Manually collecting ECCs for tau = ", tau, ", gamma = ",
                          gamma, ", at timepoint ", tpcname, " ..."))
    original_eccs <- transformed_dists %>% testEpiMelt(., tau, gamma) %>% testECCs(g_cuts, .)
    hx <- unique(original_eccs$cut)[1]
    dfx <- original_eccs %>% select(-cut) %>%
      set_colnames(c(paste(tpcname, hx, sep = "_"),
                     paste(tpcname, hx, "Size", sep = "_"),
                     paste(tpcname, hx, paste0("ECC.0.", tau, ".", gamma), sep = "_")))
    dfx[is.nan(pull(dfx, 3)),3] <- NA
    
    typing_data[[tpx]] %>% set_colnames(colnames(dfx)[1]) %>%
      rownames_to_column("Strain") %>%
      left_join(., dfx, by = intersect(colnames(dfx), colnames(.))) %>% return()
  })
  # saveRDS(results, "results/manually_gen_eccs.Rds")
  
}else {
  outputMessages("Manually prepping minimal distances necessary and sampling clusters ...\n")
  geo_dists <- strain_data %>% select(Longitude, Latitude) %>% unique() %>% 
    rownames_to_column("id") %>% distMatrix(., "geo", c("Longitude", "Latitude"))
  temp_dists <- strain_data %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
    select(Date) %>% unique() %>% rownames_to_column("id") %>% distMatrix(., "temp", "Date")
  
  extremes <- list(maxt = max(temp_dists), mint = min(temp_dists), 
                   maxg = max(geo_dists), ming = min(geo_dists))
  temps <- list(maxval = max(temp_dists), minval = min(temp_dists))
  geos <- list(maxval = max(geo_dists), minval = min(geo_dists))
  
  # Sampling 10 clusters for testing
  check_clusters <- list(tp1 %>% pull(1) %>% unique() %>% na.omit() %>% sample(., 10), 
                         tp2 %>% pull(1) %>% unique() %>% na.omit() %>% sample(., 10))
  
  tmp <- colnames(returned_eccs)
  colnames(returned_eccs)[c(2,5)] <- c("tp1_cl", "tp2_cl")
  x1 <- returned_eccs[tp1_cl %in% check_clusters[[1]] | tp2_cl %in% check_clusters[[2]]] %>% pull(Strain)
  
  while (length(x1) > 800) {
    check_clusters <- list(tp1 %>% pull(1) %>% unique() %>% na.omit() %>% sample(., 10), 
                           tp2 %>% pull(1) %>% unique() %>% na.omit() %>% sample(., 10))
    x1 <- returned_eccs[tp1_cl %in% check_clusters[[1]] | tp2_cl %in% check_clusters[[2]]] %>% pull(Strain)
  }
  
  outputMessages(paste0("Manually verifying ECCs for the following clusters: \n   ", 
                        "TP1: ", paste0(sort(check_clusters[[1]]), collapse = ", "), "\n   ", 
                        "TP2: ", paste0(sort(check_clusters[[2]]), collapse = ", "), "\n   ", 
                        length(x1), " / ", nrow(returned_eccs), " strains\n"))
  colnames(returned_eccs) <- tmp
  
  transformed_dists <- strain_data %>% filter(Strain %in% x1) %>% testDists(., temps, geos)
  
  results <- lapply(1:nrow(combos), function(i) {
    tpcname <- combos$timepoint[i]
    tpx <- tpcname %>% gsub("TP", "", .) %>% as.integer()
    g_cuts <- typing_data[[tpx]] %>% set_colnames(paste0("T", colnames(.))) %>%
      rownames_to_column("genome")
    cnames <- colnames(g_cuts)
    
    tau <- combos$tau[i] %>% as.integer()
    gamma <- combos$gamma[i] %>% as.integer()
    
    outputMessages(paste0("Manually collecting ECCs for tau = ", tau, ", gamma = ",
                          gamma, ", at timepoint ", tpcname, " ..."))
    partA <- transformed_dists %>% testEpiMelt(., tau, gamma)
    subset_cuts <- g_cuts %>% set_colnames(c("a","b")) %>% 
      filter(b %in% check_clusters[[1]]) %>% set_colnames(cnames)
    original_eccs <- subset_cuts %>% testECCs(., partA)
    
    hx <- unique(original_eccs$cut)[1]
    dfx <- original_eccs %>% select(-cut) %>%
      set_colnames(c(paste(tpcname, hx, sep = "_"),
                     paste(tpcname, hx, "Size", sep = "_"),
                     paste(tpcname, hx, paste0("ECC.0.", tau, ".", gamma), sep = "_")))
    dfx[is.nan(pull(dfx, 3)),3] <- NA
    
    subset_cuts %>% set_colnames(c("Strain", colnames(dfx)[1])) %>% 
      left_join(., dfx, by = intersect(colnames(dfx), colnames(.))) %>% return()
  })
  # saveRDS(results, "results/manually_gen_eccs.Rds")
}

original_eccs <- suppressMessages(Reduce(full_join, results)) %>% 
  as.data.table() %>% arrange(Strain)
new_eccs <- returned_eccs %>% select(colnames(original_eccs)) %>% 
  filter(Strain %in% original_eccs$Strain) %>% arrange(Strain)

dfz <- original_eccs %>% set_colnames(gsub("TP", "o_TP", colnames(.))) %>% 
  full_join(., new_eccs, by = "Strain")

outputMessages(paste0("\nChecking that the manually collected columns in the ECC ", 
                      "results file are either: \n", 
                      "   - identical to those generated by the program (cluster, size)\n", 
                      "   - or within 1e-10 to the generated values (the ECCs and average distances)\n"))

checked_all_cols <- lapply(2:length(colnames(new_eccs)), function(j) {
  cname <- colnames(new_eccs)[j]
  newname <- paste0("o_", cname)
  difs <- abs(pull(dfz, cname) - pull(dfz, newname)) %>% na.omit()
  if (grepl("ECC", cname)) {
    return(all(difs < 1e-10))
  }else {
    return(all(difs == 0))
  }
}) %>% unlist() %>% set_names(colnames(new_eccs)[-1])

if (all(checked_all_cols)) {
  outputMessages("All ECC columns match results generated manually :)\n\n\n")
}else {
  outputMessages("Not all columns of the generated ECC files match the manually generated result.\n\n\n")
}
