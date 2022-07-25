
libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source))

assert("Distances were collected and saved", file.exists("intermediate_data/TPN/extreme_dists.Rds"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ",
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))

option_list <- list(
  make_option(c("-m", "--metadata"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-t", "--trio"), metavar = "character", default = "0-1-0",
              help = "source, temporal, geographic coefficients"), 
  make_option(c("-i", "--intervaltype"), metavar = "char", default = "monthly", 
              help = "Type of intervals, choices are: weekly, monthly, multiset. If multiset, provide a time to split the dataset at."))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))
m <- read_tsv(params$metadata) %>% processedStrains()

metadata <- m$strain_data %>% 
  mutate(YearMonth = format(Date, "%Y-%m")) %>% 
  mutate(Week = strftime(Date, format = "%V")) %>% 
  select(-TP1, -TP2) %>% arrange(Week) %>% as.data.table()

tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
f2 <- tp2$filedata %>% rownames_to_column("isolate") %>% as.data.table() %>% 
  select(isolate, all_of(params$heights)) %>% arrange(isolate)

extremes <- readRDS("intermediate_data/TPN/extreme_dists.Rds")

source("scripts/interval_prep.R")

typing_data <- lapply(1:length(interval_list), function(i) {
  n1 <- as.character(interval_list[i])
  tpkstrains <- metadata[get(interval) <= n1]$Strain
  dfz <- tp2$filedata %>% rownames_to_column("isolate") %>%
    select(isolate, all_of(params$heights)) %>%
    filter(isolate %in% tpkstrains) %>% column_to_rownames("isolate")
  dfz[,hx$h[1],drop=FALSE] %>% set_colnames(hx$th[1])
}) %>% set_names(as.character(interval_list))

td <- typing_data[[length(typing_data)]] %>% rownames_to_column("Strain") %>% as.data.table()
parts <- m$dr_matches %>% filter(Strain %in% td$Strain) %>% 
  left_join(td, ., by = "Strain") %>% sectionClusters(.)
df <- parts$drs
cx <- setdiff(colnames(df), c("Strain", "dr"))
results <- parts$results

dfx <- params$trio %>% strsplit(., ",") %>% unlist() %>% 
  expand.grid(x = ., k = interval_list) %>% as.data.frame()

check_epis <- vector(mode = "list", length = nrow(dfx))

for (i in 1:nrow(dfx)) {
  cat(paste0("\nStep ", i, " / ", nrow(dfx), ":\n"))
  
  k <- as.character(dfx$k[i])
  tpkstrains <- metadata[get(interval) <= k]$Strain
  key_cls <- parts$drs[Strain %in% tpkstrains] %>% select(-Strain, -dr) %>% pull() %>% unique()
  y <- lapply(parts$results, function(x) any(key_cls %in% pull(x, 1))) %>% unlist()
  
  k_drs <- m$dr_matches %>% filter(Strain %in% tpkstrains) %>% pull(dr)
  
  c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "-") %>% unlist() %>% as.numeric()
  tau <- c1[2]
  gamma <- c1[3]
  
  fnames <- names(y[y])
  
  check_epis[[i]] <- lapply(fnames, function(f) {
    cluster_y <- df[df[[cx]] %in% pull(results[[f]], cx),-"Strain"]
    dms <- readRDS(paste0("intermediate_data/TP", k, "/dists/group", f, ".Rds"))
    
    transformed_dists <- collectTransforms(dms, extremes)
    
    selected_tp <- m$strain_data %>% filter(Strain %in% tpkstrains)
    
    strain_data <- selected_tp
    tpx <- k
    cluster_x <- cluster_y[dr %in% k_drs]
    
    epi_table <- transformed_dists %>% 
      mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
      select(dr1, dr2, Total.Dist) %>% as.data.table()
    invisible(gc())
    
    # outputMessages("      Generating epi distance matrix ...")
    epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
    epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)]) 
    rownames(epi_matrix) <- colnames(epi_matrix)
    rm(epi_table)
    invisible(gc())
    
    epi_melt <- as.matrix(1-epi_matrix) %>% 
      as.data.table(., keep.rownames = TRUE) %>% 
      melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>% 
      set_colnames(c("dr_1", "dr_2", "value")) %>% as.data.table()
    
    rm(epi_matrix)
    invisible(gc())
    
    cx <- colnames(cluster_x)[1]
    tallied_reps <- cluster_x %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
    cnames <- intersect(colnames(tallied_reps), colnames(cluster_x))
    g_cuts <- left_join(cluster_x, tallied_reps, by = cnames) %>%
      unique() %>% mutate(across(dr, as.character))
    
    dr_names <- g_cuts %>% select(dr) %>% pull() %>% unique()
    dr_assignments <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))
    
    epi_melt_joined <-
      expand_grid(dr_names, dr_names, .name_repair = function(x) {c("dr_1", "dr_2")}) %>%
      inner_join(., epi_melt, by = c("dr_1", "dr_2")) %>% as.data.table()
    
    a1 <- epi_melt %>% arrange(dr_1, dr_2)
    rm(epi_melt)
    
    b1 <- epi_melt_joined %>% arrange(dr_1, dr_2)
    rm(epi_melt_joined)
    gc()
    
    identical(a1, b1)
  }) %>% unlist() %>% all()
}

assert("epi_melt and epi_melt_joined are identical", all(unlist(check_epis)))
