libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", 
           "scripts/ECC/dist_functions.R")
invisible(sapply(files, source))

option_list <- list(
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-k", "--set"), metavar = "numeric", default = 2, help = "Timepoint distances to verify"),
  make_option(c("-t", "--type"), metavar = "character", default = "temp", help = "temp or geo, which distances to verify") 
)

params <- parse_args(OptionParser(option_list=option_list))

k <- params$set
outputMessages(paste0("Validating data for TP", k, " set, ", params$type, " distances ...\n"))

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

m <- invisible(read_tsv(params$strains)) %>% processedStrains()

if (params$type == "temp") {
  type <- c("temp", "Temp.Dist", "Date")
}else {
  type <- c("geo", "Geog.Dist", "Longitude,Latitude")
}

parts <- sectionClusters(k, typing_data, m)

outputMessages(paste0("Reading in intra-cluster distances from intermediate_data/TP", k, "/dists/"))
intra_cl <- paste0("intermediate_data/TP", k, "/dists/") %>% 
  list.files(., full.names = TRUE, pattern = "group") %>% 
  lapply(., function(fname) {readRDS(fname)})

intra_temps <- lapply(1:length(intra_cl), function(i) {
  intra_cl[[i]][[type[1]]] %>% as.data.frame() %>% 
    rownames_to_column("dr") %>% as.data.table() %>% 
    melt.data.table(id.vars = "dr") %>% 
    set_colnames(c("dr1", "dr2", type[2]))
}) %>% bind_rows() %>% unique()

# tmp <- intra_temps %>% mutate(both = paste0(dr1, "-", dr2))
# x1 <- table(tmp$both)
# x2 <- x1[x1 > 1] %>% names
# assert("No repeated dr pairwise distances, when intra-cluster distances are merged", length(x2) == 0)

outputMessages(paste0("Reading in inter-cluster distances from intermediate_data/TP", k, "/dists/"))
inter_temps <- paste0("intermediate_data/TP", k, "/dists/inter_dists_", type[1], ".Rds") %>% 
  readRDS() %>% as.data.frame() %>% rownames_to_column("dr") %>% as.data.table() %>% 
  melt.data.table(id.vars = "dr") %>% 
  set_colnames(c("dr1", "dr2", type[2]))

compiled_temps <- bind_rows(inter_temps, intra_temps) %>% unique() %>% 
  mutate(across(dr2, as.character)) %>% as.data.table()

dr_matches <- parts$drs[,c("Strain","dr")] %>% as_tibble()

outputMessages("Merging results - those collected the non-redundant way")
p1 <- left_join(compiled_temps, m$dr_matches, by = c("dr1" = "dr")) %>% rename(Strain1 = Strain)
p2 <- right_join(dr_matches, p1, by = c("dr" = "dr2")) %>% rename(Strain2 = Strain)
p3 <- p2 %>% as.data.table() %>% select(Strain1, Strain2, all_of(type[2])) %>% 
  set_colnames(c("Strain1", "Strain2", "Dist"))

outputMessages("Manually preparing the (strain) pairwise distances - note, require smallish datasets")
dm_all_rows <- m$strain_data %>% select(Strain, unlist(strsplit(type[3],","))) %>% 
  distMatrix(., type[1], unlist(strsplit(type[3],",")))

all_rows <- dm_all_rows %>% as.data.frame() %>% rownames_to_column("Strain") %>% 
  as.data.table() %>% melt.data.table(id.vars = "Strain") %>% 
  set_colnames(c("Strain1", "Strain2", "Dist")) %>% 
  mutate(across(Strain2, as.character)) %>% 
  as.data.table()

outputMessages("Comparing the results ...")
a1 <- all_rows %>% rename(Actual_Dist = Dist)
a2 <- p3 %>% rename(New_Dist = Dist)
a3 <- inner_join(a1, a2, by = c("Strain1", "Strain2"))
difs <- a3$Actual_Dist - a3$New_Dist

if (all(abs(difs) < 1e-10)) {
  outputMessages("\nSuccess!\n   - Distances generated the non-redundant way match those generated manually :)\n")
}else {
  outputMessages(paste0("\nProblem!\n   - Mismatch in comparing distances generated the ", 
                        "non-redundant way to the manual method!\n"))
}
