library(data.table); library(magrittr); library(dplyr); library(readr)
a <- read.table("inputs/processed/tp2_clusters.txt", header = TRUE, sep = "\t", #row.names = 1, 
                check.names = FALSE, quote = "", stringsAsFactors = FALSE) %>% as.data.table() %>% 
  select(Strain, "0") %>% set_colnames(c("Strain", "Heightx"))
b <- read_tsv("inputs/processed/strain_info.txt") %>% as.data.table()

d <- merge.data.table(a, b, by = "Strain")

cgms <- readRDS("results/CGM-monthly-intervals.Rds")

# a1 <- unlist(strsplit(cgms$first_tp2_flag[1], split = "_"))[2:3] %>% 
#   nchar() %>% `-`(1) %>% set_names(c("ph", "pc"))
# h_id <- padCol(as.double(hx$h), a1[['ph']], "h")

# read in AVERAGE results and use the TP column to create a TP1_ID column and an identical TP2_ID column: 
# since each TP set is both a TP1 and a TP2, depending on which two timepoints we are comparing
avgs <- readRDS("results/AVGS-monthly-intervals.Rds")
# colnames(avgs)[grep("Avg", colnames(avgs))] %<>% paste0(hx$th, "_", .)
# avgs <- avgs %>% 
#   mutate(tp1_id = paste0("TP1_", h_id, "_", padCol(!!as.symbol(hx$th), a1[['pc']], "c"))) %>% 
#   mutate(first_tp2_flag = tp1_id %>% gsub("TP1", "TP2", .))

# read in ECC results and repeat the process
eccs <- readRDS("results/ECC-monthly-intervals.Rds")# %>% 
  # mutate(tp1_id = paste0("TP1_", h_id, "_", padCol(!!as.symbol(hx$th), a1[['pc']], "c"))) %>% 
  # mutate(first_tp2_flag = tp1_id %>% gsub("TP1", "TP2", .))


cgms[tp1_id == "TP1_h000_c0968"]
cgms[tp1_id == "TP1_h000_c0968"]
