#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_epiquant.txt", open="wt")
sink(msg, type="message")

libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr","reshape2","fossil","tidyr",
          "purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/ECC/classes_ecc.R")
source("scripts/ECC/ecc_functions.R")
source("scripts/ECC/collecting_eccs.R")

# Title: "EpiQuant - Salmonella Enteritidis Project (2019-2020)"
# Authors of original work and initial modifications: Ben Hetman, Elissa Giang, Dillon Barker
# Responsible for changes made during and after merge with CGM process: Vasena Jayamanna

checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}


option_list <- list(
  make_option(c("-a", "--source"), metavar = "file", default = NULL, help = "Source data"),
  make_option(c("-b", "--strains"), metavar = "file", default = NULL, help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = NULL, help = "TP1 cluster assignments"),
  make_option(c("-d", "--tp2"), metavar = "file", default = NULL, help = "TP2 cluster assignments"),
  make_option(c("-x", "--heights"), metavar = "character", default = NULL,
              help = "Comma-delimited string of heights to collect ECCs for"),
  make_option(c("-p", "--cpus"), metavar = "numeric", default = NULL, help = "CPUs"),
  make_option(c("-t", "--trio"), metavar = "character", default = NULL,
              help = "source, temporal, geographic coefficients"))

cat(paste0("\n||", paste0(rep("-", 34), collapse = ""), " ECC metric generation ", 
           paste0(rep("-", 34), collapse = ""), "||\nStarted process at: ", Sys.time()))
stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$trio %>% strsplit(., "-") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)

typing_data <- tp1$height_list %>% append(tp2$height_list)

strain_data <- read_tsv(params$strains) %>% 
  mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")), 
         Location = paste(Country, Province, City, sep = "_")) %>% 
  filter(Strain %in% rownames(typing_data[[2]]))

source_pw <- read_tsv(params$source) %>% mutate(Source.Dist = 1- value) %>% select(-value)
  
collected_eccs <- lapply(combos, function(x) {
  c1 <- strsplit(x, split = "") %>% unlist() %>% 
    as.numeric() %>% as.list() %>% set_names(c("sigma", "tau", "gamma"))
  oneCombo(params$strains, params$source, c1$sigma, c1$tau, c1$gamma, params$cpus, typing_data)
}) %>% 
  Reduce(function(...) merge(...), .) %>% as_tibble()

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
timeTaken(pt = "ECC data collection", stopwatch) %>% outputDetails(., newcat = TRUE)

full_set <- collected_eccs %>% left_join(., tp1$proc) %>% left_join(., tp2$proc) %>% 
  mutate(TP1 = ifelse(is.na(TP1), 0, TP1)) %>% arrange(Strain)

write.table(full_set, file = "results/ECCs.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))

