#! /usr/bin/env Rscript
msg <- file("logs/logfile_epiquant.txt", open="wt")
sink(msg, type="message")

libs <- c("purrr", "magrittr", "tidyr", "dplyr", "tibble", "methods", "optparse", 
          "readr", "testit", "fossil", "reshape2", "R6")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

source("ECC/functions/classes_ecc.R")
source("ECC/functions/collecting_eccs.R")
source("ECC/functions/epi-helper.R")
source("ECC/functions/ECC-sep_singletons.R") # 010 took 9 minutes and 10 seconds
# source("ECC/ECC-helper.R") # 010 took 1 hour, 1 minute, 12 seconds

# title: "EpiQuant - Salmonella Enteritidis Project (2019-2020)"
# author: "Elissa Giang, National Microbiology Laboratory (Guelph), elissagiang6@gmail.com"

option_list <- list(
  make_option(c("-a", "--source"), metavar = "file", default = NULL, help = "Source data"),
  make_option(c("-b", "--strains"), metavar = "file", default = NULL, help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = NULL, help = "TP1 cluster assignments"), 
  make_option(c("-d", "--tp2"), metavar = "file", default = NULL, help = "TP2 cluster assignments"), 
  make_option(c("-x", "--heights"), metavar = "character", default = "0,5", 
              help = "Comma-delimited string of heights to collect ECCs for"), 
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"),
  make_option(c("-t", "--trio"), metavar = "character", default = "010-001", 
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

td <- tp1$height_list %>% append(tp2$height_list)

lapply(combos, function(x) {
  c1 <- strsplit(x, split = "") %>% unlist() %>% 
    as.numeric() %>% as.list() %>% set_names(c("sigma", "tau", "gamma"))
  oneCombo(params$strains, params$source, c1$sigma, c1$tau, c1$gamma, params$cpus, td)
}) %>% 
  Reduce(function(...) merge(...), .) %>% as_tibble() %>% 
  left_join(., tp1$proc) %>% 
  right_join(., tp2$proc) %>% 
  mutate(TP1 = ifelse(is.na(TP1), 0, TP1)) %>% arrange(Strain) %>% 
  write.table(., file = "results/ECCs.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE)

stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())
cat(timeTaken(pt = "ECC data collection", stopwatch))

cat(paste0("\n||", paste0(rep("-", 30), collapse = ""), " End of ECC metric generation ", 
           paste0(rep("-", 31), collapse = ""), "||\n"))
