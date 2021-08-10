# msg <- file("logs/logfile_datacollection.txt", open="wt")
# sink(msg, type="message")

libs <- c("R6", "tibble", "optparse", "magrittr", "dplyr", "reshape2", "progress", 
          "testit", "data.table", "readr")
y <- suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))

files <- paste0("scripts/CGM") %>% list.files(., full.names = TRUE)
invisible(sapply(files, source))

# READING IN THE INPUTS ----------------------------------------------------------------------------------------
# Change the default values to read in your own files, or feed through terminal arguments
option_list <- list(
  make_option(c("-m", "--strains"), metavar = "file", default = "inputs/processed/strain_info.txt", help = "Strain metadata file"), 
  make_option(c("-a", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "Time point 1 file name (TP1)"),
  make_option(c("-b", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "Time point 2 file name (TP2)"),
  make_option(c("-x", "--heights"), metavar = "character", default = "0",
              help = paste0("A string of comma-delimited numbers, e.g. '50,75,100' to ", 
                            "use as heights for metric generation (strain table outputs)")))
arg <- parse_args(OptionParser(option_list=option_list))

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)
