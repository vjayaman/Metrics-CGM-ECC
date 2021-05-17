#! /usr/bin/env Rscript

dir.create("logs", showWarnings = FALSE)
# This should be run first, to make sure the required packages are installed

required_packages <- c("R6", "testit", "optparse", "magrittr", "dplyr", "tibble", "readr", "reshape2", 
                       "fossil", "tidyr", "purrr", "progress", "reader", "data.table")

not_installed <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
install.packages(not_installed, quiet = TRUE)

# Testing packages were installed correctly:
x <- lapply(required_packages, require, character.only = TRUE)
names(x) <- required_packages

dir.create("results", showWarnings = FALSE)
dir.create(file.path("inputs", "processed"))

if (all(unlist(x))) {
  cat("\nEnvironment set up successful.\n")
}else {
  # if ("plotly" %in% names(which(x == FALSE))) {
  #   message("Linux libraries missing. Try: \n")
  #   message("   $ sudo apt-get update")
  #   message("   $ sudo apt-get install libcurl4-openssl-dev")
  #   message("   $ sudo apt-get install libssl-dev\n")
  #   message("Then run env_setup.R again.")
  # }
  cat("\nNot all packages were installed successfully. Please see logfile_env.txt for details.")  
}

