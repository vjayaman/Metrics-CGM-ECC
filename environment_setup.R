#! /usr/bin/env Rscript

# This should be run first, to make sure the required packages are installed

dir.create("logs", showWarnings = FALSE)
msg <- file("logs/environment_setup.txt", open="wt")
sink(msg, type="message")

cat(paste0(
  "\n||", paste0(rep("-", 36), collapse = ""), " Environment setup ", 
  paste0(rep("-", 36), collapse = ""), "||"
))

dir.create("results", showWarnings = FALSE)
dir.create("inputs/processed",  recursive = TRUE, showWarnings = FALSE)
dir.create("intermediate_data/TPN/", recursive = TRUE, showWarnings = FALSE)

required_packages <- c("R6", "testit", "optparse", "magrittr", "dplyr", "tibble", "readr", "reshape2", 
                       "fossil", "tidyr", "purrr", "progress", "reader", "data.table", "testthat", 
                       "Rcpp", "RcppArmadillo")

not_installed <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
install.packages(not_installed, quiet = TRUE)

# Testing packages were installed correctly:
x <- lapply(required_packages, require, character.only = TRUE)
names(x) <- required_packages

y <- c(dir.exists("inputs/processed/"), dir.exists("results/"), 
       dir.exists("logs/"), dir.exists("intermediate_data/TPN/"))

# cat(paste0(msg, "\n"))
# message(msg)
if (all(unlist(x)) & all(y)) {
  cat(paste0("\nR packages installed/loaded successfully; required directories created.", 
             "\nSee logfile_env.txt for details\n"))
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

cat(paste0(
  "||", paste0(rep("-", 32), collapse = ""), " End of environment setup ", 
  paste0(rep("-", 33), collapse = ""), "||"
))
