#! /usr/bin/env Rscript

# Current working directory should be Metrics-CGM-ECC/

msg <- file("logs/logfile_distmatrices.txt", open="wt")
sink(msg, type="message")

source("scripts/ecc_opener.R")

# Note: dr stands for data representative
# in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
# cat(paste0("\n\nStep 1:"))
# cat(paste0("\n   Note that the source coefficent is always 0 in this version"))
# outputMessages("   Removing redundancy (comparing date, lat-long, etc. - not every pair of strains)")
# outputMessages("   Identifying which strains match which non-redundant 'data representatives'")

# COLLECT dist matrices using TP2 clusters ---------------------------------------------------------

for (k in 1:2) {
  dir.create(paste0("results/TP", k), showWarnings = FALSE)
  dir.create(paste0("results/TP", k, "/dists"), showWarnings = FALSE)
  
  outputMessages(paste0("\nCollecting and saving distances for groups at TP", k))
  dists <- paste0("results/TP", k, "/dists/")
  extremes <- "results/"
  sectionClusters(k, typing_data, m) %>% collectDistances(k, m, ., dists, extremes)
}

outputMessages("\nFinished saving distance matrices.")
