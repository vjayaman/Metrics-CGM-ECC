#"Metadata file", "TPN data", "Analysis inputs (details)", 
#"File with Strain column and Pango_lineage column"))
arg <- data.frame(
  metadata = "../inputs/processed/strain_info.txt", 
  tpn = "../inputs/processed/allTP2.Rds", 
  details = "../inputs/form_inputs.txt"
)

params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
  set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
              "th","nsTP2", "temp_win","cnames","int_type","divs", "wklbl",
              "coeffs", "numcl", "mincl", "clustby", "trheatmaps", "lowcol", 
              "midcol", "highcol"))

# # For all scripts/ files other than 1_prepare_inputs.rscript: 
# option_list <- list(
#   make_option(c("-m", "--metadata"), metavar = "file", 
#               default = "inputs/processed/strain_info.txt", help = "Metadata file"),
#   make_option(c("-n", "--tpn"), metavar = "file", 
#               default = "inputs/processed/allTP2.Rds", help = "TPN data"), 
#   make_option(c("-d", "--details"), metavar = "file", 
#               default = "inputs/form_inputs.txt", help = "Analysis inputs (details)"))
# 
# arg <- parse_args(OptionParser(option_list=option_list)); rm(option_list)
# 
# params <- readLines(arg$details, warn = FALSE) %>% strsplit(., split = ": ") %>%
#   set_names(c("reg","cou","has_lin", "has_date","has_prov","prov",
#               "th","nsTP2", "temp_win","cnames","int_type","divs", "wklbl",
#               "coeffs", "numcl", "mincl", "clustby", "trheatmaps", "lowcol", 
#               "midcol", "highcol"))