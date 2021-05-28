checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

# Outputs the same message in two ways, one is directed to standard output and one to a log file
outputDetails <- function(msg, newcat = FALSE) {
  cat(msg)
  if (newcat) {cat("\n")}
  message(msg)
}

# Writes data to a given location, saves as tab-delimited text file
writeData <- function(df, filepath) {
  write.table(df, filepath, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

# Creates a lookup table of integer cluster (in case clusters with character names are provided)
# Returns lookup table and new integer columns, with integer thresholds as well, from 0 to # of thresholds - 1
intClusters <- function(df) {
  list_lookup <- lapply(2:ncol(df), function(j) {
    original_values <- df %>% pull(j) %>% unique()
    tibble(original_values, 1:length(original_values)) %>% 
      set_names(c(names(df)[j], paste0("new_", j-2)))
  })
  
  lookup_tbl <- suppressMessages(plyr::join_all(append(list(df), list_lookup), type = "left")) %>% as_tibble()
  cnames <- grep("Strain|new", colnames(lookup_tbl), value = TRUE)
  
  list(lookup_tbl, lookup_tbl[,cnames] %>% set_colnames(gsub("new_", "", colnames(.)))) %>% 
    set_names(c("lookup_table", "new_cols"))%>% return()
}
