checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

# Outputs the same message in two ways, one is directed to standard output and one to a log file
outputDetails <- function(msg, newcat = FALSE) {
  cat(msg)
  if (newcat) {cat("\n")}
  message(msg)
}

# function for reading raw strains and time point clusters
readData <- function(path, check_enc = TRUE) {
  if (check_enc) {
    enc <- checkEncoding(file.path(path))
  }else {
    enc <- ""
  }
  
  file.path(path) %>% 
    read.table(sep="\t", header=TRUE, fileEncoding=enc, fill=TRUE, quote="") %>% 
    as_tibble() %>% return()
}

# Writes data to a given location, saves as tab-delimited text file
writeData <- function(df, filepath) {
  write.table(df, filepath, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

# Creates a lookup table of integer cluster (in case clusters with character names are provided)
# Returns lookup table and new integer columns, with integer thresholds as well, from 0 to # of thresholds - 1
# If the clusters are already integers, simply renames the columns names to be 0 to # of thresholds
# and uses the given integer cluster names
intClusters <- function(df) {
  already_integers <- lapply(2:ncol(df), function(cx) {
    df[,cx] %>% pull() %>% is.integer() %>% all()
  }) %>% unlist() %>% all()
  
  if (!already_integers) {
    list_lookup <- lapply(2:ncol(df), function(j) {
      original_values <- df %>% pull(j) %>% unique()
      tibble(original_values, 1:length(original_values)) %>% 
        set_names(c(names(df)[j], paste0("new_", j-2)))
    })
    
    lookup_tbl <- suppressMessages(plyr::join_all(append(list(df), list_lookup), type = "left")) %>% as_tibble()
    cnames <- grep("Strain|new", colnames(lookup_tbl), value = TRUE)
    
    list(lookup_tbl, lookup_tbl[,cnames] %>% set_colnames(gsub("new_", "", colnames(.)))) %>% 
      set_names(c("lookup_table", "new_cols"))%>% return()  
  }else {
    df2 <- df[,2:ncol(df)] %>% set_colnames((1:ncol(.))-1)
    df3 <- bind_cols(df[,1], df2)
    list(lookup_table = bind_cols(df, df2), new_cols = df3) %>% return()
  }
}

updateStrains <- function(type, strain_data, time1, time2, sizes) {
  # Strains with metadata and defined lineage info at TP1
  x1 <- intersect(strain_data$Strain, time1$Strain)
  time1 <- time1 %>% filter(Strain %in% x1)
  
  # Strains with metadata and defined lineage info at TP1
  x2 <- intersect(strain_data$Strain, time2$Strain)
  time2 <- time2 %>% filter(Strain %in% x2)
  
  # Strains that have defined lineage info
  strain_data <- strain_data %>% filter(Strain %in% c(x1, x2))
  
  sizes %<>% bind_rows(tibble(type=type, a=nrow(strain_data), b=nrow(time1), d=nrow(time2)))
  list(sd = strain_data, t1 = time1, t2 = time2, sizes = sizes) %>% return()
}

strainsInSingletons <- function(df, cname) {
  singletons <- df %>% 
    select(all_of(cname)) %>% 
    group_by(!!as.symbol(cname)) %>% 
    summarise(n = n()) %>% 
    filter(n == 1) %>% 
    pull(all_of(cname))
  
  df %>% filter(!!as.symbol(cname) %in% singletons) %>% pull(Strain) %>% return()
}
