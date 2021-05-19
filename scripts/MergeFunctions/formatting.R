writeData <- function(fp, df) {write.table(df, fp, row.names = FALSE, quote = FALSE, sep = "\t")}

checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

readData <- function(fp) {
  read.table(fp, stringsAsFactors = FALSE, header = TRUE, fileEncoding = checkEncoding(fp)) %>% 
    as.data.table() %>% return()
}

getAverage <- function(df, grp, avg, newname) {
  df %>% group_by({{grp}}) %>% summarise(avg_col = mean({{avg}})) %>% rename_at(2, ~newname)
}

getDistCols <- function(cnames, x, y, returnVal = TRUE) {
  grep(paste0("(?=.*", x, ")(?=.*", y, ")"), cnames, value = returnVal, perl = TRUE) %>% return()
}

replaceDistName <- function(x) {
  a1 <- ifelse(grepl("TP1", x), "TP1", "TP2")
  if (grepl("temp", x)) {
    paste0(a1, " temp average cluster distance (days)") %>% return()
  }else {
    paste0(a1, " geo average cluster distance (km)") %>% return()
  }
}
