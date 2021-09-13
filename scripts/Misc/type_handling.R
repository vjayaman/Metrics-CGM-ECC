checkTypes <- function(df) {
  # verify types
  x <- df %>% select(tp1_cl_size, tp2_cl_size, type) %>% as_tibble() %>% 
    rownames_to_column("id") %>% mutate(across(id, as.integer))
  
  x1 <- x %>% filter(tp1_cl_size > 2, tp2_cl_size > 2, tp1_cl_size == tp2_cl_size)
  x2 <- x %>% filter(tp1_cl_size > 2, tp2_cl_size > 2, tp2_cl_size > tp1_cl_size)
  x3 <- x %>% filter(tp1_cl_size < 3, tp2_cl_size > 2)
  x4 <- x %>% filter(tp1_cl_size < 3, tp2_cl_size < 3)
  
  assert("Assigned types are correct", all(
    unique(x1$type) == "Type1", 
    unique(x2$type) == "Type2", 
    unique(x3$type) == "Type3", 
    unique(x4$type) == "Type4"))
  
  return(identical(sort(c(x1$id, x2$id, x3$id, x4$id)), x$id))
}

colsToInherit <- function(df) {
  grep("tp1|TP1", colnames(df), value = TRUE) %>% 
    grep("found_in", ., value = TRUE, invert = TRUE) %>% return()
}

inheritTP1Data <- function(typeX, clustersX) {
  cnames <- colsToInherit(typeX)
  
  for (k in 1:length(clustersX)) {
    to_inherit <- typeX[first_tp2_flag == clustersX[k] & !is.na(tp1_id), ..cnames] %>% unique()
    absent_cases <- grep("Absent", to_inherit$tp1_id, value = TRUE)
    to_inherit <- to_inherit[!(tp1_id %in% absent_cases)]
    
    if (nrow(to_inherit) == 1) {
      rowsY <- typeX[   first_tp2_flag == clustersX[k] & is.na(tp1_id), -(..cnames)] %>% cbind(to_inherit)
      typeX <- typeX[ !(first_tp2_flag == clustersX[k] & is.na(tp1_id)) ] %>% bind_rows(., rowsY)
      
    }#else {
      #stop("Number of TP1 clusters to inherit from is 0 or > 1")
    #}  
  }
  return(typeX)
}

type2Inheritance <- function(typeX) {
  clustersX <- typeX$first_tp2_flag %>% unique()
  inheritTP1Data(typeX, clustersX) %>% return()
}

# for Type 3, 
#   - in some cases we have TP1 singletons to inherit cluster info from (but no ECC), TP1 ECC := 1
#   - and in other cases no TP1 cluster to inherit from --> in this case, TP1 ECC := 1
type3Inheritance <- function(typeX) {
  
  x1 <- getDistCols(colnames(typeX), "TP1", "ECC", FALSE)
  assert("All TP1 clusters have ECCs of NA", all(is.na(typeX[,..x1])))
  set(typeX, j = x1, value = 1)

  # for those that have NA tp1_id, we know they do not exist at TP1 for Type 3, so we just set
  # the TP1 ECC to 1, there is no TP1 cluster info to inherit
  clustersX <- typeX[!is.na(tp1_id)] %>% pull(first_tp2_flag) %>% unique()

  # if 0, cluster k is entirely novel, so no TP1 data to inherit, we just force ECCs to be 1
  inheritTP1Data(typeX, clustersX) %>% return()
  
  # assert("All TP1 clusters that are either singletons or non-existent have value 1",
  #        all(typeX[,..x1] == 1))
  # check_inherited <- all(!is.na(select(typeX, grep("ECC", colnames(typeX)))))
  # if (check_inherited) {return(typeX)}else {return(NULL)}
}

# FORMATTING -------------------------------------------------------------------------------------
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
