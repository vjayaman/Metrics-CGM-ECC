# Testing statements: 
# assert("New method results == original method results (geo)", 
#        identical(geog_calc(strain_data), geog_pw$transformed))
# assert("New method results == original method results (temp)", 
#        identical(temp_calc(strain_data), temp_pw$transformed))

# input_data <- strain_data; dtype <- "temp"; cnames <- "Date"; newcol <- "Temp.Dist"
generateDistances <- function(input_data, dtype, cnames, newcol) {
  dm <- distMatrix(input_data, dtype, cnames)
  raw <- formatMatrix(input_data, dm, newcol)
  transformed <- transformData(dm, dtype) %>% formatMatrix(input_data, ., newcol)
  return(list("raw" = raw, "transformed" = transformed))
}

distMatrix <- function(input_data, dtype, cnames) {
  if (dtype == "temp") {
    dm <- input_data %>% select(all_of(cnames)) %>% pull() %>% 
      dist(diag = FALSE, upper = FALSE, method = "euclidean")
      # bigmemory::as.big.matrix() %>% return()
      # dist(diag = TRUE, upper = TRUE, method = "euclidean")
    dm %>% as.matrix(nrow = nrow(input_data), ncol = nrow(input_data)) %>% return()
    
  }else if (dtype == "geo") {
    # consider using geosphere::distm() for this
    dm <- input_data %>% select(all_of(cnames)) %>% as.data.frame() %>% earth.dist(dist = TRUE)
    dm %>% as.matrix() %>% return()
  }
}

transformData <- function(dm, dtype) {
  logdata <- dm %>% add(10) %>% log10()
  
  if (dtype == "temp") {
    
    logdata[logdata == -Inf] <- 0
    if (max(logdata) == 0) {
      logdata <- 0
    }else {
      logdata <- ((logdata - min(logdata)) / (max(logdata) - min(logdata)))
    }
    
  }else if (dtype == "geo") {
    
    if(max(logdata) == 1){
      logdata[1:nrow(logdata), 1:nrow(logdata)] <- 0
    } else {
      logdata <- ((logdata-min(logdata)) / (max(logdata)-min(logdata)))
    }
  }
  
  return(logdata)
}

formatMatrix <- function(input_data, tdata, cname) {
  #### Import row and column names from the original datafile and melt data for easy reading ####
  matrix_names <- input_data %>% select(Strain) %>% pull()
  rownames(tdata) <- matrix_names
  colnames(tdata) <- matrix_names
  
  tdata %>% 
    melt(varnames = c("Strain.1", "Strain.2"), value.name = cname) %>% 
    as_tibble() %>% return()
}
