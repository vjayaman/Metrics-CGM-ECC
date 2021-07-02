# df has colnames Strain, T0, dr (in example):
sectionTypingData <- function(df, maxsize) {
  cx <- setdiff(colnames(df), c("Strain", "dr"))
  drsizes <- df %>% select(-Strain) %>% unique() %>% group_by(!!as.symbol(cx)) %>% summarise(n = n())
  
  m <- max(drsizes$n, maxsize)
  rows_x <- c()
  total <- 0
  results <- vector("list", nrow(drsizes))
  
  for (i in 1:nrow(drsizes)) {
    row_i <- drsizes[i,]
    updated <- row_i$n + total
    
    if (updated < m) {
      rows_x <- bind_rows(rows_x, row_i)
      total <- total + row_i$n
      
    }else if (updated == m) {
      rows_x <- bind_rows(rows_x, row_i)
      total <- total + row_i$n
      
      results[[i]] <- rows_x
      rows_x <- c()
      total <- 0
      
    }else { # updated > m
      results[[i]] <- rows_x
      rows_x <- row_i
      total <- row_i$n
    }
  }
  
  if (total != 0) {results[[i]] <- rows_x}
  
  results[sapply(results, is.null)] <- NULL
  
  return(results)
}

sectionClusters <- function(k, typing_data, m) {
  df <- typing_data[[k]] %>% rownames_to_column("Strain") %>%
    as.data.table() %>% left_join(., m$dr_matches, by = "Strain")
  
  gc()
  
  results <- sectionTypingData(df, 1000)
  assert("No clusters overlooked", length(setdiff(pull(df,2), pull(rbindlist(results),1))) == 0)
  
  return(list("drs" = df, "results" = results))
}

distMatrix <- function(input_data, dtype, cnames) {
  dm_names <- pull(input_data, 1) # dr, or Strain
  if (dtype == "temp") {
    dm <- input_data %>% select(all_of(cnames)) %>% pull() %>% 
      dist(diag = FALSE, upper = FALSE, method = "euclidean")
    dm %>% as.matrix(nrow = nrow(input_data), ncol = nrow(input_data)) %>% 
      set_rownames(dm_names) %>% set_colnames(dm_names) %>% return()
    
  }else if (dtype == "geo") {
    # consider using geosphere::distm() for this
    input_data %>% select(all_of(cnames)) %>% 
      as.data.frame() %>% earth.dist(dist = TRUE) %>% 
      as.matrix() %>% set_rownames(dm_names) %>% set_colnames(dm_names) %>% return()
  }
}

collectDistances <- function(k, m, parts, fpath1, fpath2) {
  df <- parts$drs
  cx <- setdiff(colnames(df), c("Strain", "dr"))
  results <- parts$results
  p <- length(results)
  
  min_geo <- min_temp <- Inf
  max_geo <- max_temp <- -Inf
  
  for (j in 1:p) {
    outputMessages(paste0("Working through group ", j, " / ", p))
    cluster_x <- df[df[[cx]] %in% pull(results[[j]], cx),-"Strain"]
    cluster_asmts <- m$assignments[dr %in% pull(cluster_x, dr)]
    
    outputMessages("   Generating all possible date pair distances ...")
    dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
    
    outputMessages("   Generating all possible lat-long pair distances ...")
    dm_geo <- cluster_asmts %>% select(dr, Latitude, Longitude) %>% 
      distMatrix(., "geo", c("Latitude", "Longitude"))
    
    min_temp <- min(min_temp, min(dm_temp))
    max_temp <- max(max_temp, max(dm_temp))
    
    min_geo <- min(min_geo, min(dm_geo + 10))
    max_geo <- max(max_geo, max(dm_geo))
    
    list(temp = dm_temp, geo = dm_geo) %>% 
      saveRDS(., paste0(fpath1, "group", j, ".Rds"))
    
    rm(dm_temp)
    rm(dm_geo)
    gc()
  }
  
  if (k == 2) {
    tibble(temp = list("max" = max_temp, "min" = min_temp), 
           geo = list("max" = max_geo, "min" = min_geo)) %>% 
      saveRDS(., paste0(fpath2, "dist_extremes.Rds"))
  }
}
