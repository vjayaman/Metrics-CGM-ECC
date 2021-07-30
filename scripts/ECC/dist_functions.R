# df has colnames Strain, T0, dr (in example):
formatForSectioning <- function(df, maxsize) {
  cx <- setdiff(colnames(df), c("Strain", "dr"))
  drsizes <- df %>% select(-Strain) %>% unique() %>% group_by(!!as.symbol(cx)) %>% summarise(n = n())

  # we take the maximum between the largest number of drs a cluster has and the user-provided maxsize
  x1 <- max(drsizes$n)
  if (which.max(c(x1, maxsize)) == 1) {
    return(list(drsizes, "dr" = x1))
  }else {
    return(list(drsizes, "ms" = maxsize))
  }
}


# x should be a list returned by formatForSectioning()
sectionTypingData <- function(x) {
  drsizes <- x[[1]]
  m <- x[[2]] # maxsize
  rows_x <- c()
  total <- 0
  results <- vector("list", nrow(drsizes))

  # for each cluster, we add it to the dataframe rows_x, until we've reached the maxsize (# of drs)
  # if we reach the max size exactly, we add rows_x to the results[[]] list, and refresh the total
  # if the new updated total has more drs than the maxsize, then we set a new element in the results[[]] list
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

  # if we go through all clusters, but end with a dataframe with < maxsize drs,
  # we add this to the results[[]] list as a new element
  if (total != 0) {
    i <- i + 1
    results[[i]] <- rows_x
  }

  results[sapply(results, is.null)] <- NULL
  return(results)
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
    # longitude, then latitude columns, in that order
    input_data %>% select(all_of(cnames)) %>% 
      as.data.frame() %>% earth.dist(dist = TRUE) %>% 
      as.matrix() %>% set_rownames(dm_names) %>% set_colnames(dm_names) %>% return()
  }
}

# assignments <- m$assignments; fpaths <- dists
collectDistances <- function(assignments, parts, fpaths = NULL) {
  df <- parts$drs
  cx <- setdiff(colnames(df), c("Strain", "dr"))
  results <- parts$results
  p <- length(results)
  
  for (j in 1:p) {
    outputMessages(paste0("Working through group ", j, " / ", p))
    cluster_x <- df[df[[cx]] %in% pull(results[[j]], cx),-"Strain"]
    cluster_asmts <- assignments[dr %in% pull(cluster_x, dr)]
    
    outputMessages("   Generating all possible date pair distances ...")
    dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
    
    outputMessages("   Generating all possible lat-long pair distances ...")
    dm_geo <- cluster_asmts %>% select(dr, Longitude, Latitude) %>% 
      distMatrix(., "geo", c("Longitude", "Latitude"))
    
    if (!is.null(fpaths)) {
      if (length(fpaths) == 1) {
        fname <- paste0(fpaths[[1]], "group", formatC(j, width=nchar(p), format="d", flag="0"), ".Rds")
        list(temp = dm_temp, geo = dm_geo) %>% saveRDS(., fname)
      }
    }
    
    rm(dm_temp)
    rm(dm_geo)
    gc()
  }
}

# Indicates length of a process in hours, minutes, and seconds, when given a name of the process 
# ("pt") and a two-element named vector with Sys.time() values named "start_time" and "end_time"
timeTaken <- function(pt, sw) {
  if (is.null(sw[["start_time"]]) & is.null(sw[["end_time"]])) {
    paste0("Neither start nor end time were collected") %>% return()
  }else if (is.null(sw[["end_time"]]) & !is.null(sw[["start_time"]])) {
    paste0("End time was not collected.") %>% return()
  }else if (!is.null(sw[["end_time"]]) & is.null(sw[["start_time"]])) {
    paste0("Start time was not collected.") %>% return()
  }else {
    z <- difftime(sw[['end_time']], sw[['start_time']], units = "secs") %>% as.double()
    m <- 60
    h <- m^2
    
    if (z >= h) {
      hrs <- trunc(z/h)
      mins <- trunc(z/m - hrs*m)
      paste0("\nThe ", pt, " process took ", hrs, " hour(s), ", mins, " minute(s), and ", 
             round(z - hrs*h - mins*m), " second(s).") %>% return()
    }else if (z < h & z >= m) {
      mins <- trunc(z/m)
      paste0("\nThe ", pt, " process took ", mins, " minute(s) and ", round(z - mins*m), " second(s).") %>% return()
    }else {
      paste0("\nThe ", pt, " process took ", round(z), " second(s).") %>% return()
    }  
  }
}


# Outputs the same message in two ways, one is directed to standard output and one to a log file
outputDetails <- function(msg, newcat = FALSE) {
  cat(msg)
  if (newcat) {cat("\n")}
  message(msg)
}


outputMessages <- function(msgs = NULL) {
  if (!is.null(msgs)) {
    cat(paste0("\n",msgs))
  }
}

# # dr_td1 has two columns, one with the clusters, and one with the drs found in each cluster
# dr_td1 <- typing_data[[tpx]] %>% rownames_to_column("Strain") %>% as_tibble() %>%
#   left_join(., dr_matches, by = "Strain") %>% mutate(across(dr, as.character)) %>% select(-Strain)
countDataReps <- function(dr_td1) {
  # Counting data representatives (so we know how much to 
  # multiply each ECC value by to represent all strains)
  cx <- dr_td1 %>% select(-dr) %>% colnames()
  assert("Correctly formatted input dr_td1", length(cx) == 1)
  tallied_reps <- dr_td1 %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
  g_cuts <- left_join(dr_td1, tallied_reps, 
                      by = intersect(colnames(tallied_reps), colnames(dr_td1))) %>%
    unique() %>% mutate(across(dr, as.character))
  return(g_cuts)
}

# note that we sum the values for both directions e.g. (192, 346) and (346, 192)
# and then divide by the total number of pairs (counting both directions)
# this gives the same result as if we had used only one direction
#   - if we did this, we would divide by 2 in the numerator and the denominator --> would cancel
avgDists <- function(g_cuts, dm, cname, newname) {
  x <- gsub("\\.", "_", cname) %>% tolower() %>% paste0(newname, "_avg_", .)
  all_clusters <- g_cuts %>% set_colnames(c("Threshold", "dr", "n"))
  clusters <- unique(all_clusters$Threshold)
  
  lapply(clusters, function(cl) {
    f2 <- all_clusters %>% filter(Threshold == cl)
    dists <- dm[f2$dr, f2$dr]
    if (is.null(dim(dists))) {
      return(0)
    }else {
      dists <- dists * f2$n[col(dists)]
      dists <- dists * f2$n[row(dists)]
      return(sum(dists) / (sum(f2$n)^2))
    }
  }) %>% unlist() %>% tibble(clusters, .) %>% set_colnames(c(newname, x))
}
# -----------------------------------------------------------------------------------------------------
sectionClusters <- function(k, typing_data, m) {
  df <- typing_data[[k]] %>% rownames_to_column("Strain") %>%
    as.data.table() %>% left_join(., m$dr_matches, by = "Strain")
  
  gc()

  results <- formatForSectioning(df, 1000) %>% sectionTypingData(.)
  assert("No clusters overlooked", length(setdiff(pull(df,2), pull(rbindlist(results),1))) == 0)
  
  return(list("drs" = df, "results" = results))
}

avgsFromDM <- function(tpx_dists, groups, g_cuts, type, cname, tpx) {
  lapply(1:length(groups$results), function(j) {
    # clusters in group j
    group_j <- pull(groups$results[[j]], 1)
    # distances between drs in group j
    dm_type <- readRDS(tpx_dists[j])[[type]]
    # drs in the clusters in this group
    cluster_i <- groups$drs[Th %in% group_j]
    # unique drs in this cluster
    unidrs <- unique(cluster_i$dr)
    dists_i <- dm_type[unidrs, unidrs]
    
    # outputMessages(paste0("      Calculating average (not transformed) distances for timepoint ", i))
    g_cuts %>% set_colnames(c("Th", "dr", "n")) %>% 
      filter(Th %in% group_j) %>% 
      set_colnames(colnames(g_cuts)) %>% 
      avgDists(., dists_i, cname, paste0("TP", tpx, "_", colnames(g_cuts)[1])) %>% return()
  }) %>% bind_rows() %>% return()
}