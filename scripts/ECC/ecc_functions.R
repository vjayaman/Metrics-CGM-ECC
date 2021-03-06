source("scripts/ECC/epicohesionversions.R")

# -----------------------------------------------------------------------------------------------------

### Incorporating the allele data with the epidemiological data 
# strain_data <- selected_tp; tpx <- k_i; cluster_y <- cluster_x[dr %in% k_drs]
epiCollectionByClusterV2 <- function(strain_data, tau, gamma, tr_dists, tpx, cluster_y) {
  # cat(paste0("\n   Collecting ECC values for temporal = ", tau, ", geo = ", gamma))
  
  epi_melt <- epiMat(tr_dists$temp, tr_dists$geo, as.double(tau), as.double(gamma)) %>% 
    set_rownames(rownames(tr_dists$temp)) %>% set_colnames(colnames(tr_dists$temp))

  invisible(gc())
  
  # outputMessages("      Incorporating the metadata with the clusters (typing data) ...")
  # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
  cx <- colnames(cluster_y)[1]
  tallied_reps <- cluster_y %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
  cnames <- intersect(colnames(tallied_reps), colnames(cluster_y))
  g_cuts <- left_join(cluster_y, tallied_reps, by = cnames) %>%
    unique() %>% mutate(across(dr, as.character))
  
  
  td_i <- epiCohesionV4(g_cuts, epi_melt, tr_dists) %>% 
    set_colnames(c(paste0("TP", tpx, "_", colnames(.))))
  colnames(td_i) %<>% gsub("ECC", paste0("ECC.0.", tau, ".", gamma), x = .)
  
  return(td_i)
}

### Incorporating the allele data with the epidemiological data
# strain_data <- selected_tp; tpx <- k; cluster_y <- cluster_x[dr %in% k_drs]
epiCollectionByClusterV1 <- function(strain_data, tau, gamma, transformed_dists, tpx, cluster_y) {
  # cat(paste0("\n   Collecting ECC values for temporal = ", tau, ", geo = ", gamma))

  # outputMessages(paste0("      Preparing table of distances: sqrt(", tau, "*(temp^2) + ", gamma, "*(geo^2))"))
  epi_melt <- transformed_dists %>%
    mutate(Total.Distinv = 1 - sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>%
    select(dr1, dr2, Total.Distinv) %>% set_colnames(c("dr_1", "dr_2", "value")) %>% as.data.table()
  invisible(gc())

  # outputMessages("      Incorporating the metadata with the clusters (typing data) ...")
  # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
  cx <- colnames(cluster_y)[1]
  tallied_reps <- cluster_y %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
  cnames <- intersect(colnames(tallied_reps), colnames(cluster_y))
  g_cuts <- left_join(cluster_y, tallied_reps, by = cnames) %>%
      unique() %>% mutate(across(dr, as.character))

  td_i <- epiCohesionV3(g_cuts, epi_melt) %>%
    set_colnames(c(paste0("TP", tpx, "_", colnames(.))))
  colnames(td_i) %<>% gsub("ECC", paste0("ECC.0.", tau, ".", gamma), x = .)

  return(td_i)
}

# tm <- matrix(data = 2, nrow = 3, ncol = 3)
# epiMelt(tm, tm, 1, 2)
# EPI-HELPER-MODULAR -----------------------------------------------------------------------------
# # maxd <- max(logdata) # mind <- min(logdata)
transformData2 <- function(dm, dtype, min_dist, max_dist) {
  logdata <- dm %>% add(10) %>% log10()
  x1 <- min_dist %>% add(10) %>% log10()
  x2 <- max_dist %>% add(10) %>% log10()
  if (dtype == "temp") {
    logdata[logdata == -Inf] <- 0
    if (x2 == 0) {
      logdata <- 0
    }else {
      logdata <- ((logdata - x1) / (x2 - x1))
    }
  }else if (dtype == "geo") {
    if(x2 == 1){
      logdata[1:nrow(logdata), 1:nrow(logdata)] <- 0
    } else {
      logdata <- ((logdata - x1) / (x2 - x1))
    }
  }
  return(logdata)
}

collectTransforms2 <- function(dms, extremes) {
  transformed_temp <- dms$temp %>% 
    transformTempDists(., extremes$mint, extremes$maxt) %>%
    set_rownames(., rownames(dms$temp)) %>% set_colnames(., colnames(dms$temp)) 
  
  transformed_geo <- dms$geo %>%
    transformGeoDists(., extremes$ming, extremes$maxg) %>%
    set_rownames(., rownames(dms$geo)) %>% set_colnames(., colnames(dms$geo))
  
  rm(dms); gc()
  
  return(list("temp" = transformed_temp, "geo" = transformed_geo))
}

collectTransforms <- function(dms, extremes) {
  transformed_temp <- dms$temp %>%
    transformData2(., "temp", extremes$mint, extremes$maxt) %>%
    formatData(., c("dr1","dr2","Temp.Dist"))

  transformed_geo <- dms$geo %>%
    transformData2(., "geo", extremes$ming, extremes$maxg) %>%
    formatData(., c("dr1","dr2","Geog.Dist"))

  rm(dms); gc()

  transformed_dists <- merge.data.table(transformed_temp, transformed_geo)

  rm(transformed_temp); rm(transformed_geo); gc()
  return(transformed_dists)
}

formatData <- function(dm, newnames) {
  dm %>% as.data.frame() %>% rownames_to_column("dr1") %>% as.data.table() %>% 
    melt.data.table(., id.vars = "dr1", variable.name = "dr2", value.name = newnames[3]) %>%
    # as_tibble() %>% 
    mutate(dr2 = as.character(dr2)) %>% 
    # as.data.table() %>% 
    set_colnames(newnames) %>% return()
}

processedStrains <- function(base_strains) {
  loc_cols <- intersect(c("Country", "Province", "City"), colnames(base_strains)) %>% sort()
  
  if (length(loc_cols) > 0) {
    strain_data <- base_strains %>% 
      mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")), 
             Location = do.call(paste, c(base_strains[loc_cols], sep = "_")))
    
    assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
      unique() %>% rownames_to_column("dr") %>% as.data.table()
    dr_matches <- left_join(strain_data, assignments, 
                            by = c("Latitude", "Longitude", "Date", "Location")) %>% 
      select(Strain, dr)
  }else {
    strain_data <- base_strains %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")))
    assignments <- strain_data %>% select(Date, Latitude, Longitude) %>% 
      unique() %>% rownames_to_column("dr") %>% as.data.table()
    dr_matches <- left_join(strain_data, assignments, 
                            by = c("Latitude", "Longitude", "Date")) %>% select(Strain, dr)
  }
  
  list("strain_data" = strain_data, 
       "assignments" = assignments, 
       "dr_matches" = dr_matches) %>% return()
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
# -----------------------------------------------------------------------------------------------------
checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}