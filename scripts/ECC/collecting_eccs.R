### Incorporating the allele data with the epidemiological data 
epiCollection <- function(strain_data, tau, gamma, typing_data) {
  cat(paste0("\nCollecting ECC values for temporal = ", tau, ", geo = ", gamma))
  
  # Note: dr stands for data representative
  # in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
  assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
    unique() %>% rownames_to_column("dr")
    
  # Temporal distances - all possible date pair distances
  # pairwiseDists <- function(dm, type, cnames, newnames) {
  dm_temp <- assignments %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
  # formatted_temp <- dm_temp %>% formatData(., c("dr1", "dr2", "Temp.Dist"))
  transformed_temp <- dm_temp %>% pairwiseDists(., "temp", c("dr1", "dr2", "Temp.Dist"))
    
  # Geographical distances - all possible lat-long pair distances
  dm_geo <- assignments %>% select(dr, Latitude, Longitude) %>% distMatrix(., "geo", c("Latitude", "Longitude"))
  # formatted_geo <- dm_geo %>% formatData(., c("dr1", "dr2", "Geog.Dist"))
  transformed_geo <- dm_geo %>% pairwiseDists(., "geo", c("dr1", "dr2", "Geog.Dist"))
  
  # actual_dists <- merge.data.table(formatted_temp, formatted_geo)
  
  transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
  epi_table <- transformed_dists %>% 
    mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
    select(dr1, dr2, Total.Dist) %>% as.data.table()
    
  rm(transformed_geo)
  rm(transformed_temp)
  rm(transformed_dists)
  
  epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
  epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)]) 
  rownames(epi_matrix) <- colnames(epi_matrix)
  
  # create similarity values from epi distance matrix:
  epi_melt <- as.matrix(1-epi_matrix) %>% 
    as.data.table(., keep.rownames = TRUE) %>% 
    melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>% 
    set_colnames(c("Var1", "Var2", "value")) %>% as.data.table()
  
  rm(epi_table)
  rm(epi_matrix)
  
  # Identifying which strains match with which non-redundant data representatives
  dr_matches <- strain_data %>% 
    left_join(., assignments, by = c("Latitude", "Longitude", "Date", "Location")) %>% 
    select(Strain, dr)
    
  invisible(gc())
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  eccs <- lapply(1:length(typing_data), function(i) {
    
    dr_td1 <- typing_data[[i]] %>% rownames_to_column("Strain") %>% as_tibble() %>%
      left_join(., dr_matches, by = "Strain") %>%
      mutate(across(dr, as.character)) %>% select(-Strain)

    # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
    tallied_reps <- dr_td1 %>% group_by(T0) %>% count(dr) %>% ungroup()
    g_cuts <- left_join(dr_td1, tallied_reps, by = intersect(colnames(tallied_reps), colnames(dr_td1))) %>%
      unique() %>% mutate(across(dr, as.character))
    
    td_i <- epi_cohesion_new(g_cuts, epi_melt) %>% 
      set_colnames(c(paste0("TP", i, "_", colnames(.))))
    colnames(td_i) %<>% gsub("ECC", paste0("ECC.", tau, ".", gamma), x = .)
    
    a2 <- avgDists(g_cuts, dm_temp, "Temp.Dist", colnames(td_i)[1])
    b2 <- avgDists(g_cuts, dm_geo, "Geog.Dist", colnames(td_i)[1])
    
    left_join(td_i, a2, by = intersect(colnames(td_i), colnames(a2))) %>% 
      left_join(., b2, by = intersect(colnames(.), colnames(b2))) %>% return()
  })
  
  return(eccs)
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

# Indicates length of a process in hours, minutes, and seconds, when given a name of the process 
# ("pt") and a two-element named vector with Sys.time() values named "start_time" and "end_time"
timeTaken <- function(pt, sw) {
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