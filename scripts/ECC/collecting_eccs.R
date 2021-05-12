
### Incorporating the allele data with the epidemiological data 
epiCollection <- function(strain_data, tau, gamma, typing_data) {
  cat(paste0("\nCollecting ECC values for temporal = ", tau, ", geo = ", gamma))
  
  # Note: dr stands for data representative
  # in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
  assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
    unique() %>% rownames_to_column("dr")
    
  # Temporal distances - all possible date pair distances
  dm_temp <- assignments %>% select(dr, Date) %>% 
    pairwiseDists(., "temp", "Date", c("dr1", "dr2", "Temp.Dist"))
    
  # Geographical distances - all possible lat-long pair distances
  dm_geo <- assignments %>% select(dr, Latitude, Longitude) %>% 
    pairwiseDists(., "geo", c("Latitude", "Longitude"), c("dr1", "dr2", "Geog.Dist"))
    
  epi_table <- merge.data.table(dm_temp, dm_geo) %>% 
    mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
    select(dr1, dr2, Total.Dist) %>% as.data.table()
    
  rm(dm_geo)
  rm(dm_temp)
  
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
    
  gc()
  
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
    
    return(td_i)
  })
  
  return(eccs)
}

# basicAverages <- function(df, tp, avg_raw) {
#   df %>% select("Strain", all_of(tp)) %>% 
#     left_join(avg_raw, by = "Strain") %>% 
#     arrange({{tp}}) %>% 
#     group_by({{tp}}) %>% 
#     mutate(avg_date = mean(Date), avg_lat = mean(Latitude), avg_long = mean(Longitude)) %>% 
#     ungroup() %>% 
#     select(Strain, all_of(tp), grep("avg", colnames(.), value = TRUE)) %>% 
#     set_colnames(gsub("avg", paste0(as.character(tp), "_avg"), colnames(.))) %>% return()
# }

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


