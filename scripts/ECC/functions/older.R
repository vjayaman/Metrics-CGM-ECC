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


### Incorporating the allele data with the epidemiological data
epiCollection <- function(strain_data, tau, gamma, typing_data, transformed_dists,
                          dm_temp, dm_geo, dr_matches) {
  cat(paste0("\n   Collecting ECC values for temporal = ", tau, ", geo = ", gamma))
  
  outputMessages(paste0("      Preparing table of distances: sqrt(", tau, "*(temp^2) + ", gamma, "*(geo^2))"))
  epi_table <- transformed_dists %>%
    mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>%
    select(dr1, dr2, Total.Dist) %>% as.data.table()
  
  outputMessages("      Generating matrix of similarity values from epi distance matrix ...")
  epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
  epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)])
  rownames(epi_matrix) <- colnames(epi_matrix)
  
  # create similarity values from epi distance matrix:
  epi_melt <- as.matrix(1-epi_matrix) %>%
    as.data.table(., keep.rownames = TRUE) %>%
    melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>%
    set_colnames(c("dr_1", "dr_2", "value")) %>% as.data.table()
  
  rm(epi_table)
  rm(epi_matrix)
  invisible(gc())
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  eccs <- lapply(1:length(typing_data), function(i) {
    outputMessages(paste0("      Calculating ECCs for timepoint ", i, ", which has ",
                          nrow(typing_data[[i]]), " strains to consider ..."))
    dr_td1 <- typing_data[[i]] %>% rownames_to_column("Strain") %>% as_tibble() %>%
      left_join(., dr_matches, by = "Strain") %>%
      mutate(across(dr, as.character)) %>% select(-Strain)
    
    # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
    cx <- colnames(dr_td1)[1]
    tallied_reps <- dr_td1 %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
    g_cuts <- left_join(dr_td1, tallied_reps, by = intersect(colnames(tallied_reps), colnames(dr_td1))) %>%
      unique() %>% mutate(across(dr, as.character))
    
    td_i <- epi_cohesion_new(g_cuts, epi_melt) %>%
      set_colnames(c(paste0("TP", i, "_", colnames(.))))
    colnames(td_i) %<>% gsub("ECC", paste0("ECC.0.", tau, ".", gamma), x = .)
    
    return(td_i)
  })
  
  return(eccs)
}

# mergeECCs <- function(eccs, tpx, typing_data) {
#   tbl1 <- as.data.table(typing_data)
#   sapply(eccs, "[", tpx) %>% Reduce(function(...) merge(...), .) %>% as.data.table() %>%
#     merge.data.table(tbl1, ., by = intersect(colnames(tbl1), colnames(.))) %>% return()
# }