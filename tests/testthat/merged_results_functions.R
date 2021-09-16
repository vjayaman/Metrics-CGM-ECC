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

# Given a dataframe df, two column names c1 and c2 (height and cluster respectively) and a new
# ID prefix tpx (e.g. "tp1"), creates an ID column and adds to df before returning df
newID <- function(df, tpx, c1, c2, ph, pc) {
  newh <- df %>% pull(c1) %>% as.character() %>% as.integer() %>% 
    formatC(., width = max(3, ph), format = "d", flag = "0") %>% paste0("h", .)
  newc <- df %>% pull(c2) %>% as.character() %>% as.integer() %>% 
    formatC(., width = max(3, pc), format = "d", flag = "0") %>% paste0("c", .)
  df %>% add_column(id = paste0(toupper(tpx), "_", newh, "_", newc)) %>% return()
}

manualEpiMelt <- function(df, extremes) {
  geo_dists <- df %>% select(Longitude, Latitude) %>% 
    as.data.frame() %>% earth.dist(dist = TRUE) %>% as.matrix() %>% 
    set_rownames(df$Strain) %>% set_colnames(df$Strain) %>% 
    add(10) %>% log10()
  g1 <- extremes$ming %>% add(10) %>% log10()
  g2 <- extremes$maxg %>% add(10) %>% log10()
  if(g2 == 1){
    geo_dists[1:nrow(geo_dists), 1:nrow(geo_dists)] <- 0
  } else {
    geo_dists <- ((geo_dists - g1) / (g2 - g1))
  }
  
  temp_dists <- df %>% pull(Date) %>% 
    dist(diag = FALSE, upper = FALSE, method = "euclidean") %>% as.matrix() %>% 
    set_rownames(df$Strain) %>% set_colnames(df$Strain) %>% add(10) %>% log10()
  t1 <- extremes$mint %>% add(10) %>% log10()
  t2 <- extremes$maxt %>% add(10) %>% log10()
  
  temp_dists[temp_dists == -Inf] <- 0
  if (t2 == 0) {
    temp_dists <- 0
  }else {
    temp_dists <- ((temp_dists - t1) / (t2 - t1))
  }
  
  temp_dists <- temp_dists %>% as.data.frame() %>% rownames_to_column("Strain.1") %>% 
    as.data.table() %>% melt.data.table(., id.vars = "Strain.1", variable.name = "Strain.2", 
                                        value.name = "Temp.Dist")
  
  geo_dists <- geo_dists %>% as.data.frame() %>% rownames_to_column("Strain.1") %>% 
    as.data.table() %>% melt.data.table(., id.vars = "Strain.1", variable.name = "Strain.2", 
                                        value.name = "Geog.Dist")
  
  merge.data.table(temp_dists, geo_dists) %>% return()
}