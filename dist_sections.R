### Incorporating the allele data with the epidemiological data 
actualResults <- function() {
  source("zactual.R")
  library(beepr)
  # strain_data, source_file, sigma, tau, gamma, cpus, typing_data
  source_pw <- read_tsv(source_file) %>% mutate(Source.Dist = 1- value) %>% select(-value)
  
  # look into a method that runs the distances when given two separate vectors (then we can 
  # avoid the extra work that comes with overlapping)
  temp_pw <- generateDistances(strain_data, "temp", "Date", "Temp.Dist")
  geog_pw <- generateDistances(strain_data, "geo", c("Latitude", "Longitude"), "Geog.Dist")
  
  geog_temp_tr <- left_join(geog_pw$transformed, temp_pw$transformed, by = c("Strain.1", "Strain.2"))
  
  ## This generates a table of your comparisons based on your epidemiological data (source, time, geographical) 
  ## with the assigned weights of s, t and g and then computes the similarity/distance and generates a matrix
  epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, geog_temp_tr)
  # saveRDS(epi.table, "epitable.Rds")
  epi.matrix <- EpiMatrix(epi.table)
  
  typing_datum <- typing_data[[2]]
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  # eccs <- lapply(typing_data, function(typing_datum) {
  g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
  # Method with singletons merged separately, for speed/efficiency:
  eccs <- epi_cohesion_sep(g_cuts, epi.matrix, cpus) %>% as.data.table()
  beep(3)
  # Original method: epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
  # })
  
  return(eccs)
}

# ---------------------------------------------------------------------------------------------------
# Note: dr stands for data representative
testResults <- function() {
  source("ztest.R")
  library(beepr)
  # strain_data, sigma, tau, gamma, cpus, typing_data

  # in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
  assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
    unique() %>% rownames_to_column("dr")
  
  # Temporal distances - all possible date pair distances
  dm_temp <- assignments %>% select(dr, Date) %>% pairwiseDists(., "temp", "Date", c("dr1", "dr2", "Temp.Dist"))
  
  # Geographical distances - all possible lat-long pair distances
  dm_geo <- assignments %>% select(dr, Latitude, Longitude) %>% 
    pairwiseDists(., "geo", c("Latitude", "Longitude"), c("dr1", "dr2", "Geog.Dist"))
  
  epi_table <- merge.data.table(dm_temp, dm_geo) %>% 
    mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
    select(dr1, dr2, Total.Dist) %>% as_tibble()
  
  epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
  epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)]) 
  rownames(epi_matrix) <- colnames(epi_matrix)

  # create similarity values from epi distance matrix:
  epi_melt_all <- melt(as.matrix(1-epi_matrix)) %>%
    mutate(across(c(Var1, Var2), as.character)) %>% as.data.table()
  
  # Identifying which strains match with which non-redundant data representatives
  dr_matches <- strain_data %>% left_join(., assignments) %>% select(Strain, dr)
  beep(3)
  # From this part on we're testing for a single cluster - cluster 1 at TP1 (with 202 strains, but 122 drs)
  
  # dr_td1 <- td[[1]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
  dr_td1 <- td[[2]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
    left_join(., dr_matches) %>% 
    mutate(across(dr, as.character)) %>% select(-Strain)
  
  # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
  g_cuts <- dr_td1 %>% group_by(T0) %>% count(dr) %>% ungroup() %>% 
    left_join(dr_td1, .) %>% unique() %>% mutate(across(dr, as.character))
  # epi_melt <- epi_melt_all %>% filter(Var1 %in% g_cuts$dr, Var2 %in% g_cuts$dr)
  epi_melt <- epi_melt_all
  
  z2 <- epi_cohesion_new(g_cuts, epi_melt)
  beep(3)
  return(z2)
}
