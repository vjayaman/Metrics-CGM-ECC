# ---------------------------------------------------------------------------------------------------
# source("tmp.R")
# newECCS <- function(strains, source_file, sigma, tau, gamma, cpus, typing_data) {
#   strain_data <- read_tsv(strains) %>%
#     mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
#            Location = paste(Country, Province, City, sep = "_"))
# 
#   source_pw <- read_tsv(source_file) %>%
#     mutate(Source.Dist = 1- value) %>% select(-value)
# 
#   assignments <- strain_data %>% select(Source, Date, Latitude, Longitude, Location) %>%
#     unique() %>% rownames_to_column("SC")
#   assignments %<>% select(Date) %>% unique() %>% rownames_to_column("DSC") %>% left_join(assignments, .)
#   assignments %<>% select(Latitude, Longitude) %>% unique() %>% rownames_to_column("LSC") %>%
#     left_join(assignments, .) %>%
#     mutate(across(c(DSC, LSC), as.double))
# 
#   date_assignments <- assignments %>% select(DSC, Date) %>% unique()
#   dm <- date_assignments %>% distMatrix(., "temp", "Date") %>%
#     set_rownames(date_assignments$DSC) %>% set_colnames(date_assignments$DSC) %>%
#     transformData(., "temp")
#   dm2 <- dm %>% as.data.frame() %>% rownames_to_column("DSC1") %>% as.data.table()
# 
#   dm_temp <- melt.data.table(dm2, id.vars = "DSC1", variable.name = "DSC2", value.name = "Temp.Dist") %>%
#     as_tibble() %>% mutate(DSC2 = as.character(DSC2)) %>%
#     mutate(across(DSC1, as.numeric)) %>% mutate(across(DSC2, as.numeric)) %>% as.data.table()
# 
#   # Geographical distances
#   geo_assignments <- assignments %>% select(LSC, Latitude, Longitude) %>% unique()
#   dm_geo <- geo_assignments %>% distMatrix(., "geo", c("Latitude", "Longitude")) %>%
#     transformData(., "geo")
# 
#   dm_geo2 <- dm_geo %>% as.data.frame() %>% rownames_to_column("LSC1") %>% as.data.table()
#   dm_geo3 <- melt.data.table(dm_geo2, id.vars = "LSC1", variable.name = "LSC2", value.name = "Geog.Dist") %>%
#     as_tibble() %>% mutate(LSC2 = as.character(LSC2)) %>%
#     mutate(across(LSC1, as.numeric)) %>% mutate(across(LSC2, as.numeric)) %>% as.data.table()
# 
#   sc_matches <- strain_data %>% left_join(., assignments) %>%
#     select(Strain, SC)
# 
#   temp_formatted <- assignments %>% left_join(., dm_temp, by = c("DSC" = "DSC1")) %>%
#     select(SC, DSC2, Temp.Dist) %>% rename(SC1 = SC) %>%
#     left_join(assignments, ., by = c("DSC" = "DSC2")) %>%
#     rename(SC2 = SC) %>% select(SC1, SC2, Temp.Dist) %>% as.data.table()
# 
#   geo_formatted <- assignments %>% left_join(., dm_geo3, by = c("LSC" = "LSC1")) %>%
#     select(SC, LSC2, Geog.Dist) %>% rename(SC1 = SC) %>%
#     left_join(assignments, ., by = c("LSC" = "LSC2")) %>%
#     rename(SC2 = SC) %>% select(SC1, SC2, Geog.Dist) %>% as.data.table()
# 
#   matched <- merge.data.table(temp_formatted, geo_formatted) %>%
#     arrange(SC1, SC2) %>%
#     inner_join(sc_matches, ., by = c("SC" = "SC1")) %>%
#     rename(Strain.1 = Strain) %>% select(-SC) %>%
#     inner_join(sc_matches, ., by = c("SC" = "SC2")) %>%
#     rename(Strain.2 = Strain) %>% select(Strain.1, Strain.2, Temp.Dist, Geog.Dist)
# 
#   # This step brings the duplicated strain pairwise calculations back in
#   # unnecessary use of memory, esp since they're NA since we removed them earlier
#   epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, matched)# %>% filter(!is.na(Temp.Dist))
#   epi.matrix <- EpiMatrix(epi.table)
# 
#   # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
#   # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
#   eccs <- lapply(typing_data, function(typing_datum) {
#     g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
#     # Method with singletons merged separately, for speed/efficiency:
#     epi_cohesion_sep(g_cuts, epi.matrix, cpus)
#     # epi_cohes_cal_na(g_cuts, epi.matrix, cpus)
#     # Original method: epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
#   })
# 
#   a <- readRDS("eccs.Rds")
#   assert("Newly generated match the original", all(identical(eccs[[1]], a[[1]]), identical(eccs[[2]], a[[2]])))
# }

# ---------------------------------------------------------------------------------------------------
source("tmp.R")
# Note: dr stands for data representative
newECCS <- function(strains, source_file, sigma, tau, gamma, cpus, typing_data) {
  strain_data <- read_tsv(strains) %>%
    mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
           Location = paste(Country, Province, City, sep = "_"))
  
  source_pw <- read_tsv(source_file) %>%
    mutate(Source.Dist = 1- value) %>% select(-value)
  
  assignments <- strain_data %>% select(Source, Date, Latitude, Longitude, Location) %>%
    unique() %>% rownames_to_column("dr")
  assignments %<>% select(Date) %>% unique() %>% rownames_to_column("Ddr") %>% left_join(assignments, .)
  assignments %<>% select(Latitude, Longitude) %>% unique() %>% rownames_to_column("Ldr") %>%
    left_join(assignments, .) %>%
    mutate(across(c(Ddr, Ldr), as.double))
  
  # Temporal distances - all possible date pair distances
  date_assignments <- assignments %>% select(Ddr, Date) %>% unique()
  dm_temp <- pairwiseDists(date_assignments, "temp", "Date", c("Ddr1", "Ddr2", "Temp.Dist"))
  
  # Geographical distances - all possible lat-long pair distances
  geo_assignments <- assignments %>% select(Ldr, Latitude, Longitude) %>% unique()
  dm_geo <- pairwiseDists(geo_assignments, "geo", c("Latitude", "Longitude"), c("Ldr1", "Ldr2", "Geog.Dist"))
  
  temp_formatted <- assignments %>% 
    left_join(., dm_temp, by = c("Ddr" = "Ddr1")) %>%
    select(dr, Ddr2, Temp.Dist) %>% rename(dr1 = dr) %>%
    left_join(assignments, ., by = c("Ddr" = "Ddr2")) %>%
    rename(dr2 = dr) %>% select(dr1, dr2, Temp.Dist) %>% as.data.table()
  
  geo_formatted <- assignments %>% left_join(., dm_geo, by = c("Ldr" = "Ldr1")) %>%
    select(dr, Ldr2, Geog.Dist) %>% rename(dr1 = dr) %>%
    left_join(assignments, ., by = c("Ldr" = "Ldr2")) %>%
    rename(dr2 = dr) %>% select(dr1, dr2, Geog.Dist) %>% as.data.table()
  
  drs <- merge.data.table(temp_formatted, geo_formatted) %>% arrange(dr1, dr2)
  dr_matches <- strain_data %>% left_join(., assignments) %>% select(Strain, dr)

  # This step brings the duplicated strain pairwise calculations back in
  # unnecessary use of memory, esp since they're NA since we removed them earlier
  # epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, matched)
  epi.table <- drs %>% 
    mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) ), Epi.Sym = 1 - Total.Dist) %>% 
    select(dr1, dr2, Total.Dist) %>% as_tibble()
    # inner_join(dr_matches, ., by = c("dr" = "dr1")) %>%
    # rename(Strain.1 = Strain) %>% select(-dr) %>%
    # inner_join(dr_matches, ., by = c("dr" = "dr2")) %>%
    # rename(Strain.2 = Strain) %>% as_tibble() %>% select(Strain.1, Strain.2, Total.Dist)
  
  # epi.matrix <- EpiMatrix(epi.table)
  epi.matrix <- dcast(epi.table, formula = dr1 ~ dr2, value.var = "Total.Dist")
  epi.matrix <- as.matrix(epi.matrix[,2:ncol(epi.matrix)]) 
  rownames(epi.matrix) <- colnames(epi.matrix)

  drcounts <- dr_matches %>% select(dr) %>% group_by(dr) %>% count() %>% ungroup()
  
  new_td_1 <- td[[1]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
    left_join(., dr_matches) %>% 
    select(dr, T0) %>% left_join(., drcounts) %>% unique()
  
  new_td_2 <- td[[2]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
    left_join(., dr_matches) %>% 
    select(dr, T0) %>% left_join(., drcounts) %>% unique()
  
  # g_cuts <- new_td_2 %>% rename(genome = dr)
  # epi_matrix <- epi.matrix
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  eccs <- lapply(typing_data, function(typing_datum) {
    g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
    # Method with singletons merged separately, for speed/efficiency:
    epi_cohesion_sep(g_cuts, epi.matrix, cpus)
    # epi_cohes_cal_na(g_cuts, epi.matrix, cpus)
    # Original method: epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
  })

  a <- readRDS("eccs.Rds")
  assert("Newly generated match the original", all(identical(eccs[[1]], a[[1]]), identical(eccs[[2]], a[[2]])))
}

# -----------------------------------------------------------------------------------------------
# # source("tmp.R")
# newECCSRmDup <- function(strains, source_file, sigma, tau, gamma, cpus, typing_data) {
#   strain_data <- read_tsv(strains) %>%
#     mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
#            Location = paste(Country, Province, City, sep = "_"))
#   
#   source_pw <- read_tsv(source_file) %>%
#     mutate(Source.Dist = 1- value) %>% select(-value)
#   
#   assignments <- strain_data %>% select(Source, Date, Latitude, Longitude, Location) %>%
#     unique() %>% rownames_to_column("SC")
#   assignments %<>% select(Date) %>% unique() %>% rownames_to_column("DSC") %>% left_join(assignments, .)
#   assignments %<>% select(Latitude, Longitude) %>% unique() %>% rownames_to_column("LSC") %>% 
#     left_join(assignments, .) %>% mutate(across(c(DSC, LSC), as.double))
#   date_assignments <- assignments %>% select(DSC, Date) %>% unique()
#   geo_assignments <- assignments %>% select(LSC, Latitude, Longitude) %>% unique()
#   
#   
#   dm <- date_assignments %>% distMatrix(., "temp", "Date") %>%
#     set_rownames(date_assignments$DSC) %>% set_colnames(date_assignments$DSC) %>%
#     transformData(., "temp")
#   # dm[lower.tri(dm)] <- NA
#   dm2 <- dm %>% as.data.frame() %>% rownames_to_column("DSC1") %>% as.data.table()
#   
#   dm_temp <- melt.data.table(dm2, id.vars = "DSC1", variable.name = "DSC2", value.name = "Temp.Dist") %>%
#     as_tibble() %>% #na.omit(.) %>% 
#     mutate(DSC2 = as.character(DSC2)) %>%
#     mutate(across(DSC1, as.numeric)) %>% mutate(across(DSC2, as.numeric)) %>% as.data.table()
#   
#   # Geographical distances
#   dm_geo <- geo_assignments %>% distMatrix(., "geo", c("Latitude", "Longitude")) %>%
#     transformData(., "geo")
#   # dm_geo[lower.tri(dm_geo)] <- NA
#   
#   dm_geo2 <- dm_geo %>% as.data.frame() %>% rownames_to_column("LSC1") %>% as.data.table()
#   dm_geo3 <- melt.data.table(dm_geo2, id.vars = "LSC1", variable.name = "LSC2", value.name = "Geog.Dist") %>%
#     as_tibble() %>% #na.omit() %>% 
#     mutate(LSC2 = as.character(LSC2)) %>%
#     mutate(across(LSC1, as.numeric)) %>% mutate(across(LSC2, as.numeric)) %>% as.data.table()
#   
#   sc_matches <- strain_data %>% left_join(., assignments) %>%
#     select(Strain, SC)# %>% mutate(across(SC, as.integer)) %>% 
#   
#   
#   
#   temp_formatted <- assignments %>% left_join(., dm_temp, by = c("DSC" = "DSC1")) %>% 
#     select(SC, DSC2, Temp.Dist) %>% rename(SC1 = SC) %>% 
#     left_join(assignments, ., by = c("DSC" = "DSC2")) %>% 
#     rename(SC2 = SC) %>% 
#     select(SC1, SC2, Temp.Dist) %>% as.data.table()
#   
#   geo_formatted <- assignments %>% left_join(., dm_geo3, by = c("LSC" = "LSC1")) %>% 
#     select(SC, LSC2, Geog.Dist) %>% rename(SC1 = SC) %>% 
#     left_join(assignments, ., by = c("LSC" = "LSC2")) %>% 
#     rename(SC2 = SC) %>% 
#     select(SC1, SC2, Geog.Dist) %>% as.data.table()
#   
#   matched <- merge.data.table(temp_formatted, geo_formatted) %>% 
#     arrange(SC1, SC2) %>% 
#     inner_join(sc_matches, ., by = c("SC" = "SC1")) %>%
#     rename(Strain.1 = Strain) %>% select(-SC) %>%
#     inner_join(sc_matches, ., by = c("SC" = "SC2")) %>%
#     rename(Strain.2 = Strain) %>% select(Strain.1, Strain.2, Temp.Dist, Geog.Dist)
#   
#   # This part is only necessary if we remove the lower or upper triangle earlier - may not even work
#   # with_dups <- matched %>% filter(Temp.Dist != 0 & Geog.Dist != 0) %>%
#   #   select(Strain.2, Strain.1, Temp.Dist, Geog.Dist) %>%
#   #   rename(Strain.1 = Strain.2, Strain.2 = Strain.1) %>%
#   #   bind_rows(matched, .)
#   
#   # This step brings the duplicated strain pairwise calculations back in
#   # unnecessary use of memory, esp since they're NA since we removed them earlier
#   epi.table <- EpiTable2(strain_data, source_pw, sigma, tau, gamma, matched)
#   epi.matrix <- EpiMatrix(epi.table)
#   
#   # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
#   # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
#   eccs <- lapply(typing_data, function(typing_datum) {
#     g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
#     # Method with singletons merged separately, for speed/efficiency:
#     epi_cohesion_sep(g_cuts, epi.matrix, cpus)
#     # epi_cohes_cal_na(g_cuts, epi.matrix, cpus)
#     # Original method: epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
#   })
#   
#   # typing_datum <- typing_data[[1]]
#   # g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
#   # epi_matrix <- epi.matrix
#   
#   a <- readRDS("eccs.Rds")
#   assert("Newly generated match the original", all(identical(eccs[[1]], a[[1]]), identical(eccs[[2]], a[[2]])))
# }
# 
# # b <- readRDS("geog_temp.Rds") %>% mutate(across(c(Strain.1, Strain.2), as.character))
# # x <- b %>% left_join(matched, .)
# # assert("Duplicates removed early in new method", identical(x, matched))
# # # up to this point, verified
