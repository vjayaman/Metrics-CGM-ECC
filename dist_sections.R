# ---------------------------------------------------------------------------------------------------
source("ztest.R")
# Note: dr stands for data representative
newECCS <- function(strain_data, sigma, tau, gamma, cpus, typing_data) {

  # in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
  assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
    unique() %>% rownames_to_column("dr")
  
  # Temporal distances - all possible date pair distances
  dm_temp <- assignments %>% select(dr, Date) %>% 
    pairwiseDists(., "temp", "Date", c("dr1", "dr2", "Temp.Dist"))
  
  # Geographical distances - all possible lat-long pair distances
  dm_geo <- assignments %>% select(dr, Latitude, Longitude) %>% 
    pairwiseDists(., "geo", c("Latitude", "Longitude"), c("dr1", "dr2", "Geog.Dist"))
  
  drs <- merge.data.table(dm_temp, dm_geo)
  
  # Total distances - base for ECCs: 
  #     total distance = sqrt( (temporal coeff)*(temporal distances)^2 + (geo coeff)*(geographical distances)^2)
  # removed (for now): Epi.Sym = 1 - Total.Dist
  epi_table <- drs %>% 
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

  # From this part on we're testing for a single cluster - cluster 1 at TP1 (with 202 strains, but 122 drs)
  dr_td1 <- td[[1]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
    left_join(., dr_matches) %>% mutate(across(dr, as.numeric)) %>% 
    filter(T0 == 1) %>% select(-Strain)

  # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
  g_cuts <- dr_td1 %>% count(dr) %>% left_join(dr_td1, .) %>% unique() %>% mutate(across(dr, as.character))

  # # cluster1 <- epi_table %>% filter(dr1 %in% tmp$dr, dr2 %in% tmp$dr)
  # x1 <- td[[2]] %>% rownames_to_column("Strain") %>% as_tibble() %>% filter(T0 == 1)
  # cluster1 <- dr_matches %>% filter(Strain %in% x1$Strain)
  # c1 <- cluster1 %>% select(dr) %>% group_by(dr) %>% count() %>% ungroup()
  # d1 <- epi_table %>% filter(dr1 %in% c1$dr, dr2 %in% c1$dr)
  # 
  # e1 <- left_join(d1, c1, by = c("dr1" = "dr")) %>% rename(n1 = n) %>%
  #   left_join(., c1, by = c("dr2" = "dr")) %>% rename(n2 = n) %>%
  #   mutate(npairs = n1 * n2, Total.Sim = 1 - Total.Dist) %>% mutate(value = npairs * Total.Sim)
  
  epi_melt <- epi_melt_all %>% filter(Var1 %in% g_cuts$dr, Var2 %in% g_cuts$dr)
  
  
  z2 <- epi_cohesion_new(g_cuts, epi_melt)
  # epi_cohes_cal_na(g_cuts, epi_matrix, cpus)
  # Original method: epi_cohesion_calc(g_cuts, epi_matrix, cpus = cpus)
  # })
  
  a <- readRDS("eccs.Rds")
  assert("Newly generated match the original", all(identical(eccs[[1]], a[[1]]), identical(eccs[[2]], a[[2]])))
}

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
# source("tmp.R")
# # Note: dr stands for data representative
# newECCS <- function(strain_data, sigma, tau, gamma, cpus, typing_data) {
#   # strain_data <- read_tsv(strains) %>%
#   #   mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
#   #          Location = paste(Country, Province, City, sep = "_"))
#   
#   # source_pw <- read_tsv(source_file) %>%
#   #   mutate(Source.Dist = 1- value) %>% select(-value)
#   
#   assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
#     unique() %>% rownames_to_column("dr")
#   
#   assignments %<>% select(Date) %>% unique() %>% rownames_to_column("Ddr") %>% left_join(assignments, .)
#   assignments %<>% select(Latitude, Longitude) %>% unique() %>% rownames_to_column("Ldr") %>%
#     left_join(assignments, .) %>%
#     mutate(across(c(Ddr, Ldr), as.double))
#   
#   # Temporal distances - all possible date pair distances
#   date_assignments <- assignments %>% select(Ddr, Date) %>% unique()
#   dm_temp <- pairwiseDists(date_assignments, "temp", "Date", c("Ddr1", "Ddr2", "Temp.Dist"))
#   
#   # Geographical distances - all possible lat-long pair distances
#   geo_assignments <- assignments %>% select(Ldr, Latitude, Longitude) %>% unique()
#   dm_geo <- pairwiseDists(geo_assignments, "geo", c("Latitude", "Longitude"), c("Ldr1", "Ldr2", "Geog.Dist"))
#   
#   temp_formatted <- assignments %>% 
#     left_join(., dm_temp, by = c("Ddr" = "Ddr1")) %>%
#     select(dr, Ddr2, Temp.Dist) %>% rename(dr1 = dr) %>%
#     left_join(assignments, ., by = c("Ddr" = "Ddr2")) %>%
#     rename(dr2 = dr) %>% select(dr1, dr2, Temp.Dist) %>% as.data.table()
#   
#   geo_formatted <- assignments %>% left_join(., dm_geo, by = c("Ldr" = "Ldr1")) %>%
#     select(dr, Ldr2, Geog.Dist) %>% rename(dr1 = dr) %>%
#     left_join(assignments, ., by = c("Ldr" = "Ldr2")) %>%
#     rename(dr2 = dr) %>% select(dr1, dr2, Geog.Dist) %>% as.data.table()
#   
#   drs <- merge.data.table(temp_formatted, geo_formatted) %>% arrange(dr1, dr2)
# 
#   # This step brings the duplicated strain pairwise calculations back in
#   # unnecessary use of memory, esp since they're NA since we removed them earlier
#   # epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, matched)
#   epi.table <- drs %>% 
#     mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) ), Epi.Sym = 1 - Total.Dist) %>% 
#     select(dr1, dr2, Total.Dist) %>% as_tibble()
#     # inner_join(dr_matches, ., by = c("dr" = "dr1")) %>%
#     # rename(Strain.1 = Strain) %>% select(-dr) %>%
#     # inner_join(dr_matches, ., by = c("dr" = "dr2")) %>%
#     # rename(Strain.2 = Strain) %>% as_tibble() %>% select(Strain.1, Strain.2, Total.Dist)
#   
#   # epi.matrix <- EpiMatrix(epi.table)
#   epi.matrix <- dcast(epi.table, formula = dr1 ~ dr2, value.var = "Total.Dist")
#   epi.matrix <- as.matrix(epi.matrix[,2:ncol(epi.matrix)]) 
#   rownames(epi.matrix) <- colnames(epi.matrix)
# 
#   dr_matches <- strain_data %>% left_join(., assignments) %>% select(Strain, dr)
#   drcounts <- dr_matches %>% select(dr) %>% group_by(dr) %>% count() %>% ungroup()
#   
#   dr_td1 <- td[[1]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
#     left_join(., dr_matches) %>% select(-Strain) %>% group_by(T0) %>% count(dr) %>% ungroup()
#   
#   dr_td2 <- td[[2]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
#     left_join(., dr_matches) %>% select(-Strain) %>% group_by(T0) %>% count(dr) %>% ungroup()
#     # left_join(., dr_matches) %>% select(dr, T0) %>% left_join(., drcounts) %>% unique()
#   
#   g_cuts <- dr_td1 %>% rename(genome = dr)
#   epi_matrix <- epi.matrix
#   
#   # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
#   # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
#   # eccs <- lapply(typing_data, function(typing_datum) {
#   #   g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
#     # Method with singletons merged separately, for speed/efficiency:
#   # create similarity values from epi distance matrix:
#   epi_melt <- melt(as.matrix(1-epi_matrix)) %>%
#     mutate(across(c(Var1, Var2), as.character)) %>% as.data.table()
#   
#   z2 <- epi_cohesion_new(g_cuts, epi_melt)
#     # epi_cohes_cal_na(g_cuts, epi.matrix, cpus)
#     # Original method: epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
#   # })
# 
#   a <- readRDS("eccs.Rds")
#   assert("Newly generated match the original", all(identical(eccs[[1]], a[[1]]), identical(eccs[[2]], a[[2]])))
# }
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
