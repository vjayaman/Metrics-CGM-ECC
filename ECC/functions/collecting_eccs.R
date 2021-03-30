
checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

### Section 3: Incorporating the allele data with the epidemiological data 
oneCombo <- function(strains, source_file, sigma, tau, gamma, cpus, typing_data) {
  cat(paste0("\nCollecting ECC values for source = ", sigma, ", temporal = ", tau, ", geo = ", gamma))
  ### Section 2: Generating the EpiMatrix
  strain_data <- 
    read_tsv(strains) %>% 
    mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
           Location = paste(Country, Province, City, sep = "_"))
  
  source_pw <- read_tsv(source_file) %>% 
    mutate(Source.Dist = 1- value) %>% select(-value)
  
  a1 <- distMatrix(strain_data, "temp", "Date")
  a2 <- transformData(a1, "temp")
  temp_pw <- formatMatrix(strain_data, a2, "Temp.Dist")
  
  # assert("Modularized method and original method produce same result for temp data", 
  #        identical(temp_calc(strain_data), temp_pw))
  
  b1 <- distMatrix(strain_data, "geo", c("Latitude", "Longitude"))
  b2 <- transformData(b1, "geo")
  geog_pw <- formatMatrix(strain_data, b2, "Geog.Dist")
  
  # assert("Modularized method and original method produce same result for geo data", 
  #        identical(geog_calc(strain_data), geog_pw))
  
  geog_temp <- left_join(geog_pw, temp_pw, by = c("Strain.1", "Strain.2"))
  
  ## This generates a table of your comparisons based on your epidemiological data (source, time, geographical) 
  ## with the assigned weights of s, t and g and then computes the similarity/distance and generates a matrix
  epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, geog_temp)
  epi.matrix <- EpiMatrix(epi.table)
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  eccs <- lapply(typing_data, function(typing_datum) {
    g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
    # Method with singletons merged separately, for speed/efficiency
    epi_cohesion_sep(g_cuts, epi.matrix, cpus = cpus)
    # Original method:
    #   epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
  })

  # names(eccs) == c("TP1_T0", "TP2_T0")
  newnames <- sapply(strsplit(names(eccs), "_"), `[`, 1)
  
  ecc_data <- lapply(names(eccs), function(x) {
    y <- eccs[[x]]
    td[[x]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
      left_join(., y, by = colnames(y)[1]) %>% 
      set_colnames(c("Strain", paste0(unlist(strsplit(x, "_"))[1], "_", colnames(y)))) %>% 
      set_colnames(gsub("ECC", paste0("ECC_", sigma, ".", tau, ".", gamma), colnames(.)))
    
  }) %>% set_names(newnames) %>% 
    Reduce(function(...) right_join(..., by = "Strain"), .)
  
  # Adding basic average columns (Date, Longitude, Latitude)
  avg_raw <- strain_data %>% select(Strain, Date, Longitude, Latitude)
  cnames <- colnames(ecc_data) %>% grep("Size|ECC", ., value = TRUE, invert = TRUE)
  df <- ecc_data %>% select(all_of(cnames)) %>% set_colnames(c("Strain", "TP1", "TP2"))
  
  avg_tp1 <- df %>% basicAverages(., as.name("TP1"), avg_raw) %>% arrange(Strain)
  avg_tp2 <- df %>% basicAverages(., as.name("TP2"), avg_raw) %>% arrange(Strain)
  
  results <- left_join(ecc_data, avg_tp1) %>% left_join(., avg_tp2)
  
  # Adding distance averages
  # results <- geog_temp %>% 
  #   mutate(across(c(Strain.1, Strain.2), as.character)) %>% 
  #   eccAverages(., ecc_data)
  
  return(results)
}

basicAverages <- function(df, tp, avg_raw) {
  df %>% select("Strain", all_of(tp)) %>% 
    left_join(avg_raw, by = "Strain") %>% 
    arrange({{tp}}) %>% 
    group_by({{tp}}) %>% 
    mutate(avg_date = mean(Date), avg_lat = mean(Latitude), avg_long = mean(Longitude)) %>% 
    ungroup() %>% 
    select(Strain, all_of(tp), grep("avg", colnames(.), value = TRUE)) %>% 
    set_colnames(gsub("avg", paste0(as.character(tp), "_avg"), colnames(.))) %>% return()
}

# eccAverages <- function(pairwise_dists, ecc_data) {
#   
#   cnames <- colnames(ecc_data) %>% grep("Size|ECC", ., value = TRUE, invert = TRUE)
#   
#   tp1_base <- ecc_data %>% select(grep("Strain|TP1", cnames, value = TRUE)) %>% set_colnames(c("Strain", "TP1"))
#   tp1_data <- pairwise_dists %>% 
#     left_join(., tp1_base, by = c("Strain.1" = "Strain")) %>% rename(first_in_cl = TP1) %>% 
#     left_join(., tp1_base, by = c("Strain.2" = "Strain")) %>% rename(second_in_cl = TP1) %>% 
#     filter(first_in_cl == second_in_cl) %>% select(-second_in_cl) %>% rename(TP1 = first_in_cl)
#   
#   tp1_avg_eccs <- lapply(unique(tp1_data$TP1), function(x) {
#     tp1_data %>% filter(TP1 == x) %>% mutate(avg_geo = mean(Geog.Dist), avg_temp = mean(Temp.Dist))
#   }) %>% bind_rows() %>% 
#     select(TP1, avg_geo, avg_temp) %>% unique() %>% 
#     set_colnames(c(grep("TP1", cnames, value = TRUE), "TP1_avg_geo_ECC", "TP1_avg_temp_ECC"))
#   
#   
#   tp2_base <- ecc_data %>% select(grep("Strain|TP2", cnames, value = TRUE)) %>% set_colnames(c("Strain", "TP2"))
#   tp2_data <- pairwise_dists %>% 
#     left_join(., tp2_base, by = c("Strain.1" = "Strain")) %>% rename(first_in_cl = TP2) %>% 
#     left_join(., tp2_base, by = c("Strain.2" = "Strain")) %>% rename(second_in_cl = TP2) %>% 
#     filter(first_in_cl == second_in_cl) %>% select(-second_in_cl) %>% rename(TP2 = first_in_cl)
#   
#   tp2_avg_eccs <- lapply(unique(tp2_data$TP2), function(x) {
#     tp2_data %>% filter(TP2 == x) %>% mutate(avg_geo = mean(Geog.Dist), avg_temp = mean(Temp.Dist))
#   }) %>% bind_rows() %>% 
#     select(TP2, avg_geo, avg_temp) %>% unique() %>% 
#     set_colnames(c(grep("TP2", cnames, value = TRUE), "TP2_avg_geo_ECC", "TP2_avg_temp_ECC"))
#   
#   left_join(ecc_data, tp1_avg_eccs) %>% left_join(., tp2_avg_eccs) %>% return()
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


