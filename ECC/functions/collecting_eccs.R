
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
           Location = paste(Country, Province, City, sep = "_")
    )
  
  source_pw <- read_tsv(source_file) %>% 
    mutate(Source.Dist = 1- value) %>% select(-value)
  
  temp_pw <- temp_calc(strain_data)
  
  geog_pw <- geog_calc(strain_data)
  
  geog_temp <- left_join(geog_pw, temp_pw, by = c("Strain.1", "Strain.2"))
  
  ## This generates a table of your comparisons based on your epidemiological data (source, time, geographical) 
  ## with the assigned weights of s, t and g and then computes the similarity/distance and generates a matrix
  epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, geog_temp)
  epi.matrix <- EpiMatrix(epi.table)
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  
  # # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  # eccs <- lapply(typing_data, function(typing_datum) {
  #   g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
  #   epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
  # })
  
  eccs <- lapply(typing_data, function(typing_datum) {
    g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
    epi_cohesion_sep(g_cuts, epi.matrix, cpus = cpus)
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
  
  # results <- geog_temp %>% 
  #   mutate(across(c(Strain.1, Strain.2), as.character)) %>% 
  #   eccAverages(., ecc_data)
  
  avg_raw <- strain_data %>% select(Strain, Date, Longitude, Latitude)
  cnames <- colnames(ecc_data) %>% grep("Size|ECC", ., value = TRUE, invert = TRUE)
  df <- ecc_data %>% select(all_of(cnames))
  
  i <- grep("TP1", cnames, value = TRUE)
  avg_interim_tp1 <- df %>% rename("TPx" = i) %>% basicAverages(., i, avg_raw)
  
  j <- grep("TP2", cnames, value = TRUE)
  avg_interim_tp2 <- df %>% rename("TPx" = all_of(j)) %>% basicAverages(., j, avg_raw)
  
  results <- left_join(ecc_data, avg_interim_tp1) %>% left_join(., avg_interim_tp2)
  
  return(results)
}

basicAverages <- function(df, tpname, avg_raw) {
  df %>% select("Strain", "TPx") %>% 
    left_join(avg_raw, by = "Strain") %>% arrange("TPx") %>% group_by(TPx) %>% 
    mutate(TPx_avg_date = mean(Date), TPx_avg_lat = mean(Latitude), TPx_avg_long = mean(Longitude)) %>% 
    select(Strain, grep("TPx", colnames(.), value = TRUE)) %>% ungroup() %>% 
    set_colnames(gsub("TPx", tpname, colnames(.))) %>% return()
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

# z_tp1_values <- Reduce(function(...) merge(...), sapply(z, `[`, tp1$name)) %>% as_tibble()

