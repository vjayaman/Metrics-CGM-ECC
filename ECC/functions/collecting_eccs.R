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
  saveRDS(temp_pw, "inputs/ecc_averages/temp.Rds")
  
  geog_pw <- geog_calc(strain_data)
  saveRDS(geog_pw, "inputs/ecc_averages/geo.Rds")
  
  geog_temp <- left_join(geog_pw, temp_pw, by = c("Strain.1", "Strain.2"))
  
  ## This generates a table of your comparisons based on your epidemiological data (source, time, geographical) 
  ## with the assigned weights of s, t and g and then computes the similarity/distance and generates a matrix
  epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, geog_temp)
  epi.matrix <- EpiMatrix(epi.table)
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  eccs <- lapply(typing_data, function(typing_datum) {
    g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
    epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
  })
  
  # names(eccs) == c("TP1_T0", "TP2_T0")
  newnames <- sapply(strsplit(names(eccs), "_"), `[`, 1)
  
  lapply(names(eccs), function(x) {
    y <- eccs[[x]]
    td[[x]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
      left_join(., y, by = colnames(y)[1]) %>% 
      set_colnames(c("Strain", paste0(unlist(strsplit(x, "_"))[1], "_", colnames(y)))) %>% 
      set_colnames(gsub("ECC", paste0("ECC_", sigma, ".", tau, ".", gamma), colnames(.)))
    
  }) %>% set_names(newnames) %>% 
    Reduce(function(...) right_join(..., by = "Strain"), .) %>% return()
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

# z_tp1_values <- Reduce(function(...) merge(...), sapply(z, `[`, tp1$name)) %>% as_tibble()

