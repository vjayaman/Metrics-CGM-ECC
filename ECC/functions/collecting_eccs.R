
checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

### Incorporating the allele data with the epidemiological data 
oneCombo <- function(strains, source_file, sigma, tau, gamma, cpus, typing_data) {
  cat(paste0("\nCollecting ECC values for source = ", sigma, ", temporal = ", tau, ", geo = ", gamma))
  ### Generating the EpiMatrix
  strain_data <- read_tsv(strains) %>% 
    mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
           Location = paste(Country, Province, City, sep = "_"))
  
  source_pw <- read_tsv(source_file) %>% 
    mutate(Source.Dist = 1- value) %>% select(-value)
  
  temp_pw <- generateDistances(strain_data, "temp", "Date", "Temp.Dist")
  geog_pw <- generateDistances(strain_data, "geo", c("Latitude", "Longitude"), "Geog.Dist")

  geog_temp_tr <- left_join(geog_pw$transformed, temp_pw$transformed, by = c("Strain.1", "Strain.2"))
  
  ## This generates a table of your comparisons based on your epidemiological data (source, time, geographical) 
  ## with the assigned weights of s, t and g and then computes the similarity/distance and generates a matrix
  epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, geog_temp_tr)
  epi.matrix <- EpiMatrix(epi.table)
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  eccs <- lapply(typing_data, function(typing_datum) {
    g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
    # Method with singletons merged separately, for speed/efficiency:
    epi_cohesion_sep(g_cuts, epi.matrix, cpus = cpus)
    # Original method: epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
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
  dist_avgs <- left_join(geog_pw$raw, temp_pw$raw, by = c("Strain.1", "Strain.2")) %>% 
    mutate(across(c(Strain.1, Strain.2), as.character)) %>% 
    eccAverages(., ecc_data)
  
  results <- results %>% left_join(., dist_avgs)
  
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

eccAverages <- function(pw_dists, ecc_data) {
  # Given the cluster assignments at TP1, TP2, and the pairwise distances for all strains, 
  # what are the pairwise distances of strains found *within* each cluster
  clusters <- ecc_data %>% select(Strain, grep("Size|ECC", colnames(.), value = TRUE, invert = TRUE)) %>% 
    rename(TP1 = grep("TP1", colnames(.), value = TRUE), 
           TP2 = grep("TP2", colnames(.), value = TRUE))
  
  df <- left_join(pw_dists, clusters, by = c("Strain.1" = "Strain")) %>% 
    rename(first_in_TP1 = TP1, first_in_TP2 = TP2) %>% 
    left_join(., clusters, by = c("Strain.2" = "Strain")) %>% 
    rename(second_in_TP1 = TP1, second_in_TP2 = TP2)
  
  tplist <- c("TP1", "TP2")
  tpcases <- lapply(tplist, function(m) {
    df %>% select(grep(setdiff(tplist, m), colnames(.), invert = TRUE)) %>% withinClusterDists(., m)  
  }) %>% set_names(tplist)
  
  left_join(clusters, tpcases$TP1) %>% left_join(., tpcases$TP2) %>% select(-TP1, -TP2) %>% return()
}

withinClusterDists <- function(df, tp) {
  tpx_case <- df %>% 
    set_colnames(gsub(tp, "TPX", colnames(.))) %>% 
    filter(!is.na(first_in_TPX) & !is.na(second_in_TPX)) %>% 
    filter(first_in_TPX == second_in_TPX) %>% 
    select(-second_in_TPX) %>% 
    rename(TPX = first_in_TPX) %>% 
    group_by(TPX)
  
  tpx_case %>% summarise(TPX_avg_geo_dists = mean(Geog.Dist), TPX_avg_temp_dists = mean(Temp.Dist)) %>% 
    set_colnames(gsub("TPX", tp, colnames(.))) %>% return()
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


