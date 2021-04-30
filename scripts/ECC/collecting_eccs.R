
checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

### Incorporating the allele data with the epidemiological data 
oneCombo <- function(strain_data, source_file, sigma, tau, gamma, cpus, typing_data) {
  cat(paste0("\nCollecting ECC values for source = ", sigma, ", temporal = ", tau, ", geo = ", gamma))
  ### Generating the EpiMatrix
  # strain_data <- read_tsv(strains) %>% 
  #   mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
  #          Location = paste(Country, Province, City, sep = "_"))
  
  source_pw <- read_tsv(source_file) %>% 
    mutate(Source.Dist = 1- value) %>% select(-value)
  
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
  
  # typing_datum <- typing_data[[1]]

  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  eccs <- lapply(typing_data, function(typing_datum) {
    g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
    # Method with singletons merged separately, for speed/efficiency:
    epi_cohesion_sep(g_cuts, epi.matrix, cpus)
    # Original method: epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
  })
  
  return(eccs)

  # # names(eccs) == c("TP1_T0", "TP2_T0")
  # newnames <- sapply(strsplit(names(eccs), "_"), `[`, 1)
  # 
  # ecc_data <- lapply(names(eccs), function(x) {
  #   y <- eccs[[x]]
  #   td[[x]] %>% rownames_to_column("Strain") %>% as_tibble() %>% 
  #     left_join(., y, by = colnames(y)[1]) %>% 
  #     set_colnames(c("Strain", paste0(unlist(strsplit(x, "_"))[1], "_", colnames(y)))) %>% 
  #     set_colnames(gsub("ECC", paste0("ECC_", sigma, ".", tau, ".", gamma), colnames(.)))
  #   
  # }) %>% set_names(newnames) %>% 
  #   Reduce(function(...) right_join(..., by = "Strain"), .)
  # 
  # # # Adding basic average columns (Date, Longitude, Latitude)
  # # avg_raw <- strain_data %>% select(Strain, Date, Longitude, Latitude)
  # # cnames <- colnames(ecc_data) %>% grep("Size|ECC", ., value = TRUE, invert = TRUE)
  # # df <- ecc_data %>% select(all_of(cnames)) %>% set_colnames(c("Strain", "TP1", "TP2"))
  # # 
  # # avg_tp1 <- df %>% basicAverages(., as.name("TP1"), avg_raw) %>% arrange(Strain) %>% select(-TP1)
  # # avg_tp2 <- df %>% basicAverages(., as.name("TP2"), avg_raw) %>% arrange(Strain) %>% select(-TP2)
  # # 
  # # 
  # # # Adding distance averages
  # # clusters <- ecc_data %>% select(Strain, grep("Size|ECC", colnames(.), value = TRUE, invert = TRUE)) %>% 
  # #   rename(TP1 = grep("TP1", colnames(.), value = TRUE),
  # #          TP2 = grep("TP2", colnames(.), value = TRUE))
  # # 
  # # # matching to see which pairs are found in the same cluster
  # # pw_dists <- left_join(geog_pw$raw, temp_pw$raw, by = c("Strain.1", "Strain.2")) %>% 
  # #   left_join(., clusters, by = c("Strain.1" = "Strain")) %>% 
  # #   rename(first_in_TP1 = TP1, first_in_TP2 = TP2) %>% 
  # #   left_join(., clusters, by = c("Strain.2" = "Strain")) %>% 
  # #   rename(second_in_TP1 = TP1, second_in_TP2 = TP2)
  # # 
  # # dist_avgs <- lapply(c("TP1", "TP2"), function(i) {
  # #   x1 <- grep(i, colnames(pw_dists), value = TRUE)
  # #   a1 <- as.name(x1[1])
  # #   a2 <- as.name(x1[2])  
  # #   pw_dists %>% filter(all_of({{a1}}) == {{a2}}) %>% 
  # #     rename(TPX = all_of(a1)) %>% group_by(TPX) %>% 
  # #     summarise(TPX_avg_geo_dist = mean(Geog.Dist), TPX_avg_temp_dist = mean(Temp.Dist)) %>% 
  # #     set_colnames(gsub("TPX", i, colnames(.))) %>% 
  # #     left_join(clusters, .) %>% select(-TP1, -TP2)
  # # }) %>% Reduce(function(...) left_join(..., by = "Strain"), .)
  # #   
  # # results <- left_join(ecc_data, avg_tp1, by = "Strain") %>% 
  # #   left_join(., avg_tp2, by = "Strain") %>% 
  # #   left_join(., dist_avgs, by = "Strain")
  # # 
  # # return(results)
  # return(ecc_data)
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


