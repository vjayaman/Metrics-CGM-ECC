
# EPI-HELPER-MODULAR -----------------------------------------------------------------------------
pairwiseDists <- function(dm, type, newnames) {
  # dm <- assignments %>% distMatrix(., type, cnames)
  transformed <- transformData(dm, type) %>% formatData(., newnames)
    # as.data.frame() %>% rownames_to_column("dr1") %>% as.data.table() %>% 
    # melt.data.table(., id.vars = "dr1", variable.name = "dr2", value.name = newnames[3]) %>%
    # as_tibble() %>% 
    # mutate(dr2 = as.character(dr2)) %>% 
    # as.data.table() %>% set_colnames(newnames)
  return(transformed)
}

distMatrix <- function(input_data, dtype, cnames) {
  if (dtype == "temp") {
    dm <- input_data %>% select(all_of(cnames)) %>% pull() %>% 
      dist(diag = FALSE, upper = FALSE, method = "euclidean")
    dm %>% as.matrix(nrow = nrow(input_data), ncol = nrow(input_data)) %>% return()
    
  }else if (dtype == "geo") {
    # consider using geosphere::distm() for this
    dm <- input_data %>% select(all_of(cnames)) %>% as.data.frame() %>% earth.dist(dist = TRUE)
    dm %>% as.matrix() %>% return()
  }
}

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

formatData <- function(dm, newnames) {
  dm %>% 
    as.data.frame() %>% rownames_to_column("dr1") %>% as.data.table() %>% 
    melt.data.table(., id.vars = "dr1", variable.name = "dr2", value.name = newnames[3]) %>%
    as_tibble() %>% 
    mutate(dr2 = as.character(dr2)) %>% 
    as.data.table() %>% set_colnames(newnames) %>% return()
}

epi_cohesion_new <- function(g_cuts, epi_melt) {
  
  dr_names <- g_cuts %>% select(dr) %>% pull() %>% unique()
  dr_assignments <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))
  
  epi_melt_joined <-
    expand_grid(dr_names, dr_names, .name_repair = function(x) {c("Var1", "Var2")}) %>%
    left_join(., epi_melt, by = c("Var1", "Var2")) %>% as.data.table()
  
  sizes <- lapply(unique(dr_assignments$cluster), function(h) {
    dr_assignments %>% filter(cluster == h) %>% pull(n) %>% sum() %>% tibble(cluster = h, cluster_size = .)
  }) %>% bind_rows()
  
  
  calculate_s1 <- function(i) {
    
    k <- cut_cluster_members %>% filter(cluster == i) %>% pull(members) %>% unlist()
  
    matches <- dr_assignments %>% filter(cluster == i)
    
    epi_melt_joined %>% filter(Var1 %in% k & Var2 %in% k) %>% 
      left_join(., matches, by = c("Var1" = "dr")) %>% rename(n1 = n) %>% select(-cluster) %>% 
      left_join(., matches, by = c("Var2" = "dr")) %>% rename(n2 = n) %>% select(-cluster) %>% 
      mutate(value2 = value * n1 * n2) %>%
      select(value2) %>% pull() %>% sum()
  }
  
  # print("Starting Calculation")
  cut_cluster_members <-
    g_cuts %>% select(-n) %>% 
    pivot_longer(-dr, names_to = "cut", values_to = "cluster") %>%
    group_by(cut, cluster) %>%
    summarise(members = list(cur_data()$dr), .groups = "drop") %>%
    left_join(., sizes, by = "cluster")
  
  th <- names(g_cuts)[1]
  
  cut_cluster_members %>%
    mutate(
      s1 = map_dbl(cluster, calculate_s1),
      ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))
    ) %>%
    ungroup() %>% 
    select(-cut, -members, -s1) %>%
    set_colnames(c(th, paste0(th, "_Size"), paste0(th, "_ECC")))
}

# EPI-HELPER -------------------------------------------------------------------------------------
#Function to return table of all the epi-similarities and final strain similarity #
# datafile <- strain_data; source_matrix <- source_pw; source_coeff <- sigma;
# temp_coeff <- tau; geog_coeff <- gamma; geog_temp <- matched
EpiTable <- function(datafile, source_matrix = NULL, source_coeff, temp_coeff, geog_coeff, geog_temp){
  
  x <- source_coeff 
  y <- temp_coeff
  z <- geog_coeff
  
  
  #### Create the pairwise table for lookups ####
  d <- expand.grid(1:nrow(datafile), 1:nrow(datafile))
  
  d1 <- d[,1]
  d2 <- d[,2]
  
  # split into two steps, since it seems to seems to reduce memory usage  
  if (x == 0) {
    strain_sims <- tibble(
      "Strain.1" = datafile[d1, "Strain"]  %>% pull (),
      "Strain.2" = datafile[d2, "Strain"]  %>% pull (),
      "Date.1"   = datafile[d1, "Date"]  %>% pull (),
      "Date.2"   = datafile[d2, "Date"]  %>% pull (),
      "Location.1" = datafile[d1, "Location"]  %>% pull (),
      "Location.2" = datafile[d2, "Location"] %>% pull()
    )
    
    str.matrix <-
      strain_sims %>% 
      left_join(geog_temp, by = c("Strain.1", "Strain.2"))
    # This is necessary (otherwise any NA values in the Source.Dist column make Total.Dist and Epi.Sym NA as well. 
    # There is a better way to handle this; will work on that.)
    str.matrix <-
      str.matrix %>% 
      mutate(
        Total.Dist = sqrt( ((Temp.Dist^2)*y) + ((Geog.Dist^2)*z) ),
        Epi.Sym = 1 - Total.Dist
      )
    
  }else {
    strain_sims <- tibble(
      "Strain.1" = datafile[d1, "Strain"]  %>% pull (),
      "Strain.2" = datafile[d2, "Strain"]  %>% pull (),
      "Source.1" = datafile[d1, "Source"]  %>% pull (),
      "Source.2" = datafile[d2, "Source"]  %>% pull (),
      "Date.1"   = datafile[d1, "Date"]  %>% pull (),
      "Date.2"   = datafile[d2, "Date"]  %>% pull (),
      "Location.1" = datafile[d1, "Location"]  %>% pull (),
      "Location.2" = datafile[d2, "Location"] %>% pull()
    )
    
    str.matrix <-
      strain_sims %>% 
      left_join(source_matrix, by = c("Source.1", "Source.2")) %>% 
      left_join(geog_temp, by = c("Strain.1", "Strain.2"))
    
    str.matrix <-
      str.matrix %>% 
      mutate(
        Total.Dist = sqrt( (((Source.Dist^2)*x) + ((Temp.Dist^2)*y) + ((Geog.Dist^2)*z)) ),
        Epi.Sym = 1 - Total.Dist
      ) 
  }
  
  str.matrix
}  

# datafile <- strain_data; source_matrix <- source_pw; source_coeff <- sigma; temp_coeff <- tau
# geog_coeff <- gamma; geog_temp <- matched
EpiTable2 <- function(datafile, source_matrix, source_coeff, temp_coeff, geog_coeff, geog_temp){
  #### Read data into memory from previous outputs ####
  
  x <- source_coeff 
  y <- temp_coeff
  z <- geog_coeff
  
  #### Create the pairwise table for lookups ####
  d <- expand.grid(1:nrow(datafile), 1:nrow(datafile))
  d1 <- d[,1]
  d2 <- d[,2]
  
  # split into two steps, since it seems to seems to reduce memory usage  
  if (x == 0) {
    strain_sims <- tibble(
      "Strain.1" = datafile[d1, "Strain"]  %>% pull (),
      "Strain.2" = datafile[d2, "Strain"]  %>% pull (),
      "Date.1"   = datafile[d1, "Date"]  %>% pull (),
      "Date.2"   = datafile[d2, "Date"]  %>% pull (),
      "Location.1" = datafile[d1, "Location"]  %>% pull (),
      "Location.2" = datafile[d2, "Location"] %>% pull()
    )
    
    # This is necessary (otherwise any NA values in the Source.Dist column make Total.Dist and Epi.Sym NA as well. 
    # There is a better way to handle this; will work on that.)
    str.matrix <- strain_sims %>% 
      left_join(geog_temp, by = c("Strain.1", "Strain.2")) %>% 
      mutate(
        Total.Dist = sqrt( ((Temp.Dist^2)*y) + ((Geog.Dist^2)*z) ),
        Epi.Sym = 1 - Total.Dist
      )    
  }else {
    strain_sims <- tibble(
      "Strain.1" = datafile[d1, "Strain"]  %>% pull (),
      "Strain.2" = datafile[d2, "Strain"]  %>% pull (),
      "Source.1" = datafile[d1, "Source"]  %>% pull (),
      "Source.2" = datafile[d2, "Source"]  %>% pull (),
      "Date.1"   = datafile[d1, "Date"]  %>% pull (),
      "Date.2"   = datafile[d2, "Date"]  %>% pull (),
      "Location.1" = datafile[d1, "Location"]  %>% pull (),
      "Location.2" = datafile[d2, "Location"] %>% pull()
    )
    
    str.matrix <-
      strain_sims %>% 
      left_join(source_matrix, by = c("Source.1", "Source.2")) %>% 
      left_join(geog_temp, by = c("Strain.1", "Strain.2")) %>% 
      mutate(
        Total.Dist = sqrt( (((Source.Dist^2)*x) + ((Temp.Dist^2)*y) + ((Geog.Dist^2)*z)) ),
        Epi.Sym = 1 - Total.Dist
      ) 
  }
  
  str.matrix
}

#Function to return matrix of just strains and final similarity scores for building graphics
# epi.matrix <- str.matrix
EpiMatrix <- function(epi.matrix){
  
  epi.matrix <- epi.matrix %>% # I think, TODO: check
    select(Strain.1, Strain.2, Total.Dist)
  
  epi.cast <- dcast(epi.matrix, formula= Strain.1 ~ Strain.2, value.var = "Total.Dist")
  epi.cast <- as.matrix(epi.cast[,2:ncol(epi.cast)]) 
  rownames(epi.cast) <- colnames(epi.cast)
  #   epi.sym <- 1 - epi.cast
  
  epi.cast
}


