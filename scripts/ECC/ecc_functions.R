
# EPI-HELPER-MODULAR -----------------------------------------------------------------------------
# assignments <- assignments %>% select(dr, Date)
# type <- "temp"; cnames <- "Date"; newnames <- c("dr1", "dr2", "Temp.Dist")
pairwiseDists <- function(assignments, type, cnames, newnames) {
  dm <- assignments %>% distMatrix(., type, cnames)
  # formattedDM <- dm %>% 
  #   as.data.frame() %>% rownames_to_column("dr1") %>% as.data.table() %>% 
  #   melt.data.table(., id.vars = "dr1", variable.name = "dr2", 
  #                   value.name = newnames[3], variable.factor = FALSE) %>% as_tibble()
  transformed <- transformData(dm, type) %>% 
    as.data.frame() %>% rownames_to_column("dr1") %>% as.data.table() %>% 
    melt.data.table(., id.vars = "dr1", variable.name = "dr2", value.name = newnames[3]) %>%
    as_tibble() %>% 
    mutate(dr2 = as.character(dr2)) %>% 
    as.data.table() %>% set_colnames(newnames)
  return(transformed)
}

distMatrix <- function(input_data, dtype, cnames) {
  if (dtype == "temp") {
    dm <- input_data %>% select(all_of(cnames)) %>% pull() %>% 
      dist(diag = FALSE, upper = FALSE, method = "euclidean")
    # bigmemory::as.big.matrix() %>% return()
    # dist(diag = TRUE, upper = TRUE, method = "euclidean")
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

# # ECC-SEP-SINGLETONS -----------------------------------------------------------------------------
# # Vasena's edits - for speed improvement on Windows, where we have less control over cores 
# 
epi_cohesion_sep <- function(g_cuts, epi_matrix, cpus){

  genome_names <- g_cuts %>% select(genome) %>% pull()

  epi_matrix <- epi.matrix
  epi_melt <- melt(as.matrix(1-epi_matrix)) # create similarity values from epi distance matrix

  epi_melt_joined <-
    expand_grid(genome_names, genome_names, .name_repair = function(x) {c("Var1", "Var2")}) %>%
    left_join(epi_melt, by = c("Var1", "Var2")) %>% filter(!is.na(value))

  calculate_s1 <- function(k) {
    epi_melt_joined %>%
      filter(Var1 %in% k & Var2 %in% k) %>%
      select(value) %>%
      pull() %>%
      sum()
  }

  # print("Starting Calculation")

  cut_cluster_members <-
    g_cuts %>%
    pivot_longer(-genome, names_to = "cut", values_to = "cluster") %>%
    group_by(cut, cluster) %>%
    summarise(members = list(cur_data()$genome), .groups = "drop")

  # print("Part 2")
  cut_cluster_members <- cut_cluster_members %>%
    add_column(cluster_size = map_int(cut_cluster_members$members, length))

  # non-singletons:
  # k <- cut_cluster_members %>% filter(cluster_size > 1) %>% slice(1) %>% pull(members) %>% unlist()
  others <- cut_cluster_members %>% filter(cluster_size > 1) %>%
    mutate(s1 = map_dbl(members, calculate_s1))

  # print("Part 3")
  singletons <- cut_cluster_members %>% filter(cluster_size == 1)
  if (nrow(singletons) > 0) {
    singletons %<>% mutate(across(members, unlist)) %>% mutate(mem2 = members) %>%
      left_join(., epi_melt_joined, by = c("members" = "Var1", "mem2" = "Var2")) %>%
      rename(s1 = value) %>% select(-mem2) %>%
      mutate(across(members, as.list))
    full_set <- bind_rows(singletons, others)
  }else {
    full_set <- others
  }

  full_set %>% arrange(cluster) %>%
    mutate(ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))) %>%
    select(-cut, -members, -s1) %>%
    set_colnames(c(names(g_cuts)[2], paste0(names(g_cuts)[2], "_Size"),
                   paste0(names(g_cuts)[2], "_ECC")))
}

# ECC-SEP-SINGLETONS -----------------------------------------------------------------------------
# Vasena's edits - for speed improvement on Windows, where we have less control over cores 

# workingEachCluster <- function(g_cuts, epi_melt, strain_drs, cs) {
# 
#   genome_names <- g_cuts %>% select(dr) %>% pull()
# 
#   epi_melt_joined <-
#     expand_grid(genome_names, genome_names, .name_repair = function(x) {c("Var1", "Var2")}) %>%
#     left_join(., epi_melt, by = c("Var1", "Var2")) %>% filter(!is.na(value)) %>% as.data.table()
# 
#   epi_melt_joined <- left_join(epi_melt_joined, strain_drs, by = c("Var1" = "dr")) %>%
#     select(-Var1) %>% rename(Var1 = Strain) %>%
#     left_join(., strain_drs, by = c("Var2" = "dr")) %>%
#     select(-Var2) %>% rename(Var2 = Strain) %>%
#     select(Var1, Var2, value) %>% arrange(Var1, Var2)
# 
#   full_cuts <- g_cuts %>% left_join(., strain_drs) %>% select(T0, Strain)
# 
#   calculate_s1 <- function(k) {
#     epi_melt_joined %>%
#       filter(Var1 %in% k & Var2 %in% k) %>%
#       select(value) %>% pull() %>% sum()
#   }
# 
#   # print("Starting Calculation")
#   cut_cluster_members <-
#     full_cuts %>%
#     pivot_longer(-Strain, names_to = "cut", values_to = "cluster") %>%
#     group_by(cut, cluster) %>%
#     summarise(members = list(cur_data()$Strain), .groups = "drop") %>%
#     mutate(cluster_size = map_int(members, length))
# 
#   sums <- cut_cluster_members %>%
#     mutate(
#       s1 = map_dbl(members, calculate_s1),
#       ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))
#     ) %>%
#     ungroup()
# 
#   sums %>%
#     select(-cut, -members, -s1) %>%
#     set_colnames(c(names(g_cuts)[2], paste0(names(g_cuts)[2], "_Size"),
#                    paste0(names(g_cuts)[2], "_ECC")))
# }



eachCluster <- function(g_cuts, epi_melt) {
  
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
    left_join(., sizes)
  
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

# eachCluster <- function(g_cuts, epi_melt) {
#   
#   genome_names <- g_cuts %>% select(dr) %>% pull()
# 
#   epi_melt_joined <-
#     expand_grid(genome_names, genome_names, .name_repair = function(x) {c("Var1", "Var2")}) %>%
#     left_join(., epi_melt, by = c("Var1", "Var2")) %>% filter(!is.na(value)) %>% as.data.table()
#   
#   drcounts <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))
#   sizes <- lapply(unique(drcounts$cluster), function(i) {
#     drcounts %>% filter(cluster == i) %>% pull(n) %>% sum() %>% tibble(cluster = i, cluster_size = .) %>% return()
#   }) %>% bind_rows()
#   
#   # allcounts <- left_join(epi_melt, drcounts, by = c("Var1"="dr")) %>% 
#   #   mutate(value2 = value * n) %>% select(-n) %>%
#   #   left_join(., drcounts, by = c("Var2" = "dr")) %>%
#   #   mutate(value3 = value2 * n) %>% select(-n) %>%
#   #   select(Var1, Var2, value3) %>%
#   #   rename(value = value3)
# 
#   calculate_s1 <- function(i) {
#     k <- cut_cluster_members %>% filter(cluster == i) %>% pull(members) %>% unlist()
#     current <- drcounts %>% filter(cluster == i)
#     # epi_melt_joined %>%
#     #   filter(Var1 %in% k & Var2 %in% k) %>%
#     epi_melt_joined %>%
#       filter(Var1 %in% k & Var2 %in% k) %>% 
#       left_join(., current, by = c("Var1" = "dr")) %>% rename(n1 = n) %>% select(-cluster) %>% 
#       left_join(., current, by = c("Var2" = "dr")) %>% rename(n2 = n) %>% select(-cluster) %>% 
#       # mutate(n2 = ifelse(Var1 == Var2, 1, n2)) %>%
#       mutate(value2 = value * n1 * n2 /2) %>% 
#       select(value2) %>% 
#       # select(value) %>% 
#       pull() %>% 
#       sum()
#   }
#   
#   cut_cluster_members <-
#     g_cuts %>% select(-n) %>%
#     pivot_longer(-dr, names_to = "cut", values_to = "cluster") %>%
#     group_by(cut, cluster) %>%
#     summarise(members = list(cur_data()$dr), .groups = "drop") %>% 
#     left_join(., sizes)
#   
#   # non-singletons:
#   # k <- cut_cluster_members %>% filter(cluster_size > 1) %>% slice(1) %>% pull(members) %>% unlist()
#   others <- cut_cluster_members %>% filter(cluster_size > 1) %>% 
#     mutate(s1 = map_dbl(cluster, calculate_s1))
#     # mutate(s1 = map_dbl(c(members, cluster), calculate_s1))
# 
#   # print("Part 3")
#   singletons <- cut_cluster_members %>% filter(cluster_size == 1)
#   if (nrow(singletons) > 0) {
#     singletons %<>% mutate(across(members, unlist)) %>% mutate(mem2 = members) %>%
#       left_join(., epi_melt_joined, by = c("members" = "Var1", "mem2" = "Var2")) %>%
#       rename(s1 = value) %>% select(-mem2) %>%
#       mutate(across(members, as.list))
#     full_set <- bind_rows(singletons, others)
#   }else {
#     full_set <- others
#   }
# 
#   cval <- full_set$cut %>% unique()
# 
#   full_set %>% arrange(cluster) %>%
#     mutate(ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))) %>%
#     select(-cut, -members, -s1) %>%
#     set_colnames(c(cval, paste0(cval, "_Size"), paste0(cval, "_ECC")))
# }

epi_cohesion_new <- function(g_cuts, epi_melt){
  eachCluster(g_cuts, epi_melt) %>% return()
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

# COLLECTING-ECCS ------------------------------------------------------------------------------
checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}


