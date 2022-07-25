testDists <- function(strain_data, temps = NULL, geos = NULL) {
  a2 <- strain_data %>% select(Strain, Longitude, Latitude) %>% as.data.frame() %>% 
    column_to_rownames("Strain") %>% earth.dist(dist = TRUE) %>% as.matrix() %>% 
    set_rownames(strain_data$Strain) %>% set_colnames(strain_data$Strain)
  if (is.null(geos)) {
    a3 <- a2 %>% transformData2(., "geo", min(a2), max(a2)) %>% 
      formatData(., c("Strain1", "Strain2", "Geog.Dist"))
  }else {
    a3 <- a2 %>% transformData2(., "geo", geos$minval, geos$maxval) %>% 
      formatData(., c("Strain1", "Strain2", "Geog.Dist"))
  }
  
  b1 <- strain_data %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% select(Strain, Date)
  b2 <- b1 %>% pull(Date) %>% 
    dist(diag = FALSE, upper = FALSE, method = "euclidean") %>% 
    as.matrix(nrow = nrow(b1), ncol = nrow(b1)) %>% 
    set_rownames(b1$Strain) %>% set_colnames(b1$Strain)
  if (is.null(temps)) {
    b3 <- b2 %>% transformData2(., "temp", min(b2), max(b2)) %>% 
      formatData(., c("Strain1", "Strain2", "Temp.Dist"))
  }else {
    b3 <- b2 %>% transformData2(., "temp", temps$minval, temps$maxval) %>% 
      formatData(., c("Strain1", "Strain2", "Temp.Dist"))
  }
  
  transformed_dists <- merge.data.table(b3, a3)
  return(transformed_dists)
}

testEpiMelt <- function(transformed_dists, tau, gamma) {
  epi_table <- transformed_dists %>% 
    mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
    select(Strain1, Strain2, Total.Dist) %>% as.data.table()
  
  epi_matrix <- dcast(epi_table, formula = Strain1 ~ Strain2, value.var = "Total.Dist")
  epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)]) 
  rownames(epi_matrix) <- colnames(epi_matrix)
  
  epi_melt <- as.matrix(1-epi_matrix) %>% 
    as.data.table(., keep.rownames = TRUE) %>% 
    melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>% 
    set_colnames(c("Var1", "Var2", "value")) %>% as.data.table()
  return(epi_melt)
}

testECCs <- function(g_cuts, epi_melt) {
  genome_names <- g_cuts %>% select(genome) %>% pull()
  
  epi_melt_joined <- 
    expand_grid(genome_names, genome_names, .name_repair = function(x) {c("Var1", "Var2")}) %>% 
    left_join(epi_melt, by = intersect(colnames(epi_melt), colnames(.)))
  
  calculate_s1 <- function(k) {
    epi_melt_joined %>% filter(Var1 %in% k, Var2 %in% k) %>%
      select(value) %>% pull() %>% sum()
  }
  
  cut_cluster_members <- g_cuts %>%
    pivot_longer(-genome, names_to = "cut", values_to = "cluster") %>%  
    group_by(cut, cluster) %>% 
    summarise(members = list(cur_data()$genome), .groups = "keep")
  
  actual_eccs <- cut_cluster_members %>% 
    mutate(s1 = map_dbl(members, calculate_s1), 
           cluster_size = map_int(members, length), 
           ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))) %>% 
    select(-members, -s1) %>% ungroup() %>% as.data.table()
  return(actual_eccs)
}

checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

# function for reading raw strains and time point clusters
readData <- function(path, check_enc = TRUE) {
  if (check_enc) {
    enc <- checkEncoding(file.path(path))
  }else {
    enc <- ""
  }
  
  file.path(path) %>% 
    read.table(sep="\t", header=TRUE, fileEncoding=enc, fill=TRUE, quote="") %>% 
    as_tibble() %>% return()
}
