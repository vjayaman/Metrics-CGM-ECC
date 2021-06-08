
prepEpiTable <- function(transformed_dists, tau, gamma) {
  cat(paste0("\n   Collecting ECC values for temporal = ", tau, ", geo = ", gamma))
  
  outputMessages(paste0("      Preparing table of distances: sqrt(", tau, "*(temp^2) + ", gamma, "*(geo^2))"))
  transformed_dists %>% 
    mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
    select(dr1, dr2, Total.Dist) %>% as.data.table() %>% 
    return()
}

prepEpiMatrix <- function(epitable) {
  outputMessages("      Generating matrix of similarity values from epi distance matrix ...")
  epi_matrix <- dcast(epitable, formula = dr1 ~ dr2, value.var = "Total.Dist")
  epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)]) 
  rownames(epi_matrix) <- colnames(epi_matrix)
  
  # create similarity values from epi distance matrix:
  as.matrix(1 - epi_matrix) %>% return()
}

### Incorporating the allele data with the epidemiological data 
epiCollection <- function(sim_matrix, typing_data, tau, gamma, dr_matches, avgdistvals) {
  
  epi_melt <- sim_matrix %>% as.data.table(., keep.rownames = TRUE) %>% 
    melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>% 
    set_colnames(c("dr_1", "dr_2", "value")) %>% as.data.table()
  
  # ---------------------------------------------------------------------------------------------
  # epimatrix <- matrix(nrow = nrow(strain_data), ncol = nrow(strain_data)) %>% 
  #   set_colnames(strain_data$Strain) %>% set_rownames(strain_data$Strain)
  # x1 <- left_join(epi_melt, dr_matches, by = c("dr_1" = "dr"))
  # x2 <- x1 %>% rename(Strain1 = Strain) %>% select(dr_2, value, Strain1)
  # x3 <- left_join(x2, dr_matches, by = c("dr_2" = "dr"))
  
  # ---------------------------------------------------------------------------------------------
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  eccs <- lapply(1:length(typing_data), function(i) {
    outputMessages(paste0("      Calculating ECCs for timepoint ", i, ", which has ", 
                          nrow(typing_data[[i]]), " strains to consider ..."))
    
    dr_td1 <- typing_data[[i]] %>% rownames_to_column("Strain") %>% as_tibble() %>%
      left_join(., dr_matches, by = "Strain") %>%
      mutate(across(dr, as.character)) %>% select(-Strain)
    
    # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
    cx <- colnames(dr_td1)[1]
    tallied_reps <- dr_td1 %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
    g_cuts <- left_join(dr_td1, tallied_reps, by = intersect(colnames(tallied_reps), colnames(dr_td1))) %>%
      unique() %>% mutate(across(dr, as.character))
    
    td_i <- epi_cohesion_new(g_cuts, epi_melt) %>% 
      set_colnames(c(paste0("TP", i, "_", colnames(.))))
    colnames(td_i) %<>% gsub("ECC", paste0("ECC.0.", tau, ".", gamma), x = .)
    
    a2 <- avgdistvals[[i]]$temp
    b2 <- avgdistvals[[i]]$geo
    
    left_join(td_i, a2, by = intersect(colnames(td_i), colnames(a2))) %>% 
      left_join(., b2, by = intersect(colnames(.), colnames(b2))) %>% return()
  })
  
  return(eccs)
}

drCounts <- function(td_i, dr_matches) {
  dr_td1 <- td_i %>% rownames_to_column("Strain") %>% as_tibble() %>%
    left_join(., dr_matches, by = "Strain") %>%
    mutate(across(dr, as.character)) %>% select(-Strain)
  
  # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
  cx <- colnames(dr_td1)[1]
  tallied_reps <- dr_td1 %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
  tallied_reps %>% left_join(dr_td1, ., by = byNames(., dr_td1)) %>%
    unique() %>% mutate(across(dr, as.character)) %>% return()
}

# note that we sum the values for both directions e.g. (192, 346) and (346, 192)
# and then divide by the total number of pairs (counting both directions)
# this gives the same result as if we had used only one direction
#   - if we did this, we would divide by 2 in the numerator and the denominator --> would cancel
avgDists <- function(g_cuts, dm, cname, newname) {
  x <- gsub("\\.", "_", cname) %>% tolower() %>% paste0(newname, "_avg_", .)
  all_clusters <- g_cuts %>% set_colnames(c("Threshold", "dr", "n"))
  clusters <- unique(all_clusters$Threshold)
  
  lapply(clusters, function(cl) {
    f2 <- all_clusters %>% filter(Threshold == cl)
    dists <- dm[f2$dr, f2$dr]
    if (is.null(dim(dists))) {
      return(0)
    }else {
      dists <- dists * f2$n[col(dists)]
      dists <- dists * f2$n[row(dists)]
      return(sum(dists) / (sum(f2$n)^2))
    }
  }) %>% unlist() %>% tibble(clusters, .) %>% set_colnames(c(newname, x))
}

# FORMATTING ---------------------------------------------------------------------------------
checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

byNames <- function(df1, df2) {intersect(colnames(df1), colnames(df2))}

mergeECCs <- function(eccs, tpx, typing_data) {
  tbl1 <- as.data.table(typing_data)
  sapply(eccs, "[", tpx) %>% Reduce(function(...) merge(...), .) %>% as.data.table() %>% 
    merge.data.table(tbl1, ., by = intersect(colnames(tbl1), colnames(.))) %>% return()
}

hyps <- function(x) {paste0(rep("-", x), collapse = "")}

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

# EPI-HELPER-MODULAR -----------------------------------------------------------------------------
outputMessages <- function(msgs) {cat(paste0("\n",msgs))}

nonRedundantDists <- function(assignments) {
  
  outputMessages("   Generating all possible date pair distances ...")
  dm_temp <- assignments %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
  transformed_temp <- dm_temp %>% pairwiseDists(., "temp", c("dr1", "dr2", "Temp.Dist"))
  
  outputMessages("   Generating all possible lat-long pair distances ...")
  dm_geo <- assignments %>% select(dr, Latitude, Longitude) %>% 
    distMatrix(., "geo", c("Longitude", "Latitude"))
  transformed_geo <- dm_geo %>% pairwiseDists(., "geo", c("dr1", "dr2", "Geog.Dist"))
  
  transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
  dists <- list("temp" = dm_temp, "geo" = dm_geo, "tr_dists" = transformed_dists)
  return(dists)
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
    expand_grid(dr_names, dr_names, .name_repair = function(x) {c("dr_1", "dr_2")}) %>%
    left_join(., epi_melt, by = c("dr_1", "dr_2")) %>% as.data.table()
  
  sizes <- lapply(unique(dr_assignments$cluster), function(h) {
    dr_assignments %>% filter(cluster == h) %>% pull(n) %>% sum() %>% tibble(cluster = h, cluster_size = .)
  }) %>% bind_rows()
  
  cut_cluster_members <-
    g_cuts %>% select(-n) %>% 
    pivot_longer(-dr, names_to = "cut", values_to = "cluster") %>%
    group_by(cut, cluster) %>%
    summarise(members = list(cur_data()$dr), .groups = "drop") %>%
    left_join(., sizes, by = "cluster") %>% as.data.table()
  
  calculate_s1 <- function(i) {
    k <- cut_cluster_members[cluster == i, members] %>% unlist()
    matches <- dr_assignments %>% filter(cluster == i)
    
    epi_melt_joined[dr_1 %in% k][dr_2 %in% k] %>% 
      left_join(., matches, by = c("dr_1" = "dr")) %>% rename(n1 = n) %>% select(-cluster) %>% 
      left_join(., matches, by = c("dr_2" = "dr")) %>% rename(n2 = n) %>% select(-cluster) %>% 
      mutate(value2 = value * n1 * n2) %>%
      select(value2) %>% pull() %>% sum()
  }
  
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

pairwiseDists <- function(dm, type, newnames) {
  transformData(dm, type) %>% formatData(., newnames) %>% return()
}