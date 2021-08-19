# R version 1
epiCohesionV1 <- function(g_cuts, epi_melt) {
  dr_assignments <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))
  
  sizes <- lapply(unique(dr_assignments$cluster), function(h) {
    dr_assignments %>% filter(cluster == h) %>%
      pull(n) %>% sum() %>% tibble(cluster = h, cluster_size = .)
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
    
    epi_melt[dr_1 %in% k][dr_2 %in% k] %>% 
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

# C++ version 1
epiCohesionV2 <- function(g_cuts, epi_melt) {
  g_cuts <- g_cuts %>% as.data.table() %>% mutate(across(dr, as.integer))
  th <- names(g_cuts)[1]
  dr_assignments <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))
  uniclusters <- unique(pull(g_cuts, 1))
  
  sizes <- lapply(uniclusters, function(h) clusterSizes(g_cuts, h)) %>% 
    unlist() %>% data.table(cluster = uniclusters, cluster_size = .)
  
  ux <- unique(as.matrix(sizes[cluster_size > 1]$cluster))
  epi_melt <- epi_melt %>% add_column(flag1 = 0, flag2 = 0) %>% 
    mutate(across(c(dr_1, dr_2), as.integer))
  
  df <- calculateEpiV2(ux, as.matrix(dr_assignments), as.matrix(epi_melt))
  
  if (length(sizes[cluster_size == 1]$cluster) > 0) {
    epi_vals <- data.table(cluster = ux[,1], s1 = df) %>% 
      bind_rows(., data.table(cluster = sizes[cluster_size == 1]$cluster, s1 = NA))
  }else {
    epi_vals <- data.table(cluster = ux[,1], s1 = df)
  }
  
  epi_vals %>% 
    left_join(., sizes, by = "cluster") %>% 
    select(cluster, cluster_size, s1) %>%
    mutate(ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))) %>% select(-s1) %>% 
    set_colnames(c(th, paste0(th, "_Size"), paste0(th, "_ECC"))) %>%
    arrange(!!as.symbol(th)) %>% return()
}

# R version 2
epiCohesionV3 <- function(g_cuts, epi_melt) {
  dr_assignments <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))
  
  sizes <- lapply(unique(dr_assignments$cluster), function(h) {
    dr_assignments %>% filter(cluster == h) %>%
      pull(n) %>% sum() %>% data.table(cluster = h, cluster_size = .)
  }) %>% bind_rows()
  
  th <- names(g_cuts)[1]
  g_cuts <- g_cuts %>% set_colnames(c("hx", "dr", "n"))
  # clusters with more than one pair of drs
  multiclusters <- g_cuts %>% group_by(hx) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(hx)
  
  # clusters with only one unique pair of drs (either singletons or non unique clusters)
  uniclusters <- g_cuts %>% group_by(hx) %>% summarise(n = n()) %>% filter(n == 1) %>% pull(hx)
  
  eccres <- lapply(multiclusters, function(cluster_x) {
    cx <- g_cuts[hx == cluster_x]
    epitbl <- epi_melt[dr_1 %in% cx$dr & dr_2 %in% cx$dr] %>% unique()
    
    tmp <- dcast(epitbl, formula = dr_1 ~ dr_2, value.var = "value")
    tmp <- as.matrix(tmp[,2:ncol(tmp)])
    rownames(tmp) <- colnames(tmp)
    
    cx <- cx[match(colnames(tmp), cx$dr)]
    rownames(tmp) <- cx$n
    tmp2 <- sapply(1:ncol(tmp),function(x) tmp[,x]*as.integer(rownames(tmp)[x]) )
    colnames(tmp2) <- cx$n
    tmp3 <- sapply(1:nrow(tmp2),function(x) tmp2[x,]*as.integer(colnames(tmp2)[x]) )
    
    cluster_size <- sizes %>% filter(cluster == cluster_x) %>% pull(cluster_size)
    
    (sum(tmp3) - cluster_size) / (cluster_size * (cluster_size - 1))
  }) %>% unlist() %>% data.table(cluster = multiclusters, ECC = .)
  
  singletons <- data.table(sizes[cluster %in% uniclusters & cluster_size == 1]) %>% 
    add_column(ECC = NA, .after = 2)
  
  nonunique <- data.table(sizes[cluster %in% uniclusters & cluster_size != 1]) %>% 
    add_column(ECC = 1, .after = 2)
  
  eccres %>% inner_join(., sizes, by = "cluster") %>% 
    bind_rows(., singletons) %>% bind_rows(., nonunique) %>% select(cluster, cluster_size, ECC) %>% 
    set_colnames(c(th, paste0(th, "_Size"), paste0(th, "_ECC"))) %>%
    arrange(!!as.symbol(th)) %>% return()
}

# R version 3
epiCohesionV4 <- function(g_cuts, epi_melt) {
  dr_assignments <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))
  
  sizes <- lapply(unique(dr_assignments$cluster), function(h) {
    dr_assignments %>% filter(cluster == h) %>%
      pull(n) %>% sum() %>% data.table(cluster = h, cluster_size = .)
  }) %>% bind_rows()
  
  th <- names(g_cuts)[1]
  g_cuts <- g_cuts %>% set_colnames(c("hx", "dr", "n"))
  # clusters with more than one pair of drs
  multiclusters <- g_cuts %>% group_by(hx) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(hx)
  
  # clusters with only one unique pair of drs (either singletons or non unique clusters)
  uniclusters <- g_cuts %>% group_by(hx) %>% summarise(n = n()) %>% filter(n == 1) %>% pull(hx)
  
  eccres <- lapply(multiclusters, function(cluster_x) {
    cx <- g_cuts[hx == cluster_x]

    epitbl <- epi_melt[rownames(epi_melt) %in% cx$dr, colnames(epi_melt) %in% cx$dr] %>% as.matrix()
    
    cx <- cx[match(colnames(epitbl), cx$dr)]
    rownames(epitbl) <- cx$n
    tmp2 <- sapply(1:ncol(epitbl),function(x) epitbl[,x]*as.integer(rownames(epitbl)[x]) )
    colnames(tmp2) <- cx$n
    tmp3 <- sapply(1:nrow(tmp2),function(x) tmp2[x,]*as.integer(colnames(tmp2)[x]) )
    
    cluster_size <- sizes %>% filter(cluster == cluster_x) %>% pull(cluster_size)
    
    (sum(tmp3) - cluster_size) / (cluster_size * (cluster_size - 1))
  }) %>% unlist() %>% data.table(cluster = multiclusters, ECC = .)
  
  singletons <- data.table(sizes[cluster %in% uniclusters & cluster_size == 1]) %>% 
    add_column(ECC = NA, .after = 2)
  
  nonunique <- data.table(sizes[cluster %in% uniclusters & cluster_size != 1]) %>% 
    add_column(ECC = 1, .after = 2)
  
  eccres %>% inner_join(., sizes, by = "cluster") %>% 
    bind_rows(., singletons) %>% bind_rows(., nonunique) %>% select(cluster, cluster_size, ECC) %>% 
    set_colnames(c(th, paste0(th, "_Size"), paste0(th, "_ECC"))) %>%
    arrange(!!as.symbol(th)) %>% return()
}
