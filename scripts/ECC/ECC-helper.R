# Dillon's version

epi_cohesion_calc <- function(g_cuts, epi_matrix, cpus){
  
  genome_names <- g_cuts %>% select(genome) %>% pull()
  
  epi_melt <- melt(as.matrix(1-epi_matrix)) # create similarity values from epi distance matrix
  
  epi_melt_joined <- 
    expand_grid(genome_names, genome_names, .name_repair = function(x) {c("Var1", "Var2")}) %>% 
    left_join(epi_melt)
  
  calculate_s1 <- function(k) {
    
    epi_melt_joined %>% 
      filter(Var1 %in% k, Var2 %in% k) %>%
      select(value) %>% 
      pull() %>% 
      sum()
      
  }

  print("Starting Calculation")
  
  cut_cluster_members <-
    g_cuts %>%
    pivot_longer(-genome, names_to = "cut", values_to = "cluster") %>%  
    group_by(cut, cluster) %>% 
    summarise(
      members = list(cur_data()$genome)
    )
  
  cut_cluster_members %>% 
    mutate(
      s1 = map_dbl(members, calculate_s1),
      cluster_size = map_int(members, length),
      ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1)),
      W_ECC = ECC * cluster_size
      # PW_Values = Have to figure out a tidy way of doing this
    ) %>% 
    select(-members, -s1)

    # calculate_row_old <- function(k) {
    #   
    #   xselect <- expand.grid(y[[k]], y[[k]])
    #   #returns the epi similarities for each of the pairwise combinations in the cluster from the epi_sim melted list
    #   xselect1 <- merge(x = xselect, y = epi_melt, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
    #   
    #   n1 <- length(unique(xselect1[,1]))
    #   s1 <- sum(xselect1[,3])
    #   IEC <- (s1-n1) / (n1*(n1-1))
    #   W_IEC <- IEC * n1
    #   
    #   row <- tibble(
    #     Cluster_Name = names(y[k]),
    #     Cluster_Size = n1,
    #     ECC = IEC,
    #     W_ECC = W_IEC,
    #     Members = paste(unique(xselect[,2]), collapse = ","),
    #     PW_Values = paste(xselect1[,3], collapse = ",")
    #   )
    #   row
    # }
  
}
