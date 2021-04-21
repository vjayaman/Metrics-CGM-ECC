# Vasena's edits - for speed improvement on Windows, where we have less control over cores 

epi_cohesion_sep <- function(g_cuts, epi_matrix, cpus){
  
  genome_names <- g_cuts %>% select(genome) %>% pull()

  epi_matrix <- epi.matrix
  epi_matrix[lower.tri(epi_matrix)] <- NA
  # epi_melt <- melt(as.matrix(1-epi_matrix)) # create similarity values from epi distance matrix
  epi_melt <- melt(as.matrix(1-epi_matrix))
  
  # x1 <- epi_melt %>% as_tibble() %>% filter(Var1 == "Anhui/SZ005/2020") %>% arrange(Var2)
  # x2 <- epi_melt %>% as_tibble() %>% filter(Var2 == "Anhui/SZ005/2020") %>% arrange(Var1)
  
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

  full_set %>% #arrange(cluster) %>%
    mutate(ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))) %>%
    select(-cut, -members, -s1) %>% 
    set_colnames(c(names(g_cuts)[2], paste0(names(g_cuts)[2], "_Size"), 
                   paste0(names(g_cuts)[2], "_ECC")))
}

# # If we remove the duplicate pairwise data (i.e. only store one value for each pair)
# 
# epi_cohes_cal_na <- function(g_cuts, epi_matrix, cpus){
#   
#   genome_names <- g_cuts %>% select(genome) %>% pull()
#   
#   # create similarity values from epi distance matrix
#   epi_melt <- melt(as.matrix(1-epi_matrix)) %>% as_tibble()
#   
#   epi_melt_joined <- 
#     expand_grid(genome_names, genome_names, .name_repair = function(x) {c("Var1", "Var2")}) %>% 
#     left_join(epi_melt, by = c("Var1", "Var2"))
#   
#   # this is the added part (since no duplicated pairwise dists now 
#   # - for the same strains but in reverse, that is)
#   epi_melt_joined <- epi_melt_joined %>% filter(!is.na(value))
#   
#   # epi_melt_joined %<>% select(Var2, Var1, value) %>% rename(Var1 = Var2, Var2 = Var1) %>% 
#   #   bind_rows(epi_melt_joined, .)
#   
#   calculate_s1 <- function(k) {
#     epi_melt_joined %>% 
#       filter(Var1 %in% k, Var2 %in% k) %>%
#       select(value) %>% 
#       pull() %>% 
#       sum()
#   }
#   
#   # print("Starting Calculation")
#   
#   cut_cluster_members <-
#     g_cuts %>%
#     pivot_longer(-genome, names_to = "cut", values_to = "cluster") %>%  
#     group_by(cut, cluster) %>% 
#     summarise(members = list(cur_data()$genome), .groups = "drop")
#   
#   # print("Part 2")  
#   cut_cluster_members <- cut_cluster_members %>% 
#     add_column(cluster_size = map_int(cut_cluster_members$members, length))
#   
#   # non-singletons:
#   others <- cut_cluster_members %>% filter(cluster_size > 1) %>% 
#     mutate(s1 = map_dbl(members, calculate_s1))
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
#   full_set %>% arrange(cluster) %>%
#     mutate(ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))) %>%
#     select(-cut, -members, -s1) %>% 
#     set_colnames(c(names(g_cuts)[2], paste0(names(g_cuts)[2], "_Size"), 
#                    paste0(names(g_cuts)[2], "_ECC")))
# }
