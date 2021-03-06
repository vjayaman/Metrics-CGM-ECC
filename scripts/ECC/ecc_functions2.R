library(Rcpp)

cppFunction('int clusterSizes(DataFrame x, int h) {
  NumericVector cluster = x[0];
  NumericVector y = x[2];
  int nrow = x.nrow();
  int vals = 0;
  
  for (int i = 0; i < nrow; i++) {
    if (cluster[i] == h) {
      vals = vals + y[i];
    }
  }
  return vals;
}')

# cppFunction('List matchedDrs(DataFrame ccm, DataFrame dr_as) {
#   int nrow = ccm.nrow();
#   
#   NumericVector x = ccm[0]; //cluster column of ccm
#   NumericVector col1 = dr_as[0]; //cluster column of dr_assignments
#   NumericVector col2 = dr_as[1]; //dr column of dr_assignments
#   LogicalVector inds;
#   
#   for (int i = 0; i < nrow; i++) {
#     inds = (col1 == x[i]); //rows in dr_assignments with cluster x[i]
#   }
#   return testv;
# }')
# 
# inds <- matchedDrs(ccm, dr_assignments)

# NumericVector k = y[i]; //list of cluster members of cluster at row i

cppFunction('LogicalVector drsInCluster(DataFrame ccm, int i, DataFrame dr_as) {
  NumericVector x = ccm[0]; //cluster column of ccm
  
  NumericVector col1 = dr_as[0]; //cluster column of dr_assignments
  NumericVector col2 = dr_as[1]; //dr column of dr_assignments
  LogicalVector inds = (col1 == x[i]); //rows in dr_assignments with cluster x[i]
  
  return inds;
}')

cppFunction('double sumEpiVals(DataFrame epi_melt) {
  NumericVector value = epi_melt[2];
  NumericVector n1 = epi_melt[3];
  NumericVector n2 = epi_melt[4];
  
  NumericVector value2 = value * n1 * n2;
  return sum(value2);
}')

epiCohesion <- function(g_cuts, epi_melt) {
  epi_melt <- epi_melt[!is.na(value)]
  dr_assignments <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))
  
  uniclusters <- unique(pull(g_cuts, 1))
  sizes <- lapply(uniclusters, function(h) clusterSizes(g_cuts, h)) %>% 
    unlist() %>% tibble(cluster = uniclusters, cluster_size = .)
  
  cut_cluster_members <-
    g_cuts %>% select(-n) %>%
    pivot_longer(-dr, names_to = "cut", values_to = "cluster") %>%
    group_by(cut, cluster) %>%
    summarise(members = list(cur_data()$dr), .groups = "drop") %>%
    left_join(., sizes, by = "cluster") %>% as.data.table()
  
  ccm <- cut_cluster_members[,c(2,4)]
  y <- cut_cluster_members$members
  ivals <- cut_cluster_members$cluster
  
  tmp1 <- lapply(0:(nrow(cut_cluster_members)-1), function(j) {
    matches <- dr_assignments[drsInCluster(ccm, j, dr_assignments)] %>% select(-cluster)  
    i <- ivals[j + 1]
    k <- cut_cluster_members[cluster == i, members] %>% unlist()
    epi_filtered <- epi_melt[dr_1 %in% matches$dr & dr_2 %in% matches$dr] %>% 
      inner_join(., matches, by = c("dr_1" = "dr")) %>% rename(n1 = n) %>% 
      inner_join(., matches, by = c("dr_2" = "dr")) %>% rename(n2 = n)
    sumEpiVals(epi_filtered)
  }) %>% unlist() %>% bind_cols(cut_cluster_members, s1 = .)

  th <- names(g_cuts)[1]
  tmp1 %>%
    mutate(ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))) %>%
    select(-cut, -members, -s1) %>%
    set_colnames(c(th, paste0(th, "_Size"), paste0(th, "_ECC")))
  return(tmp2)
}
# New
# 0, 3
# 0, 23
# 4, 0
# 8, 33
# 4, 32
# Original
# 0, 3
# 0, 25
# 3, 57
# 8, 33
# 3, 59
# So, will need to work on cutting down the distances one - takes too long
# As well as epiCohesion() and a couple of the CGM ones