# source("scripts/part1.R")

# i <- 1
#   k <- as.character(dfx$k[i])
#   tpkstrains <- metadata[get(interval) <= k]$Strain
#   key_cls <- parts$drs[Strain %in% tpkstrains] %>% select(-Strain, -dr) %>% pull() %>% unique()
#   y <- lapply(parts$results, function(x) any(key_cls %in% pull(x, 1))) %>% unlist()
#   
#   k_drs <- m$dr_matches %>% filter(Strain %in% tpkstrains) %>% pull(dr)
#   
#   c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "-") %>% unlist() %>% as.numeric()
#   tau <- c1[2]
#   gamma <- c1[3]
#   
#   fnames <- names(y[y])
#   
#   outputMessages(paste0("Collecting ECCs for cluster groups at TP", 
#                         dfx$k[i], ", ECC coefficients ", as.character(dfx$x[i]), "\n"))
#   
#   allres <- bind_rows(results)  
#   cluster_x <- df[df[[cx]] %in% pull(allres, cx),-"Strain"]
#   
#   transformed_dists <- lapply(fnames, function(f) {
#     readRDS(paste0("intermediate_data/TP", k, "/dists/group", f, ".Rds")) %>% collectTransforms(., extremes)
#   }) %>% bind_rows() %>% unique()
#     
#   selected_tp <- m$strain_data %>% filter(Strain %in% tpkstrains)
#     
#   # eccs <- epiCollectionByCluster()
#   strain_data <- selected_tp
#   # tau <- tau; gamma <- gamma; transformed_dists <- transformed_dists
#   tpx <- k
#   cluster_y <- cluster_x[dr %in% k_drs]
#   
#   # epiCollectionByCluster <- function() {
#     epi_table <- transformed_dists %>% 
#       mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
#       select(dr1, dr2, Total.Dist) %>% as.data.table()
#     rm(transformed_dists); invisible(gc())
#     epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
#     epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)]) 
#     rownames(epi_matrix) <- colnames(epi_matrix)
#     epi_melt <- as.matrix(1-epi_matrix) %>% as.data.table(., keep.rownames = TRUE) %>% 
#       melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>% 
#       set_colnames(c("dr_1", "dr_2", "value")) %>% as.data.table()
#     rm(epi_table); rm(epi_matrix); invisible(gc())
#     cx <- colnames(cluster_y)[1]
#     tallied_reps <- cluster_y %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
#     cnames <- intersect(colnames(tallied_reps), colnames(cluster_y))
#     g_cuts <- left_join(cluster_y, tallied_reps, by = cnames) %>%
#       unique() %>% mutate(across(dr, as.character))
#     # write.table(g_cuts, "epicoh/g_cuts.txt", sep = "\t")
#     # write.table(epi_melt, "epicoh/epi_melt.txt", sep = "\t")
#     
#     # td_i <- epiCohesion(g_cuts, epi_melt) %>% 
#     #   set_colnames(c(paste0("TP", tpx, "_", colnames(.))))
#     # colnames(td_i) %<>% gsub("ECC", paste0("ECC.0.", tau, ".", gamma), x = .)

source("scripts/part1.R")

epi_melt <- read.table("epicoh/epi_melt.txt", sep = "\t", header = TRUE) %>% as.data.table() 
#%>% mutate(across(c(dr_1, dr_2), as.character))
g_cuts <- read.table("epicoh/g_cuts.txt", sep = "\t", header = TRUE) %>% as.data.table() 
#%>% mutate(across(dr, as.character))

    
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

# NumericVector k = y[i]; //list of cluster members of cluster at row i
cppFunction('LogicalVector drsInCluster(DataFrame ccm, int i, DataFrame dr_as) {
  NumericVector x = ccm[0]; //cluster column of ccm
  
  NumericVector col1 = dr_as[0]; //cluster column of dr_assignments
  NumericVector col2 = dr_as[1]; //dr column of dr_assignments
  LogicalVector inds = (col1 == x[i]); //rows in dr_assignments with cluster x[i]
  
  return inds;
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

cppFunction('double sumEpiVals(DataFrame epi_melt) {
  NumericVector value = epi_melt[2];
  NumericVector n1 = epi_melt[3];
  NumericVector n2 = epi_melt[4];
  
  NumericVector value2 = value * n1 * n2;
  return sum(value2);
}')




    # epiCohesion <- function(g_cuts, epi_melt) {
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
        epi_filtered <- epi_melt %>% 
          inner_join(., matches, by = c("dr_1" = "dr")) %>% rename(n1 = n) %>% 
          inner_join(., matches, by = c("dr_2" = "dr")) %>% rename(n2 = n)
        sumEpiVals(epi_filtered)
      }) %>% unlist() %>% bind_cols(cut_cluster_members, s1 = .)
      
      
      th <- names(g_cuts)[1]
      tmp2 <- tmp1 %>%
        mutate(ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))) %>%
        select(-cut, -members, -s1) %>%
        set_colnames(c(th, paste0(th, "_Size"), paste0(th, "_ECC")))
    
      
      
      
# RESULTS
      # epiCohesion <- function(g_cuts, epi_melt) {
      dr_assignments <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))

      uniclusters <- unique(pull(g_cuts, 1))
#     # BETTER
      sizes <- lapply(uniclusters, function(h) clusterSizes(g_cuts, h)) %>%
        unlist() %>% tibble(cluster = uniclusters, cluster_size = .)
#     # WORSE
      # sizes <- lapply(uniclusters, function(h) {
      #   dr_assignments %>% filter(cluster == h) %>%
      #     pull(n) %>% sum() %>% tibble(cluster = h, cluster_size = .)
      # }) %>% bind_rows()

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
      tmp3 <- cut_cluster_members %>%
        mutate(
          s1 = map_dbl(cluster, calculate_s1),
          ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))
        ) %>%
        ungroup() %>%
        select(-cut, -members, -s1) %>%
        set_colnames(c(th, paste0(th, "_Size"), paste0(th, "_ECC")))
      # }

      # CHECKING RESULTS
      x1 <- tmp2$T0_ECC - tmp3$T0_ECC
      any(x1[!is.na(x1)] > 1e-12)