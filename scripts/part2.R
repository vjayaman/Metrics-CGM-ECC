source("scripts/part1.R")
library(Rcpp)
sourceCpp("scripts/cppfunctions.cpp")
i <- 2

  k <- as.character(dfx$k[i])
  tpkstrains <- metadata[get(interval) <= k]$Strain
  key_cls <- parts$drs[Strain %in% tpkstrains] %>% select(-Strain, -dr) %>% pull() %>% unique()
  y <- lapply(parts$results, function(x) any(key_cls %in% pull(x, 1))) %>% unlist()

  k_drs <- m$dr_matches %>% filter(Strain %in% tpkstrains) %>% pull(dr)

  c1 <- as.character(dfx$x[i]) %>% strsplit(., split = "-") %>% unlist() %>% as.numeric()
  tau <- c1[2]
  gamma <- c1[3]

  fnames <- names(y[y])

  outputMessages(paste0("Collecting ECCs for cluster groups at TP",
                        dfx$k[i], ", ECC coefficients ", as.character(dfx$x[i]), "\n"))

  allres <- bind_rows(results)
  cluster_x <- df[df[[cx]] %in% pull(allres, cx),-"Strain"]

  transformed_dists <- lapply(fnames, function(f) {
    readRDS(paste0("intermediate_data/TP", k, "/dists/group", f, ".Rds")) %>% collectTransforms(., extremes)
  }) %>% bind_rows() %>% unique()

  selected_tp <- m$strain_data %>% filter(Strain %in% tpkstrains)

  # eccs <- epiCollectionByCluster()
  strain_data <- selected_tp
  # tau <- tau; gamma <- gamma; transformed_dists <- transformed_dists
  tpx <- k
  cluster_y <- cluster_x[dr %in% k_drs]

  # epiCollectionByCluster <- function() {
    epi_table <- transformed_dists %>%
      mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>%
      select(dr1, dr2, Total.Dist) %>% as.data.table()
    rm(transformed_dists); invisible(gc())
    epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
    epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)])
    rownames(epi_matrix) <- colnames(epi_matrix)
    epi_melt <- as.matrix(1-epi_matrix) %>% as.data.table(., keep.rownames = TRUE) %>%
      melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>%
      set_colnames(c("dr_1", "dr_2", "value")) %>% as.data.table()
    rm(epi_table); rm(epi_matrix); invisible(gc())
    cx <- colnames(cluster_y)[1]
    tallied_reps <- cluster_y %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
    cnames <- intersect(colnames(tallied_reps), colnames(cluster_y))
    g_cuts <- left_join(cluster_y, tallied_reps, by = cnames) %>%
      unique() %>% mutate(across(dr, as.character))
    # write.table(g_cuts, "epicoh/g_cuts.txt", sep = "\t")
    # write.table(epi_melt, "epicoh/epi_melt.txt", sep = "\t")

    # td_i <- epiCohesion(g_cuts, epi_melt) %>%
    #   set_colnames(c(paste0("TP", tpx, "_", colnames(.))))
    # colnames(td_i) %<>% gsub("ECC", paste0("ECC.0.", tau, ".", gamma), x = .)





# source("scripts/part1.R")
# epi_melt <- read.table("epicoh/epi_melt.txt", sep = "\t", header = TRUE) %>% as.data.table() 
epi_melt <- epi_melt %>% as.data.table() %>% add_column(flag1 = 0, flag2 = 0) %>% 
  mutate(across(c(dr_1, dr_2), as.integer))
epi_melt <- epi_melt[!is.na(value)]
#%>% mutate(across(c(dr_1, dr_2), as.character))
# g_cuts <- read.table("epicoh/g_cuts.txt", sep = "\t", header = TRUE) %>% as.data.table() 
g_cuts <- g_cuts %>% as.data.table() %>% mutate(across(dr, as.integer))



th <- names(g_cuts)[1]
dr_assignments <- g_cuts %>% set_colnames(c("cluster", "dr", "n"))

uniclusters <- unique(pull(g_cuts, 1))
sizes <- lapply(uniclusters, function(h) clusterSizes(g_cuts, h)) %>% 
  unlist() %>% tibble(cluster = uniclusters, cluster_size = .)

ux <- unique(as.matrix(dr_assignments$cluster))
df <- calculateEpi(ux, as.matrix(dr_assignments), as.matrix(epi_melt))
  
uniclusters <- unique(pull(g_cuts, 1))
sizes <- lapply(uniclusters, function(h) clusterSizes(g_cuts, h)) %>% 
  unlist() %>% tibble(cluster = uniclusters, cluster_size = .)
  
epi_vals <- tibble(cluster = unique(dr_assignments$cluster), s1 = df) %>% 
  left_join(., sizes, by = "cluster") %>% 
  select(cluster, cluster_size, s1) %>%
  mutate(ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))) %>% select(-s1) %>% 
  set_colnames(c(th, paste0(th, "_Size"), paste0(th, "_ECC"))) %>%
  arrange(!!as.symbol(th))

      
      
      

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
        set_colnames(c(th, paste0(th, "_Size"), paste0(th, "_ECC"))) %>% 
        arrange(!!as.symbol(th))

      # CHECKING RESULTS
      x1 <- epi_vals$T0_ECC - tmp3$T0_ECC
      assert("ECC results of both methods are the same", all(x1[!is.na(x1)] < 1e-12))