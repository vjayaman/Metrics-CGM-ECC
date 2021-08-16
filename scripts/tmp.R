
libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

# Current working directory should be Metrics-CGM-ECC/
files <- c("scripts/ECC/classes_ecc.R", "scripts/ECC/ecc_functions.R", "scripts/ECC/dist_functions.R")
invisible(sapply(files, source))

# # list(strain_data = selected_tp, tau = 1, gamma = 0, transformed_dists = transformed_dists,
# #      tpx = k, cluster_y = cluster_x[dr %in% k_drs]) %>% saveRDS(., "epicoh/epicollect.Rds")
# 
# a <- readRDS("epicoh/epicollect.Rds")
# 
# strain_data <- a$strain_data
# tau <- a$tau
# gamma <- a$gamma
# transformed_dists <- a$transformed_dists
# tpx <- a$tpx
# cluster_y <- a$cluster_y
# 
# epi_table <- transformed_dists %>% mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
#   select(dr1, dr2, Total.Dist) %>% as.data.table()
# 
# epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
# epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)]) 
# rownames(epi_matrix) <- colnames(epi_matrix)
# 
# epi_melt <- as.matrix(1-epi_matrix) %>% as.data.table(., keep.rownames = TRUE) %>% 
#   melt(., id.vars = "rn", variable.factor = FALSE, value.factor = FALSE) %>% 
#   set_colnames(c("dr_1", "dr_2", "value")) %>% as.data.table()
# 
# cx <- colnames(cluster_y)[1]
# tallied_reps <- cluster_y %>% group_by(!!as.symbol(cx)) %>% count(dr) %>% ungroup()
# cnames <- intersect(colnames(tallied_reps), colnames(cluster_y))
# g_cuts <- left_join(cluster_y, tallied_reps, by = cnames) %>%
#   unique() %>% mutate(across(dr, as.character))
# 
# # list(g_cuts = g_cuts, epi_melt = epi_melt) %>% saveRDS(., "epicoh/epicohesion.Rds")

# td_i <- epiCohesion(g_cuts, epi_melt) %>% 
#   set_colnames(c(paste0("TP", tpx, "_", colnames(.))))
# colnames(td_i) %<>% gsub("ECC", paste0("ECC.0.", tau, ".", gamma), x = .)

b <- readRDS("epicoh/epicohesion.Rds")
g_cuts <- b$g_cuts
# write.table(g_cuts, "epicoh/g_cuts.txt", sep = "\t")
epi_melt <- b$epi_melt
# write.table(epi_melt, "epicoh/epi_melt.txt", sep = "\t")

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
  # inner_join(., dr_assignments, by = "cluster")

calculate_s1 <- function(i) {
  k <- cut_cluster_members[cluster == i, members] %>% unlist()
  matches <- dr_assignments[cluster == i]
  
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

