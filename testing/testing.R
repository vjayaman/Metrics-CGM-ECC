libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table")
y <- lapply(libs, require, character.only = TRUE)

# FORMATTING RESULTS FROM ORIGINAL

a1 <- readRDS("results/original_9000_collected_ECC.tsv")

# source, temporal, geographic coefficients
# source = 0
# temporal = 1
# geographical = 0

a1[[1]] %<>% select(-Strain) %>% unique()
a1[[2]] %<>% select(-Strain) %>% unique()

selectCols <- function(df, keep, leave) {
  grep(keep, colnames(df), value = TRUE) %>% 
    grep(leave, ., invert = TRUE, value = TRUE) %>% return()
}

oldtp1temp <- selectCols(a1[[1]], "TP1", "avg")
oldtp2temp <- selectCols(a1[[1]], "TP2", "avg")

oldtp1geo <- selectCols(a1[[2]], "TP1", "avg")
oldtp2geo <- selectCols(a1[[2]], "TP2", "avg")

c1 <- list(
  list(unique(a1[[1]][,oldtp1temp]), unique(a1[[1]][,oldtp2temp])), 
  list(unique(a1[[2]][,oldtp1geo]), unique(a1[[2]][,oldtp2geo]))
)
  
c1[[1]][[1]] %<>% filter(!is.na(!!as.symbol(oldtp1temp[1])))
c1[[2]][[1]] %<>% filter(!is.na(!!as.symbol(oldtp1geo[1])))


# FORMATTING RESULTS FROM NEW

b1 <- readRDS("results/new_9000_collected_ECC.tsv")

# temporal, geographic coefficients
# temporal = 1
# geographical = 0

newtp1temp <- selectCols(b1[[1]][[1]], "TP1", "avg")
newtp2temp <- selectCols(b1[[1]][[2]], "TP2", "avg")

newtp1geo <- selectCols(b1[[2]][[1]], "TP1", "avg")
newtp2geo <- selectCols(b1[[2]][[2]], "TP2", "avg")

b1[[1]][[1]] %<>% as_tibble() %>% select(all_of(newtp1temp)) %>% 
  set_colnames(colnames(c1[[1]][[1]]))
b1[[1]][[2]] %<>% as_tibble() %>% select(all_of(newtp2temp)) %>% 
  set_colnames(colnames(c1[[1]][[2]]))

b1[[2]][[1]] %<>% as_tibble() %>% select(all_of(newtp1geo)) %>% 
  set_colnames(colnames(c1[[2]][[1]]))
b1[[2]][[2]] %<>% as_tibble() %>% select(all_of(newtp2geo)) %>% 
  set_colnames(colnames(c1[[2]][[2]]))

newset <- b1[[2]][[2]]
oldset <- c1[[2]][[2]]

df1 <- newset %>% set_colnames(c("Cluster", "Newsize", "NewECC"))
df2 <- oldset %>% set_colnames(c("Cluster", "Oldsize", "OldECC"))

df3 <- full_join(df1, df2, by = "Cluster")
df3[is.na(df3)] <- 1
df3$Size_dif <- abs(df3$Oldsize - df3$Newsize)
df3$ECC_dif <- abs(df3$OldECC - df3$NewECC)

notmatching <- df3 %>% filter(ECC_dif > 1e-14)

# - need to see which set has the actually correct sizes, and why
# assert("Cluster sizes match", all(abs(pull(newmethod, 2) - pull(oldmethod, 2)) < 1e-14))


# CHECK ORIGINAL METHOD CLUSTER SIZES
source("testing/test_original_clusters.R")
source("testing/test_new_clusters.R")
assert("Strain data is same for both methods", identical(sd1, sd2))
assert("TP1 clusters are the same for both methods", identical(td1[[1]], td2[[1]]))
assert("TP2 clusters are the same for both methods", identical(td1[[2]], td2[[2]]))

# arbitrary:
typing_data <- td1
strains <- sd1

tp1sizes <- typing_data[[1]] %>% rownames_to_column() %>% 
  set_colnames(c("Strain", "Cluster")) %>% as_tibble() %>% 
  group_by(Cluster) %>% summarise(n = n())

tp2sizes <- typing_data[[2]] %>% rownames_to_column() %>% 
  set_colnames(c("Strain", "Cluster")) %>% as_tibble() %>% 
  group_by(Cluster) %>% summarise(n = n())

if (grepl("TP1", colnames(newset)[1])) {
  actual_sizes <- tp1sizes %>% 
    filter(Cluster %in% notmatching$Cluster) %>% 
    left_join(notmatching, ., by = "Cluster") %>% rename(Actualsize = n)  
}else {
  actual_sizes <- tp2sizes %>% 
    filter(Cluster %in% notmatching$Cluster) %>% 
    left_join(notmatching, ., by = "Cluster") %>% rename(Actualsize = n)
}

res_new <- all(actual_sizes$Newsize - actual_sizes$Actualsize == 0)
res_old <- all(actual_sizes$Oldsize - actual_sizes$Actualsize == 0)

if (isTRUE(res_old)) {
  print("Original method had correct clusters")
}else if (isTRUE(res_new)) {
  print("New method had correct clusters")
}else if (!isTRUE(res_new) & !isTRUE(res_old)) {
  print("Neither had correct clusters")
}

# inds2 <- which(abs(b1[[1]][[2]]$TP2_T0_ECC_0.1.0 - c1[[1]][[2]]$TP2_T0_ECC_0.1.0) > 1e-14)

# for (cx in colnames(a1)) {
#   x <- all(identical(pull(a1, cx), pull(b1, cx)))
#   if (!x) {
#     print(cx)
#   }
# }  
