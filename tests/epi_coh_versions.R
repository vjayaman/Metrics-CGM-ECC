libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
          "reshape2","fossil","tidyr","purrr", "data.table", "Rcpp", "RcppArmadillo")
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))
sourceCpp("scripts/epicohversions.cpp")
source("scripts/ECC/ecc_functions.R")


g_cuts <- readRDS("intermediate_data/g_cuts.Rds")
epi_melt <- readRDS("intermediate_data/epi_melt.Rds") %>% unique()

t1 <- Sys.time()
df1 <- epiCohesionV1(g_cuts, epi_melt)
t2 <- Sys.time()
df2 <- epiCohesionV2(g_cuts, epi_melt)
t3 <- Sys.time()
df3 <- epiCohesionV3(g_cuts, epi_melt)
t4 <- Sys.time()

paste0("epiCohV1 took ", t2 - t1)
paste0("epiCohV2 took ", t3 - t2)
paste0("epiCohV3 took ", t4 - t3)

a <- df1 %>% set_colnames(c("hx", "size", "v1_ecc"))
b <- df2 %>% set_colnames(c("hx", "size", "v2_ecc"))
d <- df3 %>% set_colnames(c("hx", "size", "v3_ecc"))
alleccs <- inner_join(a, b, by = c("hx", "size")) %>% 
  inner_join(., d, by = c("hx", "size")) %>% filter(size > 1)
assert("ECC results for V1 and V2 are within tolerance 1e-10", all(abs(alleccs$v1_ecc - alleccs$v2_ecc) < 1e-10))
assert("ECC results for V1 and V3 are within tolerance 1e-10", all(abs(alleccs$v1_ecc - alleccs$v3_ecc) < 1e-10))
assert("ECC results for V2 and V3 are within tolerance 1e-10", all(abs(alleccs$v2_ecc - alleccs$v3_ecc) < 1e-10))

# # saveRDS(g_cuts, "intermediate_data/g_cuts.Rds")
# # saveRDS(epi_melt, "intermediate_data/epi_melt.Rds")



d <- inner_join(y, eccres, by = "cluster")
assert("ECCs of new R method match old R method", all(abs(d$oldECC - d$ECC) < 1e-12))

paste0("old R method took ", t2 - t1)
paste0("new R method took ", t4 - t3)