n1 <- wk_list[i]
n1strains <- metadata[num_wk <= n1]$Strain
n1clusters <- unique(pull(f2[isolate %in% n1strains], heights[1]))
outputDetails(paste0("Week ", n1, " has ", length(n1strains), " strains and ", 
                     length(n1clusters), " clusters"), newcat = TRUE)
n2 <- wk_list[i+1]
n2strains <- metadata[num_wk <= n2]$Strain
n2clusters <- unique(pull(f2[isolate %in% n2strains], heights[1]))
outputDetails(paste0("Up to week ", n2, " has ", length(c(n1strains, n2strains)), 
                     " strains and ", length(n2clusters), " clusters"), newcat = TRUE)


# f2 %>% filter(!!as.symbol(heights[1]) %in% pull(f2[isolate %in% n1strains], heights[1])) %>% pu

# tpx1 <- clusters[[i]][[1]]
# 
# b1 <- metadata[num_wk <= wk_list[i]] %>% pull(Strain)
# b2 <- f2 %>% filter(!!as.symbol(heights) %in% pull(f2[isolate %in% b1], 2))
# assert("All TP1 cluster strains <= week i are included", 
#        all(identical(pull(tpx1, "isolate"), pull(b2, "isolate")), 
#            identical(pull(tpx1, heights[1]), pull(b2, heights[1]))))
# 
# tpx2 <- clusters[[i+1]][[1]] %>% bind_rows(., tpx1) %>% unique() %>% arrange(isolate)
# assert("All TP1 strains found in TP2 as well", all(tpx1$isolate %in% tpx2$isolate))
# 
# b3 <- metadata[num_wk <= wk_list[i+1]] %>% pull(Strain)
# b4 <- f2 %>% filter(!!as.symbol(heights) %in% pull(f2[isolate %in% b3], 2))
# assert("All TP2 cluster strains <= week i+1 are included", 
#        all(identical(pull(tpx2, "isolate"), pull(b4, "isolate")), 
#            identical(pull(tpx2, heights[1]), pull(b4, heights[1]))))


# n1 <- wk_list[i]
# tpx1 <- clusters[[i]]$sofar %>% rename(isolate = Strain)
# n2 <- wk_list[i+1]
# tpx2 <- clusters[[i+1]]$sofar %>% rename(isolate = Strain)
