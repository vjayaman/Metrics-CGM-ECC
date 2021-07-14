
# Functions for test-formatting_cgm.R and test-cgm_classes.R -----------------------------
sampleSizes <- function(minval, maxval, numvals) {
  sample(seq(minval, maxval), numvals, replace = TRUE)
}

syntheticTypeSet <- function(type = c(1,2,3,4), nums = rep(25, 4)) {
  fset <- tibble()
  if (1 %in% type) {
    tp1 <- sampleSizes(3, 200, nums[1])
    fset <- fset %>% bind_rows(tibble(tp1_cl_size = tp1, tp2_cl_size = tp1))
  }
  
  if (2 %in% type) {
    tp1 <- sampleSizes(3, 200, nums[2])
    fset <- tibble(tp1_cl_size = tp1, tp2_cl_size = tp1 + sampleSizes(1, 50, nums[2])) %>% 
      bind_rows(fset, .)
  }
  
  if (3 %in% type) {
    fset <- tibble(tp1_cl_size = sampleSizes(1, 2, nums[3]), 
                   tp2_cl_size = sampleSizes(3, 200, nums[3])) %>% 
      bind_rows(fset, .)
  }
  
  if (4 %in% type) {
    tp1 <- sampleSizes(1, 2, nums[4])
    tp2 <- sampleSizes(1, 2, nums[4])
    fset <- tibble(tp1_cl_size = tp1, tp2_cl_size = ifelse(tp2 < tp1, tp1, tp2)) %>% 
      bind_rows(fset, .)
  }
  return(fset)
}

checkChanges <- function(tp, pre_h, post_h, merge_cl) {
  b1 <- tp[h == pre_h] %>% select(iso, cl) %>% set_colnames(c("iso", "pre_h"))
  b2 <- tp[h == post_h] %>% select(iso, cl) %>% set_colnames(c("iso", "post_h"))
  
  a1 <- left_join(b1, b2, by = "iso") %>% filter(pre_h != post_h) %>% 
    select(-iso) %>% unique() %>% arrange(pre_h) %>% as_tibble()
  a2 <- merge_cl %>% arrange(pre_h) %>% as_tibble()
  all_equal(a1, a2) %>% return()
}

absorbingClusters <- function(tp, pre_h, post_h) {
  tp[h == post_h]$cl <- tp[h == pre_h]$cl
  pre_clusters <- tp[h == pre_h]$cl %>% unique()
  changing <- ceiling(length(pre_clusters)/3)
  tomerge <- sample(pre_clusters, changing)
  mergeinto <- sample(setdiff(pre_clusters, tomerge), changing)
  
  x <- tp[h %in% post_h & cl %in% tomerge]$cl
  tp[h %in% post_h & cl %in% tomerge]$cl <- mergeinto[match(x, tomerge)]
  return(list("merge_clusters" = tibble(tomerge, mergeinto), "tp" = tp))
}

simulatedAssignments <- function(ntph, nisos, nclusts) {
  tp <- matrix(data = 0, nrow = nisos, ncol = ntph + 1) %>% 
    as.data.table() %>% 
    set_colnames(c("isolate", paste0("h", 1:ntph)))
  
  tp[,"isolate"] <- paste0("isolate", 1:nisos)
  tp <- melt(tp, id.vars = "isolate") %>% as.data.table() %>% 
    set_colnames(c("iso", "h", "cl")) %>% mutate(across(h, as.character))
  
  # height 1 to height 2 -----------------------------------------------
  tp[h == "h1"]$cl <- length(pull(tp[h == "h1"], "cl")) %>% 
    sample(1:nclusts, ., replace = TRUE)
  
  # height 2 to height n -----------------------------------------------
  # starts at h2 because h1-h2 is the base case
  heights <- tp$h %>% unique()
  for (i in 2:length(heights)) {
    pre_h <- heights[i-1]
    post_h <- heights[i]
    
    res_i <- absorbingClusters(tp, pre_h, post_h)
    
    tp <- res_i$tp
    merge_cl <- res_i$merge_clusters %>% rename(pre_h = tomerge, post_h = mergeinto)
    assert(paste0("Only required clusters were absorbed from height ", pre_h, " to ", 
                  post_h), checkChanges(tp, pre_h, post_h, merge_cl))
  }
  
  wide_form <- dcast(tp, iso ~ h, value.var = "cl")
  return(wide_form)
}

checkIDs <- function(id_col, vtype, vcol) {
  id_col %>% strsplit(., "_") %>% sapply(., '[[', vcol) %>% 
    gsub(vtype, "", .) %>% as.integer() %>% return()
}

# ----------------------------------------------------------------------------------------
# Functions for test-ecc_functions.R -----------------------------------------------------
fakeRaw <- function(ndrs, n, nc) {
  data.table(th = sample(1:nc, n, replace = TRUE) %>% sort(), 
             dr = rep(1:ndrs, ceiling(n/ndrs))[1:n]) %>% 
    rownames_to_column("Strain") %>% return()
}

fakeSecCluster <- function(ndrs, n, nc) {
  testcase <- fakeRaw(ndrs, n, nc)
  nunique <- length(unique(testcase$dr))
  
  assignments <- data.table(
    dr = testcase$dr %>% unique(), 
    Date = sample(seq(as.Date('2019/01/01'), as.Date('2021/01/01'), by="day"), nunique), 
    Latitude = sample(seq(25.00000, 65.00000, by = 0.00001), nunique), 
    Longitude = sample(seq(-100.00000, 125.00000, by = 0.00001), nunique))
  
  # fake sectionClusters() result object:
  parts <- formatForSectioning(testcase, 5) %>% sectionTypingData() %>% 
    list("drs" = testcase, "results" = .)
  return(list("a" = assignments, "b" = parts, "c" = testcase))
}

fakeStrains <- function(nrows, cnames = NULL) {
  sample_dates <- seq(as.Date('2019/01/01'), as.Date('2021/01/01'), by="day") %>% 
    sample(., nrows) %>% as.character() %>% strsplit(split = "-")
  
  basic_strains <- tibble(
    Strain = paste0("Strain.", 1:nrows), 
    Latitude = sample(seq(25.00000, 65.00000, by = 0.00001), size = nrows, replace = TRUE),
    Longitude = sample(seq(-100.00000, 125.00000, by = 0.00001), size = nrows, replace = TRUE), 
    Day = sapply(sample_dates, "[[", 3) %>% as.double(),
    Month = sapply(sample_dates, "[[", 2) %>% as.double(),
    Year = sapply(sample_dates, "[[", 1) %>% as.double()
  )
  
  if ("Country" %in% cnames) {
    basic_strains %<>% add_column(Country = paste0("Country.", 1:nrow(basic_strains)))
  }
  
  if ("Province" %in% cnames) {
    basic_strains %<>% add_column(Province = paste0("Province.", 1:nrow(basic_strains)))
  }
  
  if ("City" %in% cnames) {
    basic_strains %<>% add_column(City = paste0("City.", 1:nrow(basic_strains)))
  }
  
  return(basic_strains)
}


strainByStrainEpiMelt <- function(dr_matches) {
  df1 <- dr_matches %>% select(strain, v)
  
  dr_matches$strain %>% 
    expand.grid(strain1 = ., strain2 = .) %>% as_tibble() %>% 
    left_join(., df1, by = c("strain1" = "strain")) %>% rename(v1 = v) %>% 
    left_join(., df1, by = c("strain2" = "strain")) %>% rename(v2 = v) %>% 
    mutate(Total.Dist = abs(v2 - v1)) %>% 
    select(strain1, strain2, Total.Dist) %>% return()
}

strainByStrainECC <- function(g_cuts, epi_melt) {
  cut_cluster_members <- g_cuts %>% select(th, strain) %>% 
    pivot_longer(-strain, names_to = "cut", values_to = "cluster") %>%  
    group_by(cut, cluster) %>% 
    summarise(members = list(cur_data()$strain), .groups = "drop")
  
  calculate_s1 <- function(k) {
    epi_melt %>% filter(strain1 %in% k, strain2 %in% k) %>%
      select(value) %>% pull() %>% sum()}
  
  cut_cluster_members %>% 
    mutate(
      s1 = map_dbl(members, calculate_s1),
      cluster_size = map_int(members, length),
      ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))
    ) %>% ungroup() %>% 
    select(cluster, cluster_size, ECC) %>% as.data.table() %>% 
    set_colnames(c("th", "th_Size", "th_ECC")) %>% return()
}


runDrECBC <- function(cluster_asmts, strain_data, tau, gamma) {
  dm_temp <- cluster_asmts %>% select(dr, Date) %>% distMatrix(., "temp", "Date")
  transformed_temp <- transformData2(dm_temp, "temp", min(dm_temp), max(dm_temp)) %>% 
    formatData(., c("dr1","dr2","Temp.Dist"))
  
  dm_geo <- cluster_asmts %>% select(dr, Longitude, Latitude) %>% 
    distMatrix(., "geo", c("Longitude", "Latitude"))
  transformed_geo <- transformData2(dm_geo, "geo", min(dm_geo), max(dm_geo)) %>%
    formatData(., c("dr1","dr2","Geog.Dist"))
  
  transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
  # note, we calculate sizes using this, so do not make it unique
  cluster_x <- strain_data %>% select(th, dr)
  
  epiCollectionByCluster(strain_data, tau, gamma, transformed_dists, tpx = 1, cluster_x) %>% 
    set_colnames(c("th", "th_Size", "th_ECC")) %>% return()
}

runStrainECBC <- function(strain_raws, tau, gamma) {
  temp_dm <- strain_raws %>% select(Strain, Date) %>% distMatrix(., "temp", "Date")
  geo_dm <- strain_raws %>% select(Strain, Longitude, Latitude) %>% 
    distMatrix(., "geo", c("Longitude", "Latitude"))
  
  transformed_temp <- transformData2(temp_dm, "temp", min(temp_dm), max(temp_dm)) %>% 
    formatData(., c("strain1","strain2","Temp.Dist"))
  transformed_geo <- transformData2(geo_dm, "geo", min(geo_dm), max(geo_dm)) %>%
    formatData(., c("strain1","strain2","Geog.Dist"))
  
  transformed_dists <- merge.data.table(transformed_temp, transformed_geo)
  
  epi_melt <- transformed_dists %>% 
    mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
    select(strain1, strain2, Total.Dist) %>% as.data.table() %>% 
    mutate(value = 1 - Total.Dist) %>% select(-Total.Dist)
  
  g_cuts <- strain_raws %>% rename(strain = Strain)
  strainByStrainECC(g_cuts, epi_melt) %>% return()
}

# ----------------------------------------------------------------------------------------
# Functions for test-formatting_cgm.R ----------------------------------------------------
# - for checking sizes and making fakes, stubs, mocks, etc.
checkSizes <- function(res, sz) {
  lapply(1:length(res), function(i) sum(res[[i]]$n) <= sz) %>% 
    unlist() %>% all() %>% return()  
}

fakeRaw <- function(ndrs, n, nc) {
  data.table(th = sample(1:nc, n, replace = TRUE) %>% sort(), 
             dr = rep(1:ndrs, ceiling(n/ndrs))[1:n]) %>% 
    rownames_to_column("Strain") %>% return()
}

fakeSecCluster <- function(ndrs, n, nc) {
  testcase <- fakeRaw(ndrs, n, nc)
  nunique <- length(unique(testcase$dr))
  
  assignments <- data.table(
    dr = testcase$dr %>% unique(), 
    Date = sample(seq(as.Date('2019/01/01'), as.Date('2021/01/01'), by="day"), nunique), 
    Latitude = sample(seq(25.00000, 65.00000, by = 0.00001), nunique), 
    Longitude = sample(seq(-100.00000, 125.00000, by = 0.00001), nunique))
  
  # fake sectionClusters() result object:
  parts <- formatForSectioning(testcase, 5) %>% sectionTypingData() %>% 
    list("drs" = testcase, "results" = .)
  return(list("a" = assignments, "b" = parts, "c" = testcase))
}

# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------