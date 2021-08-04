
# Comma-delimited string of heights to collect ECCs for, e.g. '50,75,100' 
hx <- as.character(params$th[2]) %>% strsplit(., split = ",") %>% unlist() %>% 
  tibble(h = ., th = paste0("T", .))

combos <- params$coeffs[2] %>% strsplit(., ",") %>% unlist()

# Everything else ----------------------------------------------------------------------------------------------
cgms <- suppressMessages(read_tsv(arg$CGMs))

strain_data <- suppressMessages(read_tsv(arg$metadata)) %>% 
  mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
         Location = paste(Country, Province, City, sep = "_"))

# NON-REDUNDANT METHOD
size_details <- cgms %>% arrange(-tp2_cl_size) %>% select(Strain, tp2_cl, tp2_cl_size) %>%
  rename(TP2_cluster_size = tp2_cl_size)

top_clusters <- size_details %>% select(-Strain) %>% unique() %>%
  arrange(-TP2_cluster_size) %>% filter(TP2_cluster_size < 10000) %>%
  slice(1:as.integer(params$numcl[2])) %>% pull(tp2_cl)

clusters <- top_clusters %>% gsub("c", "", .) %>% as.integer()

tp1 <- Timepoint$new(arg$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(arg$tp2, "tp2")$Process(hx)$listHeights(hx)
typing_data <- tp1$height_list %>% append(tp2$height_list)

m <- suppressMessages(read_tsv(arg$metadata)) %>% processedStrains()

assignments <- typing_data %>% extract2(2) %>% set_colnames("tp_cl") %>% 
  rownames_to_column("Strain") %>% filter(tp_cl %in% clusters) %>% as.data.table()

cl_drs <- lapply(1:length(clusters), function(i) {
  x1 <- assignments[tp_cl %in% clusters[i]] %>% pull(Strain)
  m$dr_matches %>% filter(Strain %in% x1) %>% pull(dr) %>% unique()
}) %>% set_names(clusters)

tau <- list("set1" = getCoeff(combos[1], 2), "set2" = getCoeff(combos[2], 2))
gamma <- list("set1" = getCoeff(combos[1], 3), "set2" = getCoeff(combos[2], 3))

epi_tables <- lapply(1:length(clusters), function(j) {
  clx <- clusters[j]
  
  rawdists <- lapply(1:length(distfiles), function(i) {
    drs <- cl_drs[[as.character(clx)]]
    
    tdm <- distfiles[[i]]$temp
    ctdm <- tdm[rownames(tdm) %in% drs, colnames(tdm) %in% drs] %>% 
      transformData2(., "temp", extremes$mint, extremes$maxt) %>% 
      formatData(., c("dr1","dr2","Temp.Dist"))
    
    gdm <- distfiles[[i]]$geo
    cgdm <- gdm[rownames(gdm) %in% drs, colnames(gdm) %in% drs] %>% 
      transformData2(., "geo", extremes$ming, extremes$maxg) %>%
      formatData(., c("dr1","dr2","Geog.Dist"))
    
    merge.data.table(ctdm, cgdm)
  }) %>% bind_rows()
  
  x1 <- assignments[tp_cl %in% clx] %>% pull(Strain)
  x2 <- m$dr_matches %>% filter(Strain %in% x1) %>% as.data.table()
  x3 <- inner_join(x2, rawdists, by = c("dr" = "dr1")) %>% select(-dr) %>% rename(Strain.1 = Strain)
  
  dists <- inner_join(x3, x2, by = c("dr2" = "dr")) %>% rename(Strain.2 = Strain) %>% 
    select(Strain.1, Strain.2, Temp.Dist, Geog.Dist) %>% 
    mutate(Total.Dist.Set1 = sqrt( (((Temp.Dist^2)*tau$set1) + ((Geog.Dist^2)*gamma$set1)) ),
           Epi.Sym.Set1 = 1 - Total.Dist.Set1, 
           Total.Dist.Set2 = sqrt( (((Temp.Dist^2)*tau$set2) + ((Geog.Dist^2)*gamma$set2)) ),
           Epi.Sym.Set2 = 1 - Total.Dist.Set2) %>% unique()
  
  return(dists)
}) %>% set_names(top_clusters)

saveRDS(epi_tables, "report_specific/heatmaps/epitables_for_heatmaps.Rds")

rm(typing_data)
rm(m)
rm(cl_drs)
# epi.matrix <- lapply(1:length(clusters), function(j) EpiMatrix(epi_tables[[j]]))

for (i in length(top_clusters):1) {
  cl_strains <- size_details %>% filter(tp2_cl %in% top_clusters[i])
  if (length(cl_strains$Strain) > 1) {
    for (j in 1:2) { # the two parameter sets (e.g. 0-1-0, 0-0-1)
      epi_matrix <- epi_tables[[top_clusters[i]]] %>% 
        select(Strain.1, Strain.2, Temp.Dist, Geog.Dist, 
               paste0("Total.Dist.Set", j), paste0("Epi.Sym.Set", j)) %>% 
        set_colnames(gsub(paste0(".Set", j), "", colnames(.)))
      cl_epi <- EpiMatrix(epi_matrix)
      cl_id <- cgms %>% filter(tp2_cl == top_clusters[i]) %>% slice(1) %>% pull(first_tp2_flag)
      
      png(paste0("report_specific/heatmaps/", cl_id, ".png"))
      EpiHeatmap_pdf(cl_epi)
      dev.off()
      rm(cl_epi)
    }
  }else {
    print(paste0("TP2 cluster ", top_clusters[i], " has only one strain, no heatmap generated"))
  }
}