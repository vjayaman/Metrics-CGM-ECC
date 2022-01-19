heatmapFunction <- function(m, type = "heatmap3", heatcolor, args = TRUE, hc = NULL) {
  if (type == "heatmap.2") {
    if (args) {
      # hr <- hclust(as.dist(t(m)), method="single")
      plotx <- heatmap.2(m, col=rev(heatcolor1), 
                         Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), revC = T, 
                         scale="none", trace="none", margins = c(10,10), 
                         xlab=NULL, ylab=NULL, labRow = NA, labCol = NA, 
                         keysize = 1.3, key = T, key.title = NA, 
                         key.ylab= "Frequency", density.info = 'histogram')
    }else {
      plotx <- heatmap.2(m, col=rev(heatcolor), 
                         Rowv = T, Colv = 'Rowv', revC=T, scale = "none", 
                         trace='none',margins = c(10,10), 
                         xlab=NULL, ylab=NULL, labRow = NA, labCol = NA,
                         keysize = 1.3, key = T, key.title = NA, key.ylab=NA,
                         hclustfun = function(x) hclust(x,method = 'single'))    
    }
  }else if (type == "heatmap3") {
    # m_vec <- m %>% as.vector()
    plotx <- heatmap3(m, col=rev(heatcolor1), 
                      Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), revC = T, 
                      scale = 'none', trace='none', margins = c(10,10),
                      xlab=NULL, ylab=NULL, labRow = NA, labCol = NA, 
                      keysize = 1.3, key = FALSE, #distfun = function(x) as.dist(x), 
                      key.title = NA, key.ylab=NA, 
                      # key.ylab="Frequency", legendfun=function() hist(m_vec, freq = TRUE), 
                      method = "complete", useRaster = TRUE)
  }
  return(plotx)
}

# ecc_functions.R -----------------------------------------------------------------------------------

processedStrains <- function(base_strains) {
  loc_cols <- intersect(c("Country", "Province", "City"), colnames(base_strains)) %>% sort()

  if (length(loc_cols) > 0) {
    strain_data <- base_strains %>%
      mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")),
             Location = do.call(paste, c(base_strains[loc_cols], sep = "_")))

    assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>%
      unique() %>% rownames_to_column("dr") %>% as.data.table()
    dr_matches <- left_join(strain_data, assignments,
                            by = c("Latitude", "Longitude", "Date", "Location")) %>%
      select(Strain, dr)
  }else {
    strain_data <- base_strains %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")))
    assignments <- strain_data %>% select(Date, Latitude, Longitude) %>%
      unique() %>% rownames_to_column("dr") %>% as.data.table()
    dr_matches <- left_join(strain_data, assignments,
                            by = c("Latitude", "Longitude", "Date")) %>% select(Strain, dr)
  }

  list("strain_data" = strain_data,
       "assignments" = assignments,
       "dr_matches" = dr_matches) %>% return()
}

# Outputs the same message in two ways, one is directed to standard output and one to a log file
outputDetails <- function(msg, newcat = FALSE) {
  cat(msg)
  if (newcat) {cat("\n")}
  message(msg)
}

# epi-helper-no-source.R ----------------------------------------------------------------------------

EpiHeatmap_pdf <- function(m){
  heatcolor<- colorRampPalette(c("white","yellowgreen","darkgreen"))(512)
  # heatcolor<- colorRampPalette(c("#efedf5", "#bcbddc", "#756bb1"))(512)
  # heatcolor<- colorRampPalette(c('#eff3ff','#bdd7e7','#6baed6','#3182bd','#08519c'))(512)
  plot <- heatmap.2(m, col=rev(heatcolor), Rowv = T, Colv = 'Rowv', trace='none',
                    srtCol = 45, key.title = NA, key.ylab=NA,
                    revC=T, margins = c(10,10), keysize = 1.3, key = T,
                    xlab=NULL, ylab=NULL, 
                    labRow = NA, labCol = NA,
                    hclustfun = function(x) hclust(x,method = 'single'))
  return(plot)
  # data <- m[plot$rowInd, plot$colInd]
  # return(list(plot, data))
}

distEpiMatrix <- function(epi.matrix, dist_type) {
  epi.cast <- dcast.data.table(epi.matrix, formula = Strain.1 ~ Strain.2, value.var = dist_type)
  epi.cast <- as.matrix(epi.cast[,2:ncol(epi.cast)]) 
  rownames(epi.cast) <- colnames(epi.cast)
  return(epi.cast)
}

# cname in c("Temp.Dist", "Geog.Dist", "Temp_Geo.Dist")
prepEpiMatrix <- function(epitable, cname) {
  epi.matrix <- epitable %>% select(Strain.1, Strain.2, all_of(cname))
  mat <- distEpiMatrix(epi.matrix, cname)  
  rm(epi.matrix); gc()
  return(mat)
}

densityPlot <- function(save_to, plot_type, mat) {
  png(save_to)
  d <- density(mat, adjust = 0.25)
  plot(d, main = paste0(plot_type, " distances density"), xlab = paste0(plot_type, " pairwise distances"))
  dev.off()        
}

# type <- switch(plot_type, "temp" = "Temporal", "geo" = "Geographical", "tempgeo" = "Temp and geo")
frequencyPlot <- function(save_to, plot_type, mat) {
  png(save_to)
  half_size <- mat
  half_size[upper.tri(half_size)] <- 0
  hist(as.vector(half_size), freq = TRUE, main = paste0(plot_type, " distances counts"), 
       xlab = paste0(plot_type, " pairwise distances"))
  dev.off()
}

heatmapPlot <- function(save_to, mat, heatcolor, hc_choice, col_dir) {
  png(save_to)
  heatplot <- heatmap3(
    mat, col=rev(heatcolor), labRow = NA, labCol = NA, 
    Rowv=as.dendrogram(hc_choice), Colv=as.dendrogram(hc_choice), revC = T, scale = 'none', 
    margins = c(10,10), xlab=NULL, ylab=NULL, method = "complete", useRaster = TRUE, 
    # main = paste0("Pairwise distances for geographical data"), 
    legendfun = function() showLegend(legend=c("Low similarity", "High similarity"), col = col_dir))
  heatplot
  dev.off()  
}



