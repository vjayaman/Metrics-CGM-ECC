# 
# ##########################################################################################
# ######## Function for generating time data from input datafile ###########################
# temp_calc <- function(input_data){
#   #### import data from data file ####
#   timedata <- input_data %>% 
#     # mutate(date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% 
#     select(Date) %>% 
#     pull() %>% 
#     dist(diag=TRUE, upper=TRUE, method = 'euclidean') %>% 
#     as.matrix(nrow = nrow(input_data), ncol=nrow(input_data)) %>% 
#     add(10) %>% 
#     log10()
#     
#   #timedata <- read.delim(file="datafile.txt", header=TRUE, sep="\t")
#   
#   #### Create a new column to contain the concatenated temporal data ####
#   # timedata$date <- NA
#   # timedata$date <- as.Date(paste(timedata$Year, timedata$Month, timedata$Day, sep = "-"))
#   #### Create an empty matrix and populate it with the pairwise distances ####
#   # time_matrix <- matrix(data = NA, nrow=nrow(timedata), ncol=nrow(timedata))
#   # time_matrix <- as.matrix(dist(x=timedata$date, diag=TRUE, upper=TRUE, method = 'euclidean'), nrow=nrow(timedata), ncol=nrow(timedata))
# 
#   #Make 10 days the startpoint for dropping off in similarity - i.e. 10 days = 100% similar still
#   # time_matrix[time_matrix < 10] <- 10
#   # time_matrix <- time_matrix + 10
#   #### Convert all the distances to a log value and normalize them to [0:1] ####
#   # timedata <- log(time_matrix, base = 10) 
#   timedata[timedata == -Inf ] <- 0
#   
#   if(max(timedata) == 0){
#     timedata <- 0
#   } else {
#     timedata <- ((timedata-min(timedata)) / (max(timedata)-min(timedata)))
#   }
#   
#   #### Import row and column names from the original datafile and melt data for easy reading ####
#   matrix_names <- input_data %>% select(Strain) %>% pull()
#   rownames(timedata) <- matrix_names
#   colnames(timedata) <- matrix_names
#   
#   time_melt <- 
#     timedata %>% 
#     melt(varnames = c("Strain.1", "Strain.2"), value.name = "Temp.Dist") %>% 
#     as_tibble()
#   
#   time_melt
# }  
# 
# ##########################################################################################
# ######## Function for generating geography data from input datafile ######################
# geog_calc <- function(geogdata){
#   
#   #### Read data table from project folder - this contains locations and their GPS Coordinates ####  
# 
#   #### Create a matrix containing the pair-wise distances (in km) between all the locations using the fossil package ####
#   geog_matrix <-
#     geogdata %>% 
#     select(Latitude, Longitude) %>%
#     as.data.frame() %>% 
#     earth.dist(dist=TRUE) %>% 
#     as.matrix() %>%
#     add(10) %>%
#     log10()
#   
#   #### Calculate the maximum distance in the matrix and divide all values by it to arrive at a max distance of 1 ####
#   # geog_matrix[geog_matrix < 10] <- 10
#   # geog_matrix <- geog_matrix + 10
#   # geog_matrix <- log10(geog_matrix)
#   
#   if(max(geog_matrix) == 1){
#       geog_matrix[1:nrow(geog_matrix), 1:nrow(geog_matrix)] <- 0
#     } else {
#       geog_matrix <- ((geog_matrix-min(geog_matrix)) / (max(geog_matrix)-min(geog_matrix)))
#     }
#   #### Import row and column names from the original datafile ####
#   matrix_names <- geogdata %>% select(Strain) %>% pull()
#   colnames(geog_matrix) <- matrix_names
#   rownames(geog_matrix) <- matrix_names
#   
#   #### Create a melted pairwise distance table for easier readability ####
#   geog_melt <- 
#     geog_matrix %>% 
#     melt(varnames = c("Strain.1", "Strain.2"), value.name = "Geog.Dist") %>%
#     as_tibble()
#   
#   geog_melt
# }

##########################################################################################################
######## Function to return table of all the epi-similarities and final strain similarity ################
# datafile <- strain_data; source_matrix <- source_pw; source_coeff <- sigma; temp_coeff <- tau
# geog_coeff <- gamma; geog_temp <- geog_temp
EpiTable <- function(datafile, source_matrix, source_coeff, temp_coeff, geog_coeff, geog_temp){
  #### Read data into memory from previous outputs ####

  # geog_matrix <- geog_calc(datafile)
  # temp_matrix <- temp_calc(datafile)
  
  # geog_temp <- left_join(geog_calc(datafile), temp_calc(datafile))
  
  x <- source_coeff 
  y <- temp_coeff
  z <- geog_coeff
  
  
  #### Create the pairwise table for lookups ####
  d <- expand.grid(1:nrow(datafile), 1:nrow(datafile))
  #### Create Empty matrix ####
  #strain_sims <- matrix(ncol=8, nrow=(nrow(d)))
  #### Make concatenated strings from the date and location data for matching to the similarity scores: ####
  # datafile$Date <- as.Date(paste(datafile$Year, datafile$Month, datafile$Day, sep = "-"))
  # datafile$Location <- as.character(paste(datafile$Country, datafile$Province, datafile$City, sep = "_"))
  
  d1 <- d[,1]
  d2 <- d[,2]
  
  strain_sims <- tibble(
    "Strain.1" = datafile[d1, "Strain"]  %>% pull (),
    "Strain.2" = datafile[d2, "Strain"]  %>% pull (),
    "Source.1" = datafile[d1, "Source"]  %>% pull (),
    "Source.2" = datafile[d2, "Source"]  %>% pull (),
    "Date.1"   = datafile[d1, "Date"]  %>% pull (),
    "Date.2"   = datafile[d2, "Date"]  %>% pull (),
    "Location.1" = datafile[d1, "Location"]  %>% pull (),
    "Location.2" = datafile[d2, "Location"] %>% pull()
    )
  
  # #### Start populating matrix with pairwise strain names and their sources ####
  # strain_sims[,1] <- as.character(datafile[(d[,1]),1])
  # strain_sims[,2] <- as.character(datafile[(d[,2]),1])
  # strain_sims[,3] <- as.character(datafile[(d[,1]),2])
  # strain_sims[,4] <- as.character(datafile[(d[,2]),2])
  # 
  # 

  
  #### Populate the matrix with the date and location pairwise data and rename columns for readability: ####
  # strain_sims[,5] <- as.character(datafile[(d[,1]),11])
  # strain_sims[,6] <- as.character(datafile[(d[,2]),11])
  # strain_sims[,7] <- as.character(datafile[(d[,1]),12])
  # strain_sims[,8] <- as.character(datafile[(d[,2]),12])
  # colnames(strain_sims) <- c("Strain.1", "Strain.2", "Source.1", "Source.2", "Date.1", "Date.2", "Location.1", "Location.2")
  
  #### Lookup and merge source data : ####
  # strain_sims <- merge.data.frame(strain_sims, source_matrix, 
  #                                 by.x = c("Source.1", "Source.2"), 
  #                                 by.y= c("Var1", "Var2"))

  
  str.matrix <-
    strain_sims %>% 
    left_join(source_matrix, by = c("Source.1", "Source.2")) %>% 
    left_join(geog_temp, by = c("Strain.1", "Strain.2"))
    # left_join(temp_matrix, by = c("Strain.1", "Strain.2")) %>% 
    # left_join(geog_matrix, by = c("Strain.1", "Strain.2")) 
  
  # split into two steps, since it seems to seems to reduce memory usage  
  if (x == 0) {
    # This is necessary (otherwise any NA values in the Source.Dist column make Total.Dist and Epi.Sym NA as well. 
    # There is a better way to handle this; will work on that.)
    str.matrix <-
      str.matrix %>% 
      mutate(
        Total.Dist = sqrt( ((Temp.Dist^2)*y) + ((Geog.Dist^2)*z) ),
        Epi.Sym = 1 - Total.Dist
      )    
  }else {
    str.matrix <-
      str.matrix %>% 
      mutate(
        Total.Dist = sqrt( (((Source.Dist^2)*x) + ((Temp.Dist^2)*y) + ((Geog.Dist^2)*z)) ),
        Epi.Sym = 1 - Total.Dist
      ) 
  }
  
  # print(head(strain_sims))
  # strain_sims <- strain_sims[, c(3,4,1,2,5,6,7,8,9)]
  # colnames(strain_sims) <- c("Strain.1", "Strain.2", "Source.1", "Source.2", "Date.1", "Date.2", "Location.1", "Location.2", "Source.Dist")
  # print(head(strain_sims))
  
  
  #### Lookup and merge temporal data: ####
  # strain_sims <- merge.data.frame(strain_sims, temp_matrix, by.x= c("Strain.1", "Strain.2"), by.y = c("Var1", "Var2")) 
  # colnames(strain_sims) <- c("Strain.1", "Strain.2", "Source.1", "Source.2", "Date.1", "Date.2", "Location.1", "Location.2", "Source.Dist", "Temp.Dist")
  
  #### Lookup and merge Geography data : ####
  # strain_sims <- merge.data.frame(strain_sims, geog_matrix, by.x= c("Strain.1", "Strain.2"), by.y = c("Var1", "Var2"))
  # colnames(strain_sims) <- c("Strain.1", "Strain.2", "Source.1", "Source.2", "Date.1", "Date.2", "Location.1", "Location.2", "Source.Dist", "Temp.Dist", "Geog.Dist")
  # 
  # #### Finalize the similarity matrix and calculate the overall similarity between the strains: ####
  # # str.matrix <- strain_sims
  #str.matrix$Total.Dist <- NA
  # str.matrix$Total.Dist <- sqrt( (((str.matrix$Source.Dist^2)*x) + ((str.matrix$Temp.Dist^2)*y) + ((str.matrix$Geog.Dist^2)*z)) )
  # str.matrix$Epi.Sym <- 1 - str.matrix$Total.Dist
  str.matrix
}  

# datafile <- strain_data; source_matrix <- source_pw; source_coeff <- sigma; temp_coeff <- tau
# geog_coeff <- gamma; geog_temp <- matched
EpiTable2 <- function(datafile, source_matrix, source_coeff, temp_coeff, geog_coeff, geog_temp){
  #### Read data into memory from previous outputs ####
  
  x <- source_coeff 
  y <- temp_coeff
  z <- geog_coeff
  
  #### Create the pairwise table for lookups ####
  d <- expand.grid(1:nrow(datafile), 1:nrow(datafile))
  d1 <- d[,1]
  d2 <- d[,2]
  
  # split into two steps, since it seems to seems to reduce memory usage  
  if (x == 0) {
    strain_sims <- tibble(
      "Strain.1" = datafile[d1, "Strain"]  %>% pull (),
      "Strain.2" = datafile[d2, "Strain"]  %>% pull (),
      "Date.1"   = datafile[d1, "Date"]  %>% pull (),
      "Date.2"   = datafile[d2, "Date"]  %>% pull (),
      "Location.1" = datafile[d1, "Location"]  %>% pull (),
      "Location.2" = datafile[d2, "Location"] %>% pull()
    )
    
    # This is necessary (otherwise any NA values in the Source.Dist column make Total.Dist and Epi.Sym NA as well. 
    # There is a better way to handle this; will work on that.)
    str.matrix <- strain_sims %>% 
      left_join(geog_temp, by = c("Strain.1", "Strain.2")) %>% 
      mutate(
        Total.Dist = sqrt( ((Temp.Dist^2)*y) + ((Geog.Dist^2)*z) ),
        Epi.Sym = 1 - Total.Dist
      )    
  }else {
    strain_sims <- tibble(
      "Strain.1" = datafile[d1, "Strain"]  %>% pull (),
      "Strain.2" = datafile[d2, "Strain"]  %>% pull (),
      "Source.1" = datafile[d1, "Source"]  %>% pull (),
      "Source.2" = datafile[d2, "Source"]  %>% pull (),
      "Date.1"   = datafile[d1, "Date"]  %>% pull (),
      "Date.2"   = datafile[d2, "Date"]  %>% pull (),
      "Location.1" = datafile[d1, "Location"]  %>% pull (),
      "Location.2" = datafile[d2, "Location"] %>% pull()
    )
    
    str.matrix <-
      strain_sims %>% 
      left_join(source_matrix, by = c("Source.1", "Source.2")) %>% 
      left_join(geog_temp, by = c("Strain.1", "Strain.2")) %>% 
      mutate(
        Total.Dist = sqrt( (((Source.Dist^2)*x) + ((Temp.Dist^2)*y) + ((Geog.Dist^2)*z)) ),
        Epi.Sym = 1 - Total.Dist
      ) 
  }
  
  str.matrix
}  

##########################################################################################################
######## Function to return matrix of just strains and final similarity scores for building graphics #####
# epi.matrix <- str.matrix
EpiMatrix <- function(epi.matrix){

  # epi.matrix <- epi.matrix[,c(1,2,12)] 
  
  epi.matrix <- epi.matrix %>% # I think, TODO: check
    select(Strain.1, Strain.2, Total.Dist)
  
  epi.cast <- dcast(epi.matrix, formula= Strain.1 ~ Strain.2, value.var = "Total.Dist")
  epi.cast <- as.matrix(epi.cast[,2:ncol(epi.cast)]) 
  rownames(epi.cast) <- colnames(epi.cast)
#   epi.sym <- 1 - epi.cast
  
  epi.cast
}

# ##########################################################################################
# ######## Function to return a heatmap of the final EPIMATRIX function ####################
# EpiHeatmap_d3 <- function(m){
# #   heatcolor<- colorRampPalette(c("darkgreen","yellowgreen","white"))(512)
#   heatcolor<- colorRampPalette(c("white","yellowgreen","darkgreen"))(512)
#   d3heatmap(m, dendrogram = 'both', colors=rev(heatcolor), Rowv = T, 
#             reorderfun = function(d, w) rev(reorder(d, w)),
#             revC=TRUE, hclustfun = function(x) hclust(x,method = 'single'))
# }
# 
# EpiHeatmap_pdf <- function(m){
#   heatcolor<- colorRampPalette(c("white","yellowgreen","darkgreen"))(512)
#   # heatcolor<- colorRampPalette(c("#efedf5", "#bcbddc", "#756bb1"))(512)
#   # heatcolor<- colorRampPalette(c('#eff3ff','#bdd7e7','#6baed6','#3182bd','#08519c'))(512)
#   plot <- heatmap.2(m, col=rev(heatcolor), Rowv = T, Colv = 'Rowv', trace='none',
#             srtCol = 45, key.title = NA, key.ylab=NA,
#             revC=T, margins = c(10,10), keysize = 1.3, key = T,
#             xlab=NULL, ylab=NULL, 
#             labRow = NA, labCol = NA,
#             hclustfun = function(x) hclust(x,method = 'single'))
#   # data <- m[plot$rowInd, plot$colInd]
#   # return(list(plot, data))
# }


