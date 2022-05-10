
library(tidyverse)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(ISOweek)
library(stringr)

###  EU_China_USA_CA_2020_Jan_to_June_new  input
#EU_China_USA_CA_2020_Jan_to_June_new <- read_csv("EU_China_USA_CA_2020_Jan_to_June _new.csv")
DT <-  read_csv("GSAID770000_2020.csv")
     
DT <- DT%>%filter( Region == "Europe")%>%mutate(Date= make_date(Year, Month, Day))

DT <- DT%>%filter(Date > as.Date("2020-06-30")& Date < as.Date("2021-03-01"))%>%
   mutate(TP1= ifelse(Date < "2020-10-16","1","0"), TP2 ="1")

DT$"Year_Week" <- paste0(isoyear(DT$Date),"-", formatC(isoweek(DT$Date), format = "f", digits = 0, width = 2, flag = "0"))
EU_Month <-DT%>% group_by(Year_Month) %>% summarise(Sum = n())%>%mutate(Csum = cumsum(Sum))

EU_Week <- DT%>%group_by(Year_Week) %>% summarise(Sum = n())%>%mutate(Csum = cumsum(Sum))
write.table(DT[,1:12], file = "/Users/guzhang/Desktop/July_Feb_Eu/July_Feb/Metrics-CGM-ECC-main/inputs/strain_info.txt", sep = "\t",row.names = FALSE, quote = FALSE)
write.table(DT%>%select("Strain", "T0","T1", "T2", "T3","T4","T5"), file = "/Users/guzhang/Desktop/July_Feb_Eu/July_Feb/Metrics-CGM-ECC-main/inputs/tp2_clusters_init.txt", sep = "\t",row.names = FALSE, quote = FALSE)

### Europe 2020 data input
EU_2020 <- data%>%filter(Region =="Europe", Year == 2020)%>%na.omit()
China_2020 <- data%>%filter(Country =="China", Year == 2020)%>%na.omit()
EU_China_2020<- rbind(EU_2020,China_2020)
EU_China_2020$TP2 = "1"
EU_China_2020<- EU_China_2020 %>% mutate(TP1 = ifelse( Month < 4,"1","0"))

dir.create("Eu_china_2020_input", showWarnings = FALSE)
## srtain data
write.table(EU_China_2020[,1:12], file = "Eu_china_2020_input/strain_info.txt", sep = "\t",row.names = FALSE, quote = FALSE)
## TP1 cluster
write.table(EU_China_2020%>%filter(TP1 == 1)%>%select("Strain", "T0","T1", "T2", "T3","T4","T5"), file = "Eu_china_2020_input/tp1_clusters_init.txt", sep = "\t",row.names = FALSE, quote = FALSE)
## TP2 Cluster
write.table(EU_China_2020%>%select("Strain", "T0","T1", "T2", "T3","T4","T5"), file = "Eu_china_2020_input/tp2_clusters_init.txt", sep = "\t",row.names = FALSE, quote = FALSE)


### combine the data (192000 cases)
Total <-rbind(USA_01_06_2020,Canada_01_06_2020 )%>% rbind(China_01_06_2020)%>% rbind(EU_01_06_2020)
Total <-rbind(Total,China_01_06_2020) ## 
USA_China_Eu_Canada_01_06_2020 <-rbind(Total,EU_01_06_2020)%>%left_join(GISAID_Lineages_770000_isolates_)
USA_China_Eu_Canada_01_06_2020 <-USA_China_Eu_Canada_01_06_2020%>%group_by(Latitude,Longitude,Day, Month, Year)%>%mutate(Sum_GT=n())%>%group_by(T0.original)%>%mutate(Sum_cluster=n())
USA_China_Eu_Canada_01_06_2020 <-USA_China_Eu_Canada_01_06_2020%>% mutate(Subset = ifelse( Strain%in% USA_China_Eu_Canada_01_06_2020_subset$Strain,"1","0"))
USA_China_Eu_Canada_01_06_2020$TP2 = "1"
USA_China_Eu_Canada_01_06_2020 <- USA_China_Eu_Canada_01_06_2020 %>% mutate(TP1 = ifelse( Month < 4,"1","0"))
write.csv(USA_China_Eu_Canada_01_06_2020,"USA_China_Eu_Canada_01_06_2020.csv")

### subset the combined data (38000 cases)

USA_China_Eu_Canada_01_06_2020_subset <-USA_China_Eu_Canada_01_06_2020 %>%distinct(Latitude,Longitude,Day, Month, Year, .keep_all = T)%>%group_by(T0.original)%>%mutate(TP2_cluster=n())%>%mutate(TP2_cluster=n())%>%rowid_to_column()# 38000 cases
write.csv(USA_China_Eu_Canada_01_06_2020_subset,"USA_China_Eu_Canada_01_06_2020_subset.csv")

### TP1 March-2020

USA_China_Eu_Canada_01_06_2020_subset<- USA_China_Eu_Canada_01_06_2020_subset  %>% mutate(Date = make_datetime(Year, Month, Day))
TP1_subset <-USA_China_Eu_Canada_01_06_2020_subset%>%filter(Month < 4)%>%group_by(T0.original)%>%mutate(TP1_cluster=n())

### write text file input for ECC calculation 
data <-USA_China_Eu_Canada_01_06_2020_subset[,c(2:13,18:23)]
dir.create("U_E_C_C01_06_2020_subset_input", showWarnings = FALSE)
## srtain data
write.table(data[,1:12], file = "U_E_C_C01_06_2020_subset_input/strain_info.txt", sep = "\t",row.names = FALSE, quote = FALSE)
## TP1 cluster
write.table(data%>%filter(TP1 == 1)%>%select("Strain", "T0","T1", "T2", "T3","T4","T5"), file = "U_E_C_C01_06_2020_subset_input/tp1_clusters_init.txt", sep = "\t",row.names = FALSE, quote = FALSE)
## TP2 Cluster
write.table(data%>%select("Strain", "T0","T1", "T2", "T3","T4","T5"), file = "U_E_C_C01_06_2020_subset_input/tp2_clusters_init.txt", sep = "\t",row.names = FALSE, quote = FALSE)

### prepare Samantha-like data

library(data.table)
dt <- fread( "results_EU_CHina USA-CANADA 1-6-2020/Merged_strain_results.tsv")
dt_TP1<- dt%>%select("Strain" ,  "Country" ,"Province" , "City" ,"Year", "Month" , "Day" , "Latitude", "Longitude" ,"TP1" ,"TP1 cluster", "TP2" , "TP2 cluster", "delta_ECC_0.1.0" ,"delta_ECC_0.0.1","Type", "TP1_T0_ECC.0.0.1" ,  "TP1_T0_ECC.0.1.0", "TP1 cluster size (1)" , "TP1 temp average cluster distance (days)",  "TP2 temp average cluster distance (days)")
dt_TP2 <- dt%>%select("Strain" ,  "Country" ,"Province" , "City" ,"Year", "Month" , "Day" , "Latitude", "Longitude" ,"TP1" ,"TP1 cluster", "TP2" , "TP2 cluster", "delta_ECC_0.1.0" ,"delta_ECC_0.0.1","Type",  "TP2_T0_ECC.0.0.1" ,"TP2_T0_ECC.0.1.0","TP2 cluster size (1)", "TP1 temp average cluster distance (days)",  "TP2 temp average cluster distance (days)")
strains_long <- read_csv("strains_long.csv")
## TP1
dt_TP1$strain_latitude_jit <- 0
dt_TP1$strain_longitude_jit <- 0
dt_TP1 <- dt_TP1%>%mutate(strain_date = make_datetime(Year, Month, Day),stain_time_dif = dt_TP1$`TP2 temp average cluster distance (days)`-dt_TP1$`TP1 temp average cluster distance (days)`)
dt_TP1 <-dt_TP1%>%mutate(timepoint=ifelse (TP1== 1,1,2))
dt_TP1$cluster_size_2<-dt_TP1$`TP1 cluster size (1)`
dt_TP1 <-dt_TP1%>%mutate(single_mult =ifelse(cluster_size_2 > 1, "Multi strain clusters","Single strain clusters"))
dt_TP1$timepoint2 <- 1
TP1_data <- dt_TP1[,c(1:16,22:26,28,29,17,18,27)]
TP1_data$key <-0
TP1_data <-TP1_data %>% rename_at(vars(colnames(TP1_data)),~ colnames(strains_long))
## TP2
dt_TP2$strain_latitude_jit <- 0
dt_TP2$strain_longitude_jit <- 0
dt_TP2 <- dt_TP2%>%mutate(strain_date = make_datetime(Year, Month, Day),stain_time_dif = dt_TP1$`TP2 temp average cluster distance (days)`-dt_TP1$`TP1 temp average cluster distance (days)`)
dt_TP2 <-dt_TP2%>%mutate(timepoint=ifelse (TP1== 1,1,2))
dt_TP2$cluster_size_2<-dt_TP2$`TP2 cluster size (1)`
dt_TP2 <-dt_TP2%>%mutate(single_mult =ifelse(cluster_size_2 > 1, "Multi strain clusters","Single strain clusters"))
dt_TP2$timepoint2 <- 2
TP2_data <- dt_TP2[,c(1:16,22:26,28,29,17,18,27)]
TP2_data$key <-0
TP2_data <-TP2_data %>% rename_at(vars(colnames(TP2_data)),~ colnames(strains_long))
EU_China_USA_CA_2020_Jan_to_June <-rbind(TP1_data,TP2_data)
write.csv(EU_China_USA_CA_2020_Jan_to_June,"EU_China_USA_CA_2020_Jan_to_June.csv")

PD2<- dt[c(1:3,7876),]
PD2<-PD2%>% select("Strain", "Country","Province", "City" ,"Year", "Month" , "Day" , "Latitude", "Longitude", "delta_ECC_0.1.0" ,"delta_ECC_0.0.1", "TP1", "TP2","TP1 cluster", "TP2 cluster","TP1 cluster size (1)","TP2 cluster size (1)","TP1_T0_ECC.0.0.1" , "TP2_T0_ECC.0.0.1" ,  "TP1_T0_ECC.0.1.0","TP2_T0_ECC.0.1.0", )
PD2 <-PD2 %>% rename_at(vars(colnames(PD2)),~ colnames(PD))

PD <-PD[,-1]
df_long <- PD1 %>% pivot_longer(cols = -c(1:11,14,15),names_to = c("observation", ".value"),names_sep = "_")
df_long1 <- PD2 %>% pivot_longer(cols = -c(1:11,14,15),names_to = c("observation", ".value"),names_sep = "_")
df_long1 <- PD%>% pivot_longer(everything(),names_to = c(".value", "set"),names_pattern = "(.)(.)")
write.csv(PD, "PD.csv")
df_long1