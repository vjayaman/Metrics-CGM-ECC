
# title: "SARS_CoV2 EpiQuant Analysis Report" 
# output:   
#   bookdown::html_document2:
#   number_sections: FALSE

# ```{r setup, include=FALSE}
# # knitr::opts_chunk$set(echo = FALSE,message = FALSE,warning = FALSE,eval=TRUE)
# ```


# ```{r load-packages,warning = FALSE}
libs <- c("openxlsx","reshape2","plyr","dplyr","tidyr","grid","gridExtra", 
          "ggplot2","RColorBrewer","gplots","readr","scatterplot3d","heatmap3",
          "devtools","knitr","ggpubr","ggrepel","stringr","magrittr","tibble",
          "purrr","summarytools","lubridate","tinytex","png","patchwork",
          "tidyverse","ggridges","bookdown","foreign","data.table","forcats",
          "plotly")
y <- lapply(libs, require, character.only = TRUE); rm(y); rm(libs)

## for tables
pacman::p_load(
  rio,          # File import
  here,         # File locator
  skimr,        # get overview of data
  tidyverse,    # data management + ggplot2 graphics 
  gtsummary,    # summary statistics and tests
  rstatix,      # summary statistics and statistical tests
  janitor,      # adding totals and percents to tables
  scales,       # easily convert proportions to percents  
  flextable,    # converting tables to pretty images
  kableExtra)
# ```

# ### <span style="color:blue">Section I: Outline of data analysis</span> {#top-link -} 
# #### <span style="color:blue">1. Basic Information</span>
# *	Name of analyst: Guangzhi, Vasena, Ed, Samantha 
# *  Date of analysis: 2022-01-13
# *	Data source: [GISAID](https://www.gisaid.org)
# 
# #### <span style="color:blue">2. Analytical Parameters</span>
# 1.	Geographical area of interest
# +	Region/Countries/Provinces: 
#   2.	Period of interest
# +	TP0: 2020-07-01
# +	TP1: 2020-11-15
# +	TP2: 2022-03-01
# 3.	Subsampling information
# +	Maximum cluster size= 10000
# 4.	Variants of interest
# +	VOCs/Prevalent lineages/Mutations of interest
# ### <span style="color:blue">Section II: Data analysis</span> {.unnumbered}



## loading data 
ECC_long <- "results_TP/Merged_strain_results/" %>% 
  list.files(path = ., full.names = TRUE, pattern = "tsv") %>% map(read_tsv) %>% reduce(rbind)

# Add Lineage information
Lineages <- read_csv("report_specific/input/GISAID Lineages (770000 isolates).csv") %>% 
  filter (T0 != "#N/A") %>% select("Strain", "T0.original","Pango_lineage")

ECC_long <- merge(x = ECC_long, y = Lineages, by = "Strain", all.x=TRUE)

# create new varibles Timepoints, Multistrains, Week and Running day

ECC_long <- ECC_long %>% 
  mutate(Week = week(Date), 
         Year_Week = strftime(Date, format = "%G-%V"), 
         Year_Month = strftime(Date, format = "%Y-%m"), 
         Running_Day = difftime(Date, min(Date), units = "days"), 
         Multistrain_Cluster = ifelse(`TP2 cluster size (1)` > 2,"Multistrain","Singleton"),
         Timepoint = ifelse(interval == "set0-set1","TP1","TP2" ), 
         Cluster_Size = `TP2 cluster size (1)`,  
         Geo_ECC = TP2_ECC.0.0.1,
         Temp_ECC = TP2_ECC.0.1.0,
         Delta_Geo_ECC = delta_ECC_0.0.1,
         Delta_Temp_ECC = delta_ECC_0.1.0, 
         Cluster_Size1 = `TP1 cluster size + 1 (2)`,
         Geo_ECC1 = TP1_ECC.0.0.1,
         Temp_ECC1 = TP1_ECC.0.1.0,          
         "Cluster growth by size" = `Actual cluster growth (TP2 size - TP1 size)`, 
         "Cluster growth by rate" = `Novel growth = (TP2 size) / (TP2 size - number of novels)`
         )

Top20_clusters <- ECC_long %>% group_by(Timepoint, T0.original) %>% dplyr::summarise(n=n()) %>% top_n(20)
ECC_long <- ECC_long %>% 
  mutate(PANGO_lineage = ifelse(T0.original %in% Top20_clusters$T0.original, T0.original,"Others")) %>% 
  mutate_if(is.numeric, ~ round(., digit = 2))

## loading monthly data 
### need to replace '2020-07' with a variable
# Add Lineage information
ECC_Month <- "results_Monthly/Merged_strain_results/" %>% 
  list.files(path = ., full.names = TRUE, pattern = "tsv") %>% map(read_tsv) %>% reduce(rbind)


# create new varibles Timepoints, Multistrains, Week and Running day
ECC_Month <- ECC_Month %>% 
  mutate(interval = gsub('set0', '2020-07', interval)) %>% 
  separate(interval, into = c('Month_1', 'Month_2'), sep = 8) %>% 
  merge(x = ., y = Lineages, by = "Strain", all.x = TRUE) %>% 
  mutate(Multistrain_Cluster = ifelse(`TP2 cluster size (1)` > 2,"Multistrain","Singleton"), 
         Cluster_Size = `TP2 cluster size (1)`,
         Geo_ECC = TP2_ECC.0.0.1,
         Temp_ECC = TP2_ECC.0.1.0)

## Create figure caption variables 
Fig.1.1.A <- "Accumulated cases in countries"
Fig.1.1.B <- "Predominant lineages by genome count (top20)"

Fig.2.1.A <- "Geospatial and temporal ECC distribution (histogram)"
Fig.2.1.B <- "Geospatial and temporal ECC 3D surface density at TP1 and TP2"

Fig.2.2.A <- "Bubble plots with directionality vectors"
Fig.2.2.B <- "Summary of directionality vectors by TP1 cluster size"
Fig.2.2.C <- "Summary of directionality vectors by cluster growth rate"
Tabble.2.1  <- "ECC change vector classes and possible transmission modes"
Fig.2.2.D <- "Summary of directionality vectors by ECC and Delta ECC level "


Fig.2.3.A <- "Top 20 growth clusters (by size)"
Fig.2.3.B <- "Average geographic and temporal pairwise distance for top 20 clusters (by size)"
Fig.2.3.C <- "Top 20 growth clusters (by growth rate)"
Fig.2.3.D <- "Average geographic and temporal pairwise distance for top 20 clusters (by growth rate)"

Fig.2.4.A <- "Monthly Geo and Temp ECC (average)"
Fig.2.4.B <- "Top 10 clusters (by size) Geo ECC, Geo-Temp ECC,Temp ECC"
Fig.2.4.C <- "Top 10 clusters (by growth rate) Geo ECC,Geo-Temp ECC,Temp ECC"

Fig.2.5.A <- "Weekly Geo and Temp ECC (average)"
Fig.2.5.B <- "Top 10 clusters (by size) Geo ECC, Geo-Temp ECC,Temp ECC"
Fig.2.5.C <- "Top 10 clusters (by growth rate) Geo ECC,Geo-Temp ECC,Temp ECC"

# num_top_clusters <- function()
Fig.3.1.A <- "GEO + 50-50 + TEMP heatmaps for top 5 cluster"
Fig.3.1.B <- "GEO + 50-50 + TEMP heatmaps for top 4 cluster"
Fig.3.1.C <- "GEO + 50-50 + TEMP heatmaps for top 3 cluster"
Fig.3.1.D <- "GEO + 50-50 + TEMP heatmaps for top 2 cluster "
Fig.3.1.E <- "GEO + 50-50 + TEMP heatmaps for top 1 cluster"

Fig.3.1.F <- "GEO + 50-50 + TEMP frequency for top 5 cluster"
Fig.3.1.G <- "GEO + 50-50 + TEMP frequency for top 4 cluster"
Fig.3.1.H <- "GEO + 50-50 + TEMP frequency for top 3 cluster"
Fig.3.1.I <- "GEO + 50-50 + TEMP frequency for top 2 cluster "
Fig.3.1.J <- "GEO + 50-50 + TEMP frequency for top 1 cluster"

Fig.3.2.A <- "GEO + 50-50 + TEMP heatmaps for top 5 cluster"
Fig.3.2.B <- "GEO + 50-50 + TEMP heatmaps for top 4 cluster"
Fig.3.2.C <- "GEO + 50-50 + TEMP heatmaps for top 3 cluster"
Fig.3.2.D <- "GEO + 50-50 + TEMP heatmaps for top 2 cluster "
Fig.3.2.E <- "GEO + 50-50 + TEMP heatmaps for top 1 cluster"

Fig.3.2.F <- "GEO + 50-50 + TEMP pairwise distance frequency for top 5 cluster"
Fig.3.2.G <- "GEO + 50-50 + TEMP pairwise distance frequency for top 4 cluster"
Fig.3.2.H <- "GEO + 50-50 + TEMP pairwise distance frequency for top 3 cluster"
Fig.3.2.I <- "GEO + 50-50 + TEMP pairwise distance frequency for top 2 cluster "
Fig.3.2.J <- "GEO + 50-50 + TEMP pairwise distance frequency for top 1 cluster"

Fig.3.3.A <- "Pairwise geospatial distance for top20 clusters (by size)"
Fig.3.3.B <- "Pairwise temporal distance for top20 clusters (by size) "

Tab.4.1.a <- "SARS-Cov2 genomes by country and TP/months"
Tab.4.1.b <- "SARS-Cov2 genomes by country and week"

Tab.4.2.a <- "SARS-Cov2 genomes by PANGO lineages and country"
Tab.4.2.b <- "SARS-Cov2 genomes by PANGO lineages and TP/month"
Tab.4.2.c <- "SARS-Cov2 genomes by PANGO lineages and week"
Tab.4.2.d <- "Change vector directionality"

Fig.5.1.A.I <-  "TP1 vs TP2 genome counts"
Fig.5.1.A.II <- "Cumulative genome counts by month, faceted by country"
Fig.5.1.A.III <-  "Genome counts, faceted by month"
Fig.5.1.A.IV <- "Cumulative genome counts by week, faceted by country "
Fig.5.1.A.V  <- "Genome counts, faceted by week "

Fig.5.1.B.I  <-  "TP1 vs TP2 genome counts"
Fig.5.1.B.II <- "Cumulative genome counts by month, faceted by country"
Fig.5.1.B.III <-  "Genome counts, faceted by month"
Fig.5.1.B.IV <- "Cumulative genome counts by week, faceted by country "
Fig.5.1.B.V  <- "Genome counts, faceted by week "

Fig.5.2.A.I <-  "Cumulative genome counts for most prevalent lineages, by country"
Fig.5.2.B.I <-  "Cumulative genome counts for most prevalent lineages, by country"

Fig.5.3.A.I <-  "TP1 vs TP2  of countries detected "
Fig.5.3.A.II <- "Countries detected: most prevalent lineages, by month"
Fig.5.3.A.III <-  "Countries detected: most prevalent lineages, by week"

Fig.5.3.B.I <-  "TP1 vs TP2 of	Provinces detected"
Fig.5.3.B.II <- "Provinces detected: most prevalent lineages, by month"
Fig.5.3.B.III <-  "Provinces detected: most prevalent lineages, by week"

Fig.5.4 <- "ECC direction classes vs compass rose"

#### <span style="color:blue">1. Overview of cases</span>  {.unnumbered}
##### &nbsp;<span style="color:blue"> 1.1. Cases and lineages </span> {.unnumbered}
 # Figure 1.1.A. **`r paste0(Fig.1.1.A)`** {#Fig.1.1.A-link -}

## Define the theme for plotting

Plot_Theme = theme(plot.margin = unit(c(1,1,1,1), "cm"),
                   axis.text = element_text( size = 12,face = "bold"),
                   axis.title = element_text( size = 16, face = "bold" ),
                   #legend.position="none",
                   strip.text = element_text(size = 12,face = "bold"),
                   strip.background = element_rect(fill = "lightblue", colour = "black",size = 1),
                   plot.title = element_text(size = 22),
                   plot.tag = element_text(size = 22,face = "bold"))

Tag_Theme_1 = theme(plot.tag = element_text(size =24,face = "bold"),
                    plot.margin = unit(c(1,0,0,1), "cm"),
                    plot.tag.position = c(0.00,1),
                    plot.title = element_text(hjust= -0.03))

Legend_Theme =theme(legend.title = element_text(size=12,face = "bold"),
                    legend.text = element_text(size=10,face = "bold"))

## define font for plotly 
t <- list(family = "Courier New", size = 16, color = "blue")
t1 <- list(size = 16, family = "Arial", color = "black")
t2 <- list(family = "Courier New", size = 16, color = "green")
t3 <- list(family = 'Arial')



## Monthly accumulated case line plot

Country_sum <- ECC_Month %>% group_by(Country,Month_2) %>% summarise(Count = n()) %>% 
  group_by(Country) %>% mutate(Count2 = cumsum(Count))
Country_Order <- ECC_Month %>% group_by(Country) %>% summarise(Count = n()) %>% arrange(desc(Count))
P1 <-ggplot(data = Country_sum, aes(Month_2, y = Count2,  group = Country )) + 
  theme_bw()+
  facet_wrap(~ fct_relevel(Country, Country_Order$Country))+
  geom_line(color= "#0077b6", size =1)+
  labs(y = "Cases", x = "Month" ) +  
  theme(axis.text.x = element_text(size = 10,face = "bold",angle = 90, ,hjust = 1, vjust = 0.1)) +
  Plot_Theme+  Legend_Theme
# P1


##### &nbsp;&nbsp;Figure 1.1.B. **`r paste0(Fig.1.1.B)`**  {#Fig.1.1.B-link -} 
Country_Order <- ECC_long %>% filter(Timepoint == "TP2") %>% group_by(Country) %>% 
  dplyr::summarise(n=n()) %>% arrange(desc(n)) %>% top_n(20)
Country_Order <-Country_Order[1:20,]
arranged_clusters <- ECC_long %>% filter(Timepoint == "TP2") %>% group_by(T0.original) %>% 
  dplyr::summarise(n=n()) %>% arrange(desc(n))
T20_clusters <- arranged_clusters %>% top_n(20)
T10_clusters <- arranged_clusters %>% top_n(10)
T5_clusters  <- arranged_clusters %>% top_n(5)
Variant_prevelance <- ECC_long %>% 
  select ("Strain","Timepoint", "Country","Province", "Year_Month", "Year_Week","PANGO_lineage")
Dt20 <- Variant_prevelance %>% group_by(Country, Timepoint,PANGO_lineage) %>% 
  filter (PANGO_lineage %in% T20_clusters$T0.original, 
          Country %in% Country_Order$Country) %>% 
  mutate(Genome_count=n())

## create full lineages in top 20 countries
mat <- matrix(0, 20, 20)
rownames (mat) <- unique(Dt20$PANGO_lineage)
colnames(mat) <- unique(Dt20$Country)
mat1 <- reshape2::melt(mat)
mat2 <- mat1
mat1$timepoint2 <- "TP1"
mat2$timepoint2 <- "TP2"
mat3 <- rbind(mat1,mat2)
colnames(mat3) <-c("PANGO_lineage","Country","Count","Timepoint")
DT20 <- Dt20 %>% group_by(Country, Timepoint,PANGO_lineage) %>% dplyr::summarise(n=n())
DT20_FUll <- left_join(mat3,DT20)
DT20_FUll[is.na(DT20_FUll)] <- 0
Country_Order <- DT20_FUll %>% filter(Timepoint == "TP2") %>% group_by(Country) %>% 
  dplyr::summarise(count= sum(n)) %>% arrange(desc(count)) %>% 
  rowid_to_column() %>% mutate(group = ifelse(rowid<6, "1",ifelse(rowid <11, "2",ifelse(rowid <16, "3","4"))))
DT20_FUll <- left_join(DT20_FUll,Country_Order) %>% 
  mutate(Country = factor(Country, levels = Country_Order$Country))

# -Faceted histogram by Country
P1 <- ggplot(DT20_FUll %>% filter(group == 1), 
             aes(x = fct_relevel(PANGO_lineage, T20_clusters$T0.original), y = n, fill= Timepoint)) + 
  theme_bw() +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#0077b6","#FCA311"))+
  facet_wrap(~ Country, ncol = 5)+
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.1)) +
  # scale_x_discrete(position = "top") +
  labs(y = "Genome Count", x = " ", fill = "Time point")
#ggtitle("Top 20 clusters-Genome count") 
# P1 + Plot_Theme + Legend_Theme




P2 <- ggplot(data = DT20_FUll %>% filter(group == 2), 
             aes(x = fct_relevel(PANGO_lineage, T20_clusters$T0.original), y = n, fill = Timepoint)) + 
  theme_bw() +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#0077b6","#FCA311"))+
  facet_wrap(~ Country, ncol = 5)+
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.1))+
  labs(y = "Genome Count", x = " ", fill = "Time point")+
  Plot_Theme + Tag_Theme_1 + Legend_Theme
# P2


P3 <- ggplot(data = DT20_FUll %>% filter(group == 3), 
             aes(x = fct_relevel(PANGO_lineage, T20_clusters$T0.original), y = n, fill = Timepoint)) + 
  theme_bw() +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values = c("#0077b6","#FCA311"))+
  facet_wrap(~ Country, ncol = 5)+
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.1)) +
  # scale_x_discrete(position = "top") +
  labs(y = "Genome Count", x = " ", fill = "Time point")+
  Plot_Theme + Tag_Theme_1 + Legend_Theme
# P3


P4 <- ggplot(data = DT20_FUll %>% filter(group == 4), 
             aes(x = fct_relevel(PANGO_lineage, T20_clusters$T0.original), y = n, fill = Timepoint)) + 
  theme_bw() +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values = c("#0077b6","#FCA311"))+
  facet_wrap(~ Country, ncol = 5)+
  theme(axis.text.x = element_text(size = 10,face = "bold", angle = 90, hjust = 1, vjust = 0.1)) +
  labs(y = "Genome Count", x = "PANGO_lineage", fill = "Time point")+
  Plot_Theme + Tag_Theme_1 + Legend_Theme
P4




# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;go back to [Overview of cases](#plot-link1)
#### <span style="color:blue">2. Epidemiological cluster cohesion (ECC) analysis</span> {.unnumbered}
##### &nbsp;<span style="color:blue"> 2.1. ECC histogram and contour plots</span> {.unnumbered}
    
##### &nbsp;&nbsp; Figure 2.1.A. **`r paste0(Fig.2.1.A)`**  {#Fig.2.1.A-link -}
## GEO ECC histogram: 
TP1_ECC <- ECC_long %>% filter(Timepoint=="TP1") %>% select(Cluster_Size, Geo_ECC, Temp_ECC) %>% drop_na() 
TP2_ECC <- ECC_long %>% filter(Timepoint=="TP2") %>% select(Cluster_Size, Geo_ECC, Temp_ECC) %>% drop_na() 

P1 <- ggplot(TP1_ECC, aes(x = Geo_ECC)) +
  theme_bw()+
  geom_histogram(binwidth = 0.02, fill = "#0077b6")+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.00), 
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), expand = c(0,0), lim = c(0,1))+
  ggtitle("TP1")+
  labs(y = "Genome Count", x = "Geo_ECC") + Plot_Theme + Legend_Theme

# TEMP ECC histogram: 
# -Faceted histogram by time interval
P2 <- ggplot(TP1_ECC, aes(x = Temp_ECC)) +
  theme_bw()+
  geom_histogram(binwidth = 0.02,fill = "#0077b6")+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.00), 
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), expand = c(0,0), lim = c(0,1)) +
  ggtitle(" ") +
  labs(y = "Genome Count", x = "Temp_ECC") + Plot_Theme + Legend_Theme

P1 + P2 + plot_layout(widths = c(1, 1))



## GEO ECC histogram: 
# TP2_ECC <- ECC_long%>%filter (Timepoint=="TP2")%>% select(Cluster_Size, Geo_ECC, Temp_ECC)%>%drop_na() 

P3 <- ggplot(TP2_ECC, aes(x = Geo_ECC)) +
  theme_bw()+
  geom_histogram(binwidth = 0.02, fill = "#FCA311")+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.00), 
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), expand = c(0,0), lim = c(0,1)) +
  ggtitle("TP2") +
  labs(y = "Genome Count", x = "Geo_ECC") + Plot_Theme + Legend_Theme

# TEMP ECC histogram: 
# -Faceted histogram by time interval
P4 <- ggplot(TP2_ECC, aes(x = Temp_ECC)) +
  theme_bw() +
  geom_histogram(binwidth = 0.02, fill = "#FCA311")+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.00), 
                     labels=c("0.00", "0.25", "0.50", "0.75", "1.00"), expand = c(0,0), lim = c(0,1))+
  ggtitle(" ") +
  labs(y = "Genome Count", x = "Temp_ECC") + Plot_Theme + Legend_Theme

P3 + P4 + plot_layout(widths = c(1, 1))


##### &nbsp;&nbsp;Figure 2.1.B. **`r paste0(Fig.2.1.B)`** {#Fig.2.1.B-link -}

## TP1 Contour
TP1_ECC <- ECC_long %>% filter (Timepoint=="TP1") %>% select(Geo_ECC, Temp_ECC) %>% drop_na() 

kd1 <- with(TP1_ECC, MASS::kde2d(Geo_ECC, Temp_ECC, n = 50))

P1 <- plot_ly(x = kd1$x, y = kd1$y, z = kd1$z, scene='scene1') %>% 
  add_surface(colorbar = list(title ='TP1 Density (z)'), reversescale = T) %>%
  layout(title = list(text = "</b> TP1 (x= Geo_ECC, y= Temp_ECC, z= Density)", font = t1, y = 1.1, x = 0.05),
         scene = list(camera = list(x=1.87, y=0.88, z=0.2)))

## TP2 Contour
TP2_ECC <- ECC_long %>% filter(Timepoint=="TP2") %>% select(Geo_ECC, Temp_ECC) %>% drop_na() 
kd2 <- with(TP2_ECC, MASS::kde2d(Geo_ECC, Temp_ECC, n = 50))

P2 <- plot_ly(x = kd2$x, y = kd2$y, z = kd2$z, scene = 'scene2') %>% 
  add_surface(colorbar = list(title ='TP2 Density (z)'), reversescale = T) %>%
  layout(title= list(text = "</b> TP1 (left) vs TP2 (Geo_ECC(x),Temp_ECC (y), Density(z))",font = t1, y = 1.1, x = 0.05),
         scene = list(camera = list(x=1.87, y=0.88, z=-0.2)))
# 
# subplot(P1,P2)#  %>%
#   layout(showlegend = TRUE, legend = list(font = list(size = 30)))
# 



# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;go back to [Epidemiological cluster cohesion (ECC) analysis](#plot-link2)
##### &nbsp;<span style="color:blue"> 2.2.	ECC Change vector directionality analysis</span> {.unnumbered}

### create bubble plot data
Bubble_Plot_data2 <- ECC_long %>% group_by(Timepoint,T0.original) %>% 
  summarise_at(vars(Geo_ECC,Temp_ECC,Cluster_Size,Cluster_Size1, Geo_ECC1,
                    Temp_ECC1,`Cluster growth by size`,`Cluster growth by rate`), mean)
df3 <- Bubble_Plot_data2 %>% filter(Timepoint == "TP2", Cluster_Size > 1) %>% 
  mutate(Timepoint = "TP1", Cluster_Size = Cluster_Size1-1) %>% 
  drop_na() %>% rename_at(vars(c("Geo_ECC","Temp_ECC", "Geo_ECC1","Temp_ECC1")), 
                          ~c("Geo_ECC1","Temp_ECC1", "Geo_ECC","Temp_ECC")) %>% 
  mutate(TP1_Size = ifelse (Cluster_Size == 0,"TP1_size = 0", ifelse (Cluster_Size == 1,"TP1_size = 1",ifelse (Cluster_Size> 1 & Cluster_Size < 11 ,"1< TP1_size < 11", 
                                                                                                               ifelse (Cluster_Size > 100,"TP1_size > 100", "10< TP1_size < 101"))))) %>% as.tibble()
    
df2 <- df3 %>% select("T0.original","TP1_Size")
TP2_Bubble <- Bubble_Plot_data2 %>% filter(Timepoint=="TP2") %>% drop_na() %>% left_join(df2)

Bubble_Plot_data1 <- rbind(df3, TP2_Bubble) 

## TP1+TP2 plot
P3 <- ggplot(Bubble_Plot_data1, 
             aes(Geo_ECC, y = Temp_ECC,lable1 = T0.original,size = Cluster_Size, 
                 lable12 = `Cluster growth by size`,lable13 = `Cluster growth by rate`, fill = Timepoint)) +
  theme_bw() +
  geom_point(alpha = 0.8, shape = 21, color = "black",position = "jitter") +
  scale_fill_manual(values = c("#0077b6","#FCA311")) +
  scale_size_continuous(limits = c(0, max(Bubble_Plot_data1$Cluster_Size)), 
                        range = c(1,17), breaks =  c(10,100,500,1000,3000,6000))+
  xlim(0, 1)+
  ylim(0,1)+
  # geom_text(aes(label = cluster_size_2),size=2,  hjust = 0.5)+
  labs(y = "Temp_ECC", x = "Geo_ECC", tag = "C", fill = " ", size =" ")
layout_ggplotly <- function(gg, x = -0.04, y = -0.04){
  # The 1 and 2 goes into the list that contains the options for the x and y axis labels respectively
  gg[['x']][['layout']][['annotations']][[1]][['y']] <- x
  gg[['x']][['layout']][['annotations']][[2]][['x']] <- y
  gg
}
ggplotly(P3, height = 800, width = 1000) %>%layout( margin = list(l = 40),title= list(text = "TP1 & TP2",font = t1, y = 1.1, x = 0.05), font=t3,xaxis = list(title=list(text='<b>Geo_ECC</b>', font = list(size = 16), standoff = 10)),yaxis = list(title=list(text='<b>Temp_ECC</b>',font = list(size = 16), standoff = 10)),legend = list(title=list(text = "<b>Time point</b>",y = 0.8, x = 1.0), x=1.0,y =0.75))

    
    ##### &nbsp;&nbsp; Figure 2.2.A. **`r paste0(Fig.2.2.A)`**  {#Fig.2.2.A-link -}
    
    ```{r error=FALSE, warning=FALSE, message=FALSE,fig.height= 9, fig.width= 18}
    ## TP1  Cluster size < 1
    df <- df3%>%filter(Cluster_Size < 2)
    P2data <-Bubble_Plot_data1%>%filter(T0.original%in%df$T0.original)%>% mutate(Color = ifelse(TP1_Size== "TP1_size = 0"  ,"blue", "#ff00ee" ))
    
    P1 <- ggplot(P2data,aes(x=Geo_ECC, y= Temp_ECC,lable1=T0.original,size = Cluster_Size, lable12= `Cluster growth by size`,lable13= `Cluster growth by rate`,fill = TP1_Size)) +
      theme_bw()+
      geom_point(alpha=0.8, shape=21, color= "black")+
      scale_fill_manual(values = c("#0077b6", "#ff00ee"))+
      scale_size_continuous(limits = c(0, max(Bubble_Plot_data1$Cluster_Size)), range = c(1,17), breaks =  c(10,100,500,1000,3000,6000))+
      xlim(0, 1)+
      ylim(0,1)+
      geom_point(data = df,aes(Geo_ECC, y= Temp_ECC))+ 
      labs(y = "Temp_ECC", x = "Geo_ECC",fill = " ",size =" ")
    
    
    ggplotly(P1,height = 800, width = 1000)%>%
      layout( margin = list(l = 10), title= list(text = "</b> TP1 cluster size < 2",font = t1, y = 1.1, x = 0.05),xaxis = list(title=list(text='<b>Geo_ECC</b>', font = list(size = 16), standoff = 10)),yaxis = list(title=list(text='<b>Temp_ECC</b>',font = list(size = 16), standoff = 10)),legend = list(title=list(text = "<b>TP1_size</b>",y = 0.8, x = 1.0), x=1.0,y =0.75))
    ```
    
    
    ```{r error=FALSE, warning=FALSE, message=FALSE}
    ##  1< TP1_size < 11
    
    df<- df3%>%filter(TP1_Size =="1< TP1_size < 11")
    P3data <-Bubble_Plot_data1%>%filter(T0.original%in%df$T0.original)
    
    
    P2 <- ggplot(P3data,aes(x=Geo_ECC, y= Temp_ECC,lable1=T0.original,size = Cluster_Size, lable12= `Cluster growth by size`,lable13= `Cluster growth by rate`,fill= Timepoint)) +
      theme_bw()+
      geom_point(alpha=0.8, shape=21, color= "black")+
      scale_fill_manual(values = c("#0077b6","#FCA311"))+
      scale_size_continuous(limits = c(0, max(Bubble_Plot_data1$Cluster_Size)), range = c(1,17), breaks =  c(10,100,500,1000,3000,6000))+
      xlim(0, 1)+
      ylim(0,1)+
      geom_point(data = df, aes(Geo_ECC, y= Temp_ECC))+ 
      labs(y = "Temp_ECC", x = "Geo_ECC",,fill = " ",size =" ")
    
    ggplotly(P2,height = 800, width = 1000)%>%
      add_annotations( x = df$Geo_ECC1,y = df$Temp_ECC1,
                       xref = "x", 
                       yref = "y",
                       text = "",
                       showarrow = T,
                       arrowcolor='grey',
                       arrowhead = 3,
                       arrowsize = 1,
                       arrowwidth= 1,
                       ax = ~ df$Geo_ECC,
                       ay = ~ df$Temp_ECC,
                       axref = "x", 
                       ayref = "y",
                       alpha =0.02) %>%layout( margin = list(l = 10),title= list(text = "</b> 1< TP1_size < 11",font = t1, y = 1.1, x = 0.05),xaxis = list(title=list(text='<b>Geo_ECC</b>', font = list(size = 16), standoff = 10)),yaxis = list(title=list(text='<b>Temp_ECC</b>',font = list(size = 16), standoff = 10)),legend = list(title=list(text = "<b>Time point</b>",y = 0.8, x = 1.0), x=1.0,y =0.75))
    
    ```
    
    
    

















