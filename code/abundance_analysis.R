library(ggplot2)
library(dplyr)
library(gridExtra)

#Read raw dataset for isolates
Isolate_data <- read.csv("data/clean/Isolates_11-12_final.csv")

#Subsetting by year
Isolate_11 <-subset(Isolate_data, Isolate_data$Year==2011)
Isolate_12 <-subset(Isolate_data, Isolate_data$Year==2012)

#Summarise data by year
Data_2011 <- ddply(Isolate_11, "Species", summarise,
                   N = as.numeric(length(qDef)),
                   freq = (N/length(Isolate_11$Year))*100
)

Data_2012 <- ddply(Isolate_12, "Species", summarise,
                   N = as.numeric(length(qDef)),
                   freq = (N/length(Isolate_12$Year))*100
)
#Re-merge summarized data by year and replace NAs by 0
Data_11_12 <- full_join(Data_2011, Data_2012, by="Species") %>% 
  dplyr::rename(N_11=N.x, freq_11=freq.x, N_12=N.y, freq_12=freq.y)

Data_11_12[is.na(Data_11_12)] <- 0

#Labeling reported pathogens
spp <- c("Phytophthora megasperma","Phytophthora sansomeana",
         "Phytophthora sojae","Pythium aphanidermatum",
         "Pythium attrantheridium","Pythium debaryanum","Pythium delawarense",
         "Pythium dissotocum","Pythium echinulatum","Pythium helicoides",
         "Pythium inflatum","Pythium irregulare","Pythium myriotylum",
         "Pythium oligandrum","Pythium sylvaticum","Pythium torulosum",
         "Pythium ultimum","Pythium ultimum var. sporangiiferum",
         "Pythium ultimum var. ultimum")

for(i in spp) {
  levels(Data_11_12$Species)[levels(Data_11_12$Species)==i] <- paste(i,"*",sep="")
}

Data_11_12 <- Data_11_12[with(Data_11_12, order(-N_11, Species)),] 