## Map projection

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, 
#then load them into the R session.
# Source: https://gist.github.com/stevenworthington/3178163
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#Loading packages
packages <- c("maps", "mapproj", "maptools", "ggcounty", "ggplot2",
              "magrittr", "plyr", "dplyr", "raster", "grid")

ipak(packages)

#Activate map permit
gpclibPermit()

par(mfrow = c(1,1))

#Read files  with location data for survey and load the file with soybean area planted (Acres)
fields <- read.csv2("data/clean/US_GIS_data.txt", 
                    header = TRUE,sep = "\t", dec = ".")
soy_data2 <- read.table("data/clean/soybean_data.csv", 
                        sep = ",", header = TRUE, colClasses=c(rep("character",19)))

#Reformating data for soy data
soy_data2$counties <- paste(soy_data2$State,soy_data2$County, sep = ",")
soy_data2$region <- paste(soy_data2$State.ANSI,soy_data2$County.ANSI, sep = "")
soy_data2$counties <- tolower(soy_data2$counties)
soy_data2$value <- as.numeric(gsub(",", "", as.character(soy_data2$Value)))
#Remove other data and extracolumn
soy_data2 <- soy_data2[soy_data2$County != "OTHER (COMBINED) COUNTIES",]

soy_data2 <- soy_data2 %>% group_by(region) %>% summarise(value = mean(value))

#Vectors with data for 2011 and 2012
fields$Year <- as.factor(fields$Year)
fields.2011 <- fields[fields$Year==2011,]
fields.2012 <- fields[fields$Year==2012,]

#ggcounty
soy_data2$brk <- cut(soy_data2$value, breaks=c(0, 1, 50000, 100000,200000,300000,700000), 
                     labels=c("Not Estimated","0-50K","50-100K","100-200K","200-300K","300-700K"))

us_state <- map_data(map = "state")
us <- ggcounty::ggcounty.us()

##US map
gg <- us$gg
gg <- gg + geom_map(data=soy_data2, map=us$map, aes(map_id=region, fill=brk), color="gray10", size=0.125) + 
  coord_quickmap(xlim = c(-106,-80), ylim=c(22,55)) + 
  scale_fill_brewer(palette = "Greens", name="Soybean planted acres") + 
  geom_path(data = us_state, colour="gray15", aes(x=long, y=lat, group=group))

gg2 <- gg + geom_point(data = fields, aes(x=Long, y=Lat, shape=Year), 
                       stat = "identity", size=3, fill = "gray70",
                       position=position_jitter(h=0.1)) + 
  scale_shape_manual(values = c(21,24)) +
  theme(legend.key = element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(c(-20,-10,-10,-10), "mm"))

#Obtain canada map
can2 <- getData('GADM', country="CAN", level=2) # counties
ca.province <- fortify(can2[can2$NAME_1 %in% "Ontario",], region = "NAME_2")

#Read data for locations
Ontario <- read.csv2("data/clean/CA_GIS_data.txt", 
                     header = TRUE,sep = "\t", dec = ".")
Ontario$Lat <- as.numeric(Ontario$Lat)
Ontario$Long <- as.numeric(Ontario$Long)

#Read data for soy pdx
On_soy <- read.csv2("data/clean/Ontario_soybean.txt", 
                    sep = "\t", header = TRUE)

On_soy$brk <- cut(On_soy$Acres, 
                  breaks=c(0, 1, 50000, 100000,200000,300000,700000), 
                  labels=c("Not Estimated", "0-50K","50-100K","100-200K","200-300K","300-700K"))

#plot for data
gg3 <- ggplot(data=ca.province) + geom_path(colour="#7f7f7f", 
                                            aes(x=long, y=lat, group=group)) +
  coord_quickmap(xlim=c(-95.15, -74.34), ylim=c(41.67, 56.86)) +
  geom_map(data=On_soy, map=ca.province, aes(map_id=County, fill=brk), 
           color="gray10", size=0.125) +
  scale_fill_brewer(palette = "Greens", name="Soybean planted acres") + 
  guides(fill="none") +
  geom_point(data = Ontario, aes(x=Long, y=Lat, shape=as.factor(Year)), 
             stat = "identity", size=3, fill = "gray70", position = "jitter") +
  scale_shape_manual(values = c(24)) +
  guides(shape = FALSE) + 
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        legend.key = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

