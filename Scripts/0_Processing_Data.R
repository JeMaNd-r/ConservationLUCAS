#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Comparison soil parameters in protected           #
#     vs. non-protected areas, LUCAS data              #
#                                                      #
#               author: Romy Zeiss                     #
#                                                      #
#     0. Preprocessing and defining protected sites    #
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

# set working directories
work.wd <- "~/LUCAS"
data.wd <- "I:/eie/==PERSONAL/RZ SoilBON/RZ LUCAS"
figu.wd <- "H:/Figures"
setwd(work.wd); getwd()

# load packages
library(ggplot2)
library(dplyr)
library(MBESS)   # to calculate CI of Cohens d effect size
library(reshape2)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load LUCAS data ####
setwd(data.wd); getwd()
lucas <- read.csv("Full_LUCAS2018-iDiv_dataset_02072021.csv", sep=",")
str(lucas)

data.old <- read.csv("Full_LUCAS2018-iDiv_dataset_01032021.csv", sep=",")
str(data.old)

# quick look at coordinates
par(mfrow=c(1,2))
plot(data.old$X_LAEA, data.old$Y_LAEA)
plot(lucas$TH_LONG, lucas$TH_LAT)
par(mfrow=c(1,1))

# replace spelling mistake in Country column
unique(data$Country)
data[data$Country=="germany","Country"] <- "GERMANY"
unique(data$Country)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Handle coordinates ####
library(raster)
library(rgdal)  #to load shapefile
library(rgeos)  #to crop shapefile
#library(sf)     #to manage large shapefile

# select rows with coordinates and land cover type
lucas.coor <- lucas[,c(1,8:12,57)]
colnames(lucas.coor) <- replace(colnames(lucas.coor), colnames(lucas.coor)=="LC_5", "LC")
head(lucas.coor)

# # Transform coordinates
# d <- data.frame(lon=lucas.coor$TH_LONG, lat=lucas.coor$TH_LAT)
# coordinates(d) <- c("lon", "lat")
# proj4string(d) <- CRS("+proj=longlat +datum=WGS84") # ETRS89-extended / LAEA Europe
# CRS.new <- CRS("+init=epsg:4326")
# d.new <- spTransform(d, CRS.new)
# unclass(d.new)
# 
# # Plot the results
# par(mfrow=c(1,3))
# plot.default(lucas.coor$TH_LONG,lucas.coor$TH_LAT, main="Raw data", cex.axis=.95)
# plot(d, axes=TRUE, main="Original lat-lon", cex.axis=.95)
# plot(d.new, axes=TRUE, main="Projected", cex.axis=.95)
# par(mfrow=c(1,1))
# 
# # add transformed coordinated to data
# lucas.coor$lon <- d.new$lon
# lucas.coor$lat <- d.new$lat
# head(lucas.coor)

rm(d)
#rm(d.new)

# create SpatialPointsDataFrame
d <- data.frame(Longitude = lucas.coor$TH_LONG,
                Latitude = lucas.coor$TH_LAT,
                names = lucas.coor$LUCAS_ID)

coordinates(d) <- ~ Longitude + Latitude

# save shapefile with LUCAS coordinates
writeOGR(d, layer="LUCAS_all_points", data.wd, driver="ESRI Shapefile")

