#- - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load LUCAS data set
##   Romy Zeiss 

# Set working directory
getwd()
setwd("~/Studium/Master BEE/2021 SoSe/Master thesis/Data")

# Load data into R
raw.dat <- read.csv("LUCAS/Full_LUCAS2018-iDiv_dataset_15072020_data.csv", sep=";")
str(raw.dat)

#- - - - - - - - - - - - - - - - - - - - - - - - - - 
## Select and rename colums ####
colnames(raw.dat)
dat <- raw.dat[,c(1:7, 9, 13, 21:29, 33:37, 74:75, 77:79)]
colnames(dat) <- replace(colnames(dat), colnames(dat)=="LC_4", "LC")

# select relevant soil property measure if it was a repeated measure in 2015
colnames(raw.dat) <- replace(colnames(raw.dat), colnames(raw.dat)=="CaCO2_2009", "CaCO3_2009")
soil.props <- sub("\\_.*", "", colnames(raw.dat[,c(44:55)]))
soil.props

for(i in soil.props){
  a <- paste(i, "_2009", sep="")
  b <- paste(i, "_2015", sep="")
  dat[,i] <- raw.dat[,b]
  for(j in 1:nrow(raw.dat)){
   if(is.na(raw.dat[j,b])) dat[j,i] <- raw.dat[j,a]
  }
}
dat$Electrical_conductivity <- raw.dat$Electrical_conductivity_2015
str(dat)

#- - - - - - - - - - - - - - - - - - - - - - - 
## Transform coordinates ####
x <- c(7.173500, 7.172540, 7.171636, 7.180180, 7.178070, 7.177229, 7.175240, 7.181409, 7.179299)
y <- c(45.86880, 45.86887, 45.86924, 45.87158, 45.87014, 45.86923, 45.86808, 45.87177, 45.87020)

## Define the coordinate systems
library(rgdal)
d <- data.frame(lon=dat$X_LAEA, lat=dat$Y_LAEA)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:3035") # ETRS89-extended / LAEA Europe
CRS.new <- CRS("+init=epsg:4326")
d.new <- spTransform(d, CRS.new)

# Plot the results.
par(mfrow=c(1,3))
plot.default(dat$X_LAEA,dat$Y_LAEA, main="Raw data", cex.axis=.95)
plot(d, axes=TRUE, main="Original lat-lon", cex.axis=.95)
plot(d.new, axes=TRUE, main="Projected", cex.axis=.95)
unclass(d.new)
dat$lon <- d.new$lon
dat$lat <- d.new$lat

#- - - - - - - - - - - - - - - - - - - - - - - - - - 
## Save selected data ####
write.csv(dat, "LUCAS/test_SDM_data.csv")

#- - - - - - - - - - - - - - - - - - - - - - - - - - 
## Load untransformed data from Linnea ####
SDM.dat <- read.csv("LUCAS/untransformed_SEM_data.csv")
str(SDM.dat)
