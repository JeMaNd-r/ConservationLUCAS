#- - - - - - - - - - - - - - - - - - - - - - - - - - 
# Analyse enzyme data
##   Romy Zeiss 

# Set working directory
getwd()
setwd("~/Studium/Master BEE/2021 SoSe/Master thesis/Data")

## Load helpful functions
source("../column_description.R") #to get overview of datasets
source("../diag_sem.R")           #to make easy & nice(r) looking SEM path plots

# Load data into R
raw.dat <- read.csv("LUCAS/Full_LUCAS2018-iDiv_dataset_15072020_data.csv", sep=";")
str(raw.dat)

dat <- read.csv("LUCAS/test_SDM_data.csv")
str(dat)

# use function column_description that gives summary for dataset based on columns
col.inf <- col_descr(dat)
View(col.inf)

#- - - - - - - - - - - - - - - - - - - - - - - - - - 
## GLMMs ####

#- enzymes ~ factors
#- think about models first!
#- check for methods with raw.data (water content, device, monthly weather)




#- - - - - - - - - - - - - - - - - - - - - - - - - - 
## SDMs ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - 
library(raster)
library(rgdal)

## Extract bioclim rasters ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - 
## Convert .tif files into .asc files
for(i in c("bio_1", "bio_4", "bio_12", "bio_15","elev")){
  f <- paste0("./BioClim/", "wc2.1_10m_", i,".tif")
  r <- raster(f)
  ra <- aggregate(r, fact=2)  ## By default aggregates using mean, but see fun=
  writeRaster(ra, paste0("./BioClim/", i,".asc"), format="ascii")
}

## Load .asc files with environmental data
bio1 <- raster("BioClim/bio_1.asc")  # BIO1 = Annual Mean Temperature
bio4 <- raster("BioClim/bio_4.asc")  # BIO4 = Temperature Seasonality (standard deviation ×100)
bio12 <- raster("BioClim/bio_12.asc")  # BIO12 = Annual Precipitation
bio15 <- raster("BioClim/bio_15.asc")  # BIO15 = Precipitation Seasonality (Coefficient of Variation)
bio.e <- raster("BioClim/elev.asc") # elevation

## Crop Europe
template <- extent(-39,90, 30,80)
#shapefile <- readOGR("./BioClim", layer = "ne_50m_admin_0_countries_lakes" ) 

# crop env layers to Europe 
for(i in c(bio1, bio4, bio12, bio15, bio.e)){
  n <- names(i)
  temp.crop <- crop(i, template, snap="out")
  crop <- setValues(temp.crop, NA)
  template.r <- rasterize(template, crop)
  m <- mask(x=temp.crop, mask=template.r)
  assign(paste0(n, ".masked"), m)
  plot(m, main=n)
  writeRaster(m, filename=(paste0("./BioClim/cropped/",n,".asc")))
}

### Make List-Object with new cropped layers
lst <- list(bio_1.masked,bio_4.masked,bio_12.masked,bio_15.masked)
bio <- stack(lst)

### Plot layers
mar <- c(1,1,1,1)
par(mfrow=c(2,2))
plot(bio[[1]], main = "bio1")
plot(bio[[2]], main = "bio4")
plot(bio[[3]], main = "bio12")
plot(bio[[4]], main = "bio15") #!!!
names(bio)

# - - - - - - - - - - - - - - - - - - - - - - #
## Build SDMs ####
library(sdm)
library(raster)
library(rgdal)
library(sp)
library(usdm)

## check correlations (not done)
# sp <- read.csv('occurrence_records_clean_10.csv')
# sp <- sp[,c('decimalLongitude','decimalLatitude')]
# spx <- extract(bio, sp)
# spx <- data.frame(spx)
# v <- vifstep(spx)
# v
# bio <- exclude(bio, v)
# bio

setwd('SDMs/')
wd = getwd()
mar <- c(1,1,1,1)
enzymes <- c("Xylosidase", "Cellulase", "N.actylglucosaminidase","beta.glucosidase","Acid.phosphatase")

#for(i in enzymes){
  enzyme <- dat[,i]
  sp <- cbind(dat[,c("lon", "lat")], enzyme)
  head(sp)
  coordinates(sp) <- ~ lon + lat # have to be transformed
  d <- sdmData(enzyme~., train=sp, predictors= bio, bg=list(n=1000))
  methods <- c('maxlike','bioclim', 'bioclim.dismo', 'glm') # which methods required? - e.g. 'bioclim','bioclim.dismo','gam','glm','maxlike','rpart'
  m <- sdm(enzyme ~ . , d, methods=methods, replication='sub', test.p=30, n=1) # data partitioning/ test percentage?
  #ensemble_model <- ensemble(m,bio, setting=list(id=1:2, method='weighted', stat='AUC', opt=2)) # how to combine models?
  
  writeRaster(ensemble_model, filename=paste0(wd,'/SDM_',i,'.asc'), format="ascii", overwrite=TRUE)
  plot(m, main = paste0("SDM: ", i), col=colorRampPalette(c('navy','lightyellow','orange','red'))(50)) # use zlim=c(0,1) if you want to scale between 0 and 1
  lst1 <- list.files(path=getwd(),pattern='grd$',full.names = T) 
  lst2 <- list.files(path=getwd(),pattern='gri$',full.names = T) 
  file.remove(lst1)
  file.remove(lst2)
}

plot(sp, col=sp$enzyme)

