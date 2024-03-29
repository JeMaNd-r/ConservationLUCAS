#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Comparison soil parameters in protected           #
#     vs. non-protected areas, LUCAS data              #
#                                                      #
#               author: Romy Zeiss                     #
#                                                      #
#     1. Pairing and comparing sample sites            #
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

# set working directories
work.wd <- "I:/eie/==PERSONAL/RZ_SoilBON/ConservationLUCAS"
data.wd <- "I:/eie/==PERSONAL/RZ_SoilBON/ConservationLUCAS/Data/Data_v3"
figu.wd <- "I:/eie/==PERSONAL/RZ_SoilBON/ConservationLUCAS/Figures/Figures_v3"
setwd(work.wd); getwd()

# load packages
library(ggplot2)
library(tidyverse)
library(reshape2)
library(psych)   # to calculate CI of Cohens d effect size

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## load shape file(s) for protected areas ####
# based on World Database on Protected Areas (WDPA) (protectedplanet.net)

#... done in ArcMap

# crop Europe 
#template <- extent(-15,35, 30,75) #longitude left - right, latitude min - max

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load LUCAS data ####
setwd(data.wd); getwd()
lucas <- read.csv("Full_LUCAS2018-iDiv_dataset_02072021.csv", sep=",")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Add location in or outside of protected area ####
# after checking for location with ArcMap
lucas

# load table created in ArcMap with points in PA only
lucas.pa <- read.table("LUCAS_points_withPA.txt", sep=",", header=T)
colnames(lucas.pa) <- c("FID", "LUCAS_ID", "PA")

# add column with information about protected and non-protected sites
lucas <- full_join(lucas, lucas.pa[,2:3], by="LUCAS_ID") #....
#lucas[is.na(lucas$PA), "PA"] <- 0  #replace NA by 0
summary(as.factor(lucas$PA))

# number of observations (raw)
nrow(lucas); nrow(lucas[lucas$PA,])  #893 with 102 PAs

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## pair one PA with most similar nonPA sample point ####
# based on physical properties and distance

# define functions to be compared
fns = c("Basal_respiration","Cmic","Xylosidase","Cellulase","beta.glucosidase",
        "N.actylglucosaminidase","Acid.phosphatase","MWD","WSA")

## define variables to be compared between PA and nonPA, & threshold
mahal.vars <- c("TH_LAT", "TH_LONG", "Elevation", "AnnualPrecip", "AnnualTemp", 
                "MonthlyPrecipSum","MonthlyMeanTemp", "Organic_carbon",
                "Clay_content", "Silt_content",#"sand_yearCombine", #remove one because colinear
                "pH_H2O")

# keep complete cases only
lucas <- lucas[complete.cases(lucas[,c(mahal.vars, fns)]),]
nrow(lucas); nrow(lucas[lucas$PA==1,])  # n = 876 including 101 PAs

mahal.thres <- qchisq(.975, df=length(mahal.vars)) #21.92005

# scale variables for mahalanobis distance
for(i in mahal.vars){
  lucas[,paste0(i,".z")] <- as.numeric(scale(lucas[,i]))
}

mahal.vars <- paste0(mahal.vars, ".z")

# define each land cover type
lc.names <- c("Cropland", "Woodland", "Grassland", "Other")

# Note: There were no nonPA sites at all with mahalanobis distance below threshold
# for three protected sites with LUCAS_ID. They have to been removed.
lucas <- lucas[lucas$LUCAS_ID!=31783686,]
lucas <- lucas[lucas$LUCAS_ID!=28302282,]
lucas <- lucas[lucas$LUCAS_ID!=29542352,]

# split data
nonpa <- lucas[lucas$PA==0,c("LUCAS_ID", "LC_5", "PA", mahal.vars, fns)]
pa <- lucas[lucas$PA==1,c("LUCAS_ID", "LC_5", "PA", mahal.vars, fns)]

# randomization
p.list.total <- vector("list", length = 1000)
effect.size.d <- vector("list", length=1000)
pa.pairs <- data.frame(LUCAS_ID=NULL, nonPA=NULL, mahal.min=NULL, LC=NULL, run=NULL)
missing.pa <- data.frame(run=NA, pa.site=NA)[0,]
times <- 0; times.with.error <- 0; set.seed(1) 
repeat {
  times <- times+1; if(times > 1000) {break} #stop loop if reached 1000 trails
  times.with.error <- times.with.error + 1
  
  check.error <- try({  # for doing run again if error occurs (e.g. in pairing)
  
    # add columns to store pairing temporary
    pa[,c("nonPA", "mahal.min")] <- NA  
    
    for(i in 1:length(lc.names)){
      data_PA <- pa[pa$LC_5==lc.names[i],] 
      data_nonPA  <- nonpa[nonpa$LC_5==lc.names[i],]
      
      # select environmental data only as matrix, remove rows with NAs
      data_PA <- data_PA[complete.cases(data_PA[,mahal.vars]),c("LUCAS_ID", mahal.vars)]
      data_nonPA <- data_nonPA[complete.cases(data_nonPA[,mahal.vars]),c("LUCAS_ID", mahal.vars)]
      
      sigma <- cov(data_nonPA[,mahal.vars]) # covariance/correlation between variables
      
      # random order of PA sites
      data_PA <- data_PA[order(sample(data_PA$LUCAS_ID)),]
      
      # add empty columns
      data_nonPA[,as.character(data_PA$LUCAS_ID)] <- NA
      
      # calculate Mahalanobis distance
      for(j in 1:nrow(data_PA)){
        mu = as.numeric(data_PA[j,mahal.vars])
        data_nonPA[,as.character(data_PA[j,"LUCAS_ID"])] <- 
          mahalanobis(data_nonPA[,mahal.vars], mu, sigma, tol=1e-30)
        #print(j)
      }
      
      # add column to PA data with respective nonPA sites 
      # based on minimal mahalanobis distance
      data_PA[, c("nonPA", "mahal.min")]  <- NA
      #temp.col <- 0
        
      for(k in data_PA$LUCAS_ID){
        
        # if Mahal. distances compared to all relevant nonPA are above threshold...
        if(min(data_nonPA[,as.character(k)])>=mahal.thres){
          #... stop and re-do run
          missing.pa <- rbind(missing.pa,c(times, k))  # to know why we've stopped
          print(paste0("Not all PA sites paired. Check PA ", k, "."))
          
          pa[pa$LUCAS_ID==k, c("nonPA", "mahal.min")] <-
           c(NA, min(data_nonPA[,as.character(k)]))
          data_PA <- data_PA[data_PA$LUCAS_ID!=k,] # remove respective PA site
        }
        
        # select (max. 10) nonPA sites with mahalanobis distance below threshold
        nonPA.pair <- data_nonPA[data_nonPA[,as.character(k)]<=mahal.thres,c("LUCAS_ID", as.character(k))]
        nonPA.pair[nrow(nonPA.pair)+1:nrow(nonPA.pair)+10,] <- NA #add empty rows for when there are less than 10 sites
        nonPA.pair <- nonPA.pair[order(nonPA.pair[,2]),]
        nonPA.pair <- nonPA.pair[c(1:10),1]
      
        # sample one of the top 10 nonPA site that isn't NA 
        data_PA[data_PA$LUCAS_ID==k,"nonPA"] <- sample(as.character(nonPA.pair[!is.na(nonPA.pair)]),1)
        
        # add value of distance
        data_PA[data_PA$LUCAS_ID==k,"mahal.min"] <- unique(data_nonPA[data_nonPA$LUCAS_ID==data_PA[data_PA$LUCAS_ID==k,"nonPA"],as.character(k)])
        data_nonPA <- data_nonPA[data_nonPA$LUCAS_ID!=data_PA[data_PA$LUCAS_ID==k,"nonPA"],]
      }
      
      
      # add to result table
      pa[pa$LUCAS_ID %in% data_PA$LUCAS_ID, c("nonPA", "mahal.min")] <-
        data_PA[order(as.numeric(rownames(data_PA))),c("nonPA", "mahal.min")]
      
      # add to result table to analyse all at once below
      pa.pairs <- rbind(pa.pairs, cbind(data_PA[,c("LUCAS_ID", "nonPA", "mahal.min")],
                                        lc.names[i], times.with.error, times))
    }
    
#})}  # stop here if only calculating p-value for all runs together...
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Compare difference between PA and nonPA sites ####
  
  # merge pairs of PA and nonPA in one table ####
  nonpa.subset <- nonpa[unique(nonpa$LUCAS_ID) %in% pa$nonPA,]
  lucas.paired <- full_join(pa[order(pa$nonPA),], nonpa.subset[order(nonpa.subset$LUCAS_ID),])
  #head(lucas.paired)
  
  # Perform t tests
  p.list <- list()
  effect.size.d[[times]] <- list()
  for(l in lc.names){
    p.list[[l]] <- vector("list", length(fns))
    names(p.list[[l]]) <- fns
    
    # subset data based on land cover type
    temp.PA <- lucas.paired[lucas.paired$LC_5==l & lucas.paired$PA==1,]
    temp.nonPA <- lucas.paired[lucas.paired$LC_5==l & lucas.paired$PA==0,]
    temp.nonPA <- temp.nonPA[order(match(temp.nonPA$LUCAS_ID, temp.PA$nonPA)), ]
    
    temp.cohens <- psych::cohen.d(rbind(temp.PA, temp.nonPA)[,c("PA",fns)], "PA")
    effect.size.d[[times]][[l]] <- cbind(lc=l, data.frame(temp.cohens$cohen.d), run=times)
    
    for(no.fns in 1:(length(fns))){
      # unpaired t test
      p.list[[l]][[no.fns]] <- t.test(temp.PA[,fns[no.fns]],temp.nonPA[,fns[no.fns]])
    }
    print(summary(temp.PA[,"nonPA"] == temp.nonPA[temp.nonPA$LUCAS_ID %in% temp.PA[,"nonPA"],"LUCAS_ID"])) # make sure that the sites are pairing properly
    
  }
  effect.size.d[[times]] <- do.call(rbind, effect.size.d[[times]])
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Make table of all p values for all functions and LC types
  p.table = as.data.frame(matrix(ncol=length(fns),nrow=5))
  colnames(p.table) = fns
  rownames(p.table) = c("Other","Cropland","Grassland", "Woodland","total")
  
  for(n in 1:ncol(p.table)){
    p.table[1,n] = p.list[["Cropland"]][[n]]["p.value"]
    p.table[2,n] = p.list[["Grassland"]][[n]]["p.value"]
    p.table[3,n] = p.list[["Woodland"]][[n]]["p.value"]
    p.table[4,n] = p.list[["Other"]][[n]]["p.value"]
    #... total missing
  }
  
  # add to overall result list
  p.list.total[[times]] <- p.table  #p-values from Chi-squared test
  })
  
  # do run again if there is an error (i.e. if no nonPA site with distance lower than threshold)
  if(is(check.error,"try-error")) {times <- times-1} else {print(check.error)}

} #end of randomization

save(effect.size.d,  file="d_1000_trails.RData")
load(file="d_1000_trails.RData")

# show total count of unpaired (and removed) PAs and compare with number of paired sites
table(missing.pa[,2])  #should be 0
table(pa.pairs$LUCAS_ID) # all counts should be number of runs (i.e. times)

# check for runs that failed (i.e. count < number of PA sites), and remove the respective pairs
# note: other result objects are not effected as they were overwriten
nrow(pa.pairs)  # more than 97 * 1000
pa.pairs <- pa.pairs %>% add_count(times.with.error) %>% 
  filter(n==length(lucas[lucas$PA==1,"LUCAS_ID"]))
nrow(pa.pairs)  # exactly 97 protected sites * 1000

# look what non-protected sites have (not) been paired to any PA
hist(table(pa.pairs$nonPA))  # frequency distribution of the use of sites from all runs
length(setdiff(lucas[lucas$PA==0,]$LUCAS_ID, pa.pairs$nonPA))  # sites never used: 269

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Save total list with p tables
save(p.list.total, file="p_1000_trails.RData")
write.csv(p.list.total, file="p_1000_trails.csv")

# 2. do a p-value test on all runs at once (i.e. collapsing all 1000 
# pairs into one table and perform a pair-wise comparison) per LC-type
pa.pairs <- merge(pa.pairs, lucas[,c("LUCAS_ID", fns)], by="LUCAS_ID")
pa.pairs <- merge(pa.pairs, lucas[,c("LUCAS_ID", fns)], by.x="nonPA", by.y="LUCAS_ID")
colnames(pa.pairs)[4] <- "LC"

save(pa.pairs, file="Pairs_paNonpa_1000trails.RData")
write.csv(pa.pairs, file="Pairs_paNonpa_1000trails.csv", row.names=F)
load("Pairs_paNonpa_1000trails.RData") #pa.pairs

# - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Extract Cohen's d effect size for each run ####
d.table <- do.call(rbind, effect.size.d)
d.table$fns <- rep(fns)
head(d.table)

save(d.table, file="Effect_size_d_table.RData")
#write.csv(d.table, file="Effect_size_d_table.csv", row.names = F)

#load(file="Effect_size_d_table.RData") #d.table

## Summarize per land use and function (i.e. get mean and confidence intervals)
d.table.mean <- d.table %>% group_by(lc, fns) %>% summarize_all(mean, na.rm=T)
d.table.mean
write.csv(d.table.mean, file="Effect_size_d_mean.csv", row.names = F)

# Merge mean, max and min into one data frame
d.table.summary <- d.table.mean %>% rename("mean.D"=effect, "CI.upper"=upper, "CI.lower"=lower)
head(d.table.summary)

write.csv(d.table.summary, file="Effect_size_d_summary.csv", row.names = F)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# combine different results in one table
setwd(data.wd); load("p_1000_trails.RData") #p.list.total
head(p.list.total)

pvals.agri <- data.frame(run=1:1000); pvals.agri[,fns] <- NA
pvals.gras <- data.frame(run=1:1000); pvals.gras[,fns] <- NA
pvals.fore <- data.frame(run=1:1000); pvals.fore[,fns] <- NA
pvals.othe <- data.frame(run=1:1000); pvals.othe[,fns] <- NA

for(i in 1:length(p.list.total)){
  try({
  pvals.agri[i,-1] <- p.list.total[[i]]["Cropland",]
  pvals.gras[i,-1] <- p.list.total[[i]]["Grassland",]
  pvals.fore[i,-1] <- p.list.total[[i]]["Woodland",]
  pvals.othe[i,-1] <- p.list.total[[i]]["Other",]
  })
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plot p-values of all runs ####
pvals.agri <- reshape2::melt(pvals.agri, id.vars=1, na.rm=T)
pvals.gras <- reshape2::melt(pvals.gras, id.vars=1, na.rm=T)
pvals.fore <- reshape2::melt(pvals.fore, id.vars=1, na.rm=T)
pvals.othe <- reshape2::melt(pvals.othe, id.vars=1, na.rm=T)

## Extract mean p-values
p.vals.mean <- data.frame(variable=NA, run=NA, value=NA, lc=NA)[0,] 
p.vals.mean <- rbind(p.vals.mean, cbind(aggregate(.~variable, pvals.agri, mean), lc="Cropland"))
p.vals.mean <- rbind(p.vals.mean, cbind(aggregate(.~variable, pvals.gras, mean), lc="Grassland"))
p.vals.mean <- rbind(p.vals.mean, cbind(aggregate(.~variable, pvals.fore, mean), lc="Woodland"))
p.vals.mean <- rbind(p.vals.mean, cbind(aggregate(.~variable, pvals.othe, mean), lc="Other"))
p.vals.mean

write.csv(p.vals.mean, file="p_1000_trails_mean.csv", row.names = F)

## Extract maximum/upper quartile of p-values
p.vals.max <- data.frame(variable=NA, run=NA, value=NA, lc=NA)[0,] 
p.vals.max <- rbind(p.vals.max, cbind(aggregate(.~variable, pvals.agri, 
                                                function(x) {quantile(x,0.75)}), lc="Cropland"))
p.vals.max <- rbind(p.vals.max, cbind(aggregate(.~variable, pvals.gras, 
                                                function(x) {quantile(x,0.75)}), lc="Grassland"))
p.vals.max <- rbind(p.vals.max, cbind(aggregate(.~variable, pvals.fore, 
                                                function(x) {quantile(x,0.75)}), lc="Woodland"))
p.vals.max <- rbind(p.vals.max, cbind(aggregate(.~variable, pvals.othe, 
                                                function(x) {quantile(x,0.75)}), lc="Other"))
p.vals.max

# minimum/ lower quartile values
p.vals.min <- data.frame(variable=NA, run=NA, value=NA, lc=NA)[0,] 
p.vals.min <- rbind(p.vals.min, cbind(aggregate(.~variable, pvals.agri, 
                                                function(x) {quantile(x,0.25)}), lc="Cropland"))
p.vals.min <- rbind(p.vals.min, cbind(aggregate(.~variable, pvals.gras, 
                                                function(x) {quantile(x,0.25)}), lc="Grassland"))
p.vals.min <- rbind(p.vals.min, cbind(aggregate(.~variable, pvals.fore, 
                                                function(x) {quantile(x,0.25)}), lc="Woodland"))
p.vals.min <- rbind(p.vals.min, cbind(aggregate(.~variable, pvals.othe, 
                                                function(x) {quantile(x,0.25)}), lc="Other"))
p.vals.min

# Merge mean, max and min into one data frame
colnames(p.vals.mean)[3] <- "mean.P"
colnames(p.vals.max)[3] <- "quartile75"
colnames(p.vals.min)[3] <- "quartile25"

p.vals.summary <- merge(p.vals.mean[,-2], p.vals.max[,-2], 
                        by=c("variable", "lc"))
p.vals.summary <- merge(p.vals.summary, p.vals.min[,-2], by=c("variable", "lc"))
head(p.vals.summary)

write.csv(p.vals.summary, file="p_1000_trails_summary.csv", row.names = F)

# merge with d values
statistics <- merge(p.vals.summary, d.table.mean, 
                    by.x=c("variable", "lc"), by.y=c("fns", "lc"))
head(statistics)
write.csv(statistics, file="1000_trails_statistics_p+d.csv", row.names = F)

## violin plot per land use type
# agri
aplot <- ggplot(data=pvals.agri, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.05, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  ggtitle("Cropland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))

# gras
gplot <- ggplot(data=pvals.gras, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.05, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Grassland") +
  stat_summary(fun = "mean",geom = "point",color = "black")+
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.text.y = element_blank(),  axis.text.x = element_text(size=20))

# Woodland
fplot<- ggplot(data=pvals.fore, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.05, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Woodland") +
  theme_bw() + # use a white background
  stat_summary(fun = "mean",geom = "point",color = "black")+
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.text.y = element_blank(),  axis.text.x = element_text(size=20))

# others
oplot<- ggplot(data=pvals.othe, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.05, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Others") +
  theme_bw() + # use a white background
  stat_summary(fun = "mean",geom = "point",color = "black")+
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.text.y = element_blank(),  axis.text.x = element_text(size=20))

library(ggpubr)
setwd(figu.wd); pdf(file="Pvals_violin_1000trails_ttest.pdf", width=15)
ggarrange(aplot, gplot, fplot, oplot, ncol = 4, nrow = 1, align = "none",
          widths=c(2,1,1,1))
dev.off()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Effect size as violin plot ####
setwd(data.wd)
load(file="Effect_size_d_table.RData")

dvals.all <- reshape2::melt(d.table, id.vars=c("lc", "run", "fns"))
dvals.agri <- reshape2::melt(d.table[d.table$lc=="Cropland",], id.vars = c("lc", "run", "fns"))
dvals.gras <- reshape2::melt(d.table[d.table$lc=="Grassland",], id.vars = c("lc", "run", "fns"))
dvals.fore <- reshape2::melt(d.table[d.table$lc=="Woodland",], id.vars = c("lc", "run", "fns"))
dvals.othe <- reshape2::melt(d.table[d.table$lc=="Other",], id.vars = c("lc", "run", "fns"))

## violin plot per land use type
setwd(figu.wd); pdf(file="Dvals_violin_1000trails_in1.pdf", height=15)
ggplot(data=dvals.all[dvals.all$lc!="Other" & dvals.all$variable=="effect",], 
       aes(x=value, y=fns, fill=lc, width=0.7))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  #facet_grid(lc~.) +
  theme(#legend.position = "none", 
    axis.title.y =element_blank(), 
        axis.text = element_text(size=10))
dev.off()

## or as separate plots merged at the end
# agri
aplot <- ggplot(data=dvals.agri, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Cropland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))+
  scale_x_continuous(limits = c(-1.5,0.8))

# gras
gplot <- ggplot(data=dvals.gras, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Grassland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))+
  scale_x_continuous(limits = c(-1.5,0.8))

# Woodland
fplot<- ggplot(data=dvals.fore, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Woodland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))+
  scale_x_continuous(limits = c(-1.5,0.8))

# others
oplot<- ggplot(data=dvals.othe, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Others") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))+
  scale_x_continuous(limits = c(-1.5,0.8))

library(ggpubr)
setwd(figu.wd); pdf(file="Dvals_violin_1000trails.pdf", height=15)
ggarrange(aplot, gplot, fplot, ncol = 1, nrow = 3, align = "none")
dev.off()

# ## with confidence intervals ####
setwd(data.wd)
d.table.summary <- read.csv("Effect_size_d_summary.csv")
head(d.table.summary)

str(d.table.summary)
# sort labels
d.table.summary$label <- factor(d.table.summary$fns, 
                                c("Basal_respiration", "Cmic", "Acid.phosphatase",
                                  "beta.glucosidase", "Cellulase", "N.actylglucosaminidase",
                                  "Xylosidase", "MWD", "WSA"))
name.labels <- c("Basal respiration", "Microbial biomass", "Acid-Phosphatase",
  "Beta-Glucosidase", "Cellulase", "N-Actylglucosaminidase",
  "Xylosidase", "Mean weight diameter", "Water stable aggregates")

aplot <- ggplot(data=d.table.summary[d.table.summary$lc=="Cropland",], 
                aes(y=mean.D, x=label, ymin=CI.lower,ymax=CI.upper))+
  geom_pointrange() +
  geom_hline(yintercept=0, linetype="dashed")+
  ggtitle("Cropland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),plot.margin = unit(c(1, 1, 1, 2), "cm"),
        axis.text.y = element_text(size=20),  axis.text.x = element_blank())+
  scale_y_continuous(limits = c(-1,1.1))

gplot <- ggplot(data=d.table.summary[d.table.summary$lc=="Grassland",], 
                aes(y=mean.D, x=label, ymin=CI.lower,ymax=CI.upper))+
  geom_pointrange() +
  ggtitle("Grassland") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),plot.margin = unit(c(1, 1, 1, 2), "cm"),
        axis.text.y = element_text(size=20),  axis.text.x = element_blank())+
  scale_y_continuous(limits = c(-1,1.1))

fplot <- ggplot(data=d.table.summary[d.table.summary$lc=="Woodland",], 
                aes(y=mean.D, x=label, ymin=CI.lower,ymax=CI.upper))+
  geom_pointrange() +
  ggtitle("Woodland") +
  geom_hline(yintercept=0, linetype="dashed")+
  scale_x_discrete(labels= name.labels) +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(), plot.margin = unit(c(1, 1, 1, 2), "cm"),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20, angle=45, hjust=1))+
  scale_y_continuous(limits = c(-1,1.1))

library(ggpubr)
setwd(figu.wd); pdf(file="Dvals_pointrange_1000trails.pdf", width=5, height=11)
ggarrange(aplot, gplot, fplot, ncol = 1, nrow = 3, align = "none", heights=c(1,1,1.7))
dev.off()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Boxplots ####

# Standardize
std.data <- lucas[lucas$LUCAS_ID %in% unique(c(pa.pairs$LUCAS_ID, pa.pairs$nonPA)),]
std.data[,fns] <- scale(std.data[,fns])

# Labels for functions
labeling <- c("BAS", "Cmic", "Xylo", "Cellu", "bGluco", "NAG", 
              "Phosp", "MWD", "WSA")

# Split up by land cover type
std.Woodland = std.data[which(std.data$LC_5 == "Woodland"),]
std.grass = std.data[which(std.data$LC_5 == "Grassland"),]
std.agri = std.data[which(std.data$LC_5 == "Cropland"),]
std.other = std.data[which(std.data$LC_5 == "Other"),]

# melt data to get one colume with both PA and nonPA variables
amelt <- melt(std.agri, id.vars=c("LUCAS_ID", "PA"), measure.vars = fns, na.rm=T)
gmelt <- melt(std.grass, id.vars=c("LUCAS_ID", "PA"), measure.vars = fns, na.rm=T)
fmelt <- melt(std.Woodland, id.vars=c("LUCAS_ID", "PA"), measure.vars = fns,na.rm=T)
omelt <- melt(std.other, id.vars=c("LUCAS_ID", "PA"), measure.vars = fns, na.rm=T)

## Plot
library(ggplot2)
library(ggpubr)

#Woodland
fplot <- ggplot(fmelt,aes(variable,value)) +
  geom_boxplot(aes(fill=as.factor(PA)))+#,fatten = NULL) + # activate to remove the median lines
  scale_fill_manual(values=c("orange","darkgreen"),name = NULL, labels = c("Non-Protected","Protected")) +
  #scale_y_continuous(limits=c(-3,13)) + 
  geom_hline(yintercept = 0, lty="dashed") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(),text = element_text(size=20),
        legend.position = c(0.9,0.9), legend.text = element_text(size=20), axis.text.y = element_text(size=15)) +
  labs(x=NULL,y="Woodland") #+
  #geom_text(x=2,y=12.5,label="*",size=9,color="red3")

#Cropland
aplot <- ggplot(amelt,aes(variable,value)) +
  geom_boxplot(aes(fill=as.factor(PA)))+#,fatten = NULL) + # activate to remove the median lines
  scale_fill_manual(values=c("orange","darkgreen"),name = NULL, labels = c("Non-Protected","Protected")) +
  #scale_y_continuous(limits=c(-3,13)) + 
  geom_hline(yintercept = 0, lty="dashed") +
  theme_classic() +
  scale_x_discrete(labels=labeling) +
  theme(text = element_text(size=20),legend.position = "none", axis.text.y = element_text(size=15)) +
  labs(x=NULL,y="Cropland") #+
  #geom_text(x=1,y=2.25,label="*",size=9,color="red3")

#grassland
gplot <- ggplot(gmelt,aes(variable,value)) +
  geom_boxplot(aes(fill=as.factor(PA)))+#,fatten = NULL) + # activate to remove the median lines
  scale_fill_manual(values=c("orange","darkgreen"),name = NULL, labels = c("Non-Protected","Protected")) +
  #scale_y_continuous(limits=c(-3,13)) + 
  geom_hline(yintercept = 0, lty="dashed") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(),text = element_text(size=20),
        legend.position = "none", axis.text.y = element_text(size=15)) +
  labs(x=NULL,y="Grassland")

#setwd(figu.wd); pdf(file="LUCAS_Boxplot_paNonpa.pdf", width=12)
ggarrange(fplot, gplot, aplot, ncol = 1, nrow = 3, align = "v")
dev.off()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Boxplot mahalanobis distance ####
pa.pairs <- merge(pa.pairs, lucas[,c("LUCAS_ID", "Country")])

# get sample sizes
no.sample <- pa.pairs %>% count(Country, LC)
#write.csv(no.sample, file="LUCAS_Boxplot_mahal.distance_samplesize.csv", row.names=F)
#no.sample <- as.character(no.sample$n)

#setwd(figu.wd); pdf(file="LUCAS_Boxplot_mahal.distance_wide.pdf", width=15) 
ggplot(pa.pairs[pa.pairs$LC!="Other",],aes(Country,mahal.min, fill=LC)) +
  geom_boxplot()+#,fatten = NULL) + # activate to remove the median lines
  theme_classic() +
  labs(x="Country",y="Mahalanobis distance") +
  theme(axis.text.x=element_text(size=15, angle=45, hjust=1),text = element_text(size=20),  
        legend.position = c(0.6,0.8), axis.text.y = element_text(size=15), legend.title = element_blank())+
  scale_fill_manual(values=c("gold3", "forestgreen", "limegreen"))
dev.off()

# # manually check Cropland & Cmic because looks significant
# std.agri <- std.data[which(std.data$LC_5 == "Cropland") & !is.na(std.data$nonPA),]
# std.agri$pair <- rep(1:20, 2)[1:39]
# amelt = melt(std.agri,id.vars = c("PA", "pair"),measure.vars=fns)
# amelt = amelt[amelt$variable=="Cmic",]
# ggplot(amelt,aes(pair,value)) +
#   geom_point(aes(col=as.factor(PA))) + # activate to remove the median lines
#   scale_fill_manual(values=c("black",rgb(76,255,0,maxColorValue = 255)),name = NULL, labels = c("Non-Protected","Protected")) +
#   theme_classic() +
#   theme(text = element_text(size=13),legend.position = "bottom",axis.text.x=element_blank()) +
#   labs(x=NULL,y="Cropland") +
#   geom_line(aes(group = pair))
# 
# #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # check if paired sites are uniquely paired (i.e. no nonPA site doubled) 
# doubled <- as.data.frame(table(pa$nonPA.1)[table(pa$nonPA.1)>1]); doubled #those aren't
# sum(doubled$Freq) #previously: 30 sites are not uniquely paired, now 0
# 
# head(pa)
# 
# #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ## merge pairs of PA and nonPA in one table ####
# nonpa.subset <- nonpa[nonpa$LUCAS_ID %in% pa$nonPA.1,]
# lucas.paired <- full_join(pa, nonpa.subset)
# head(lucas.paired)
# 
# ## save pairs of PA and nonPA sites
# setwd(data.wd)
# #write.csv(lucas.paired, file="LUCAS_PA_pairs.csv")
# lucas.paired <- read.csv("LUCAS_PA_pairs.csv")
# 
# ## Make shapefile of selected sample points
# # Transform coordinates
# d <- data.frame(lon=lucas.paired$X_LAEA, lat=lucas.paired$Y_LAEA)
# coordinates(d) <- c("lon", "lat")
# proj4string(d) <- CRS("+init=epsg:3035") # ETRS89-extended / LAEA Europe
# CRS.new <- CRS("+init=epsg:4326")
# d.new <- spTransform(d, CRS.new)
# 
# # add transformed coordinated to data
# lucas.paired$lon <- d.new$lon
# lucas.paired$lat <- d.new$lat
# head(lucas.paired)
# 
# rm(d); rm(d.new)
# 
# # create SpatialPointsDataFrame
# d <- lucas.paired
# 
# coordinates(d) <- ~ lon + lat
# 
# # save shapefile with LUCAS coordinates
# #writeOGR(d, layer="LUCAS_PA_pairs_points", data.wd, driver="ESRI Shapefile", overwrite_layer = T)
# shapefile(d, "LUCAS_PA_pairs_points.shp", overwrite=T)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Make shapefile of sampling points (that have been used for pairing)
library(raster)

# select actual pairs only
lucas.clean <- lucas[lucas$LUCAS_ID %in% unique(c(pa.pairs$LUCAS_ID, pa.pairs$nonPA)),]
lucas.clean <- lucas.clean[lucas.clean$LC_5!="Other",]

# Transform coordinates
d <- cbind(lucas.clean, lon=lucas.clean$TH_LONG, lat=lucas.clean$TH_LAT)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:4326") 
# "+init=epsg:3035" for ETRS89-extended / LAEA Europe
# "+init=epsg:4326" for WGS84
# note: ignore possible warning message with "Discarded datum..."

# use different projection
#CRS.new <- CRS("+init=epsg:4326")  # just stick with European one
#d[,c("lon", "lat")] <- spTransform(d[,c("lon", "lat")], CRS.new)

head(d)

# save shapefile with LUCAS coordinates
setwd(data.wd)
#writeOGR(d, layer="LUCAS_PA_pairs_afterRun", data.wd, driver="ESRI Shapefile", overwrite_layer = T)
shapefile(d, "LUCAS_points_pairs_afterRun.shp", overwrite=T)

## number of observations per land-use type
lucas.clean %>% group_by(LC_5, PA) %>% count() 


