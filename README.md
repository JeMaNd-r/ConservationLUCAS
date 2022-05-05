# Text mining in German management plans

## Table of contents

* [General info](#1-general-info)
* [Technologies](#2-technologies)
* [Setup](#3-setup)
* [References](#references)


## 1. General info

Belowground biodiversity has been linked to many ecosystem functions and services related again to human health and wellbeing, but faces threats similar to the aboveground system. Acting to halt biodiversity loss requires transformative change and effective conservation of all ecosystems. Indeed, soil biodiversity has been neglected in several global biodiversity assessments and conservation actions. To better understand how and why soil biodiversity and related ecosystem functions have lesser priority for nature conservation, it is important to get an overview of how society, and particularly policy-makers, have addressed it in the past and present. 

With this analysis, we aimed to review the current status of soil protection by investigating the effect of current conservation efforts in Europe as represented by nature conservation areas.

## 2. Technologies

Analysis are done in:
* R version 3.6.3
* RStudio version 1.2.5033
* ArcMap version 10.7.1

Attached R packages:

* dplyr_1.0.2  
* ggplot2_3.3.3      
* ggpubr_0.4.0  
* MBESS_4.8.0  
* psych_2.0.12   
* raster_3.4-5
* reshape2_1.4.4    
* sp_1.4-5  

## 3. Setup

If R and ArcGIS are already installed, you only have to download this repository.

### 3.1 Installation guide

After installing and loading all packages necessary, you can start with the actual analysis (see [Instructions](#33-instructions-for-use)). 
To run the analysis, you have to open the script `0_define_PA_nonPA.R`.

Typical install time on a "normal" desktop computer: ~5min


### 3.2 Demo

To see if the installation works, it is recommended to run the analysis on a subset of randomized runs (i.e., 10 or 100 instead of 1000) first. This can be done by replacing `if(times > 1000)` with `if(times > 10)` or something similar.

To fix error messages among others, you can try to run the code within loops line by line with i (or j or k) defined by hand according to the necessary type of object.

The output should be saved according to the R script. The expected run time for demo on a "normal" desktop computer depends on the subset of data used, but should not exceed 2h.

### 3.3 Instructions for use

The analysis includes the following steps, which are also descriped in the script `0_define_PA_nonPA.R`: 

#### Dataset

If you want to run the analysis right from the beginning, you need to have the data file called `Full_LUCAS2018-iDiv_dataset_02072021.csv`. This raw dataset from the Land use/cover area survey (LUCAS) soil project 2018 is available on [request](romy.zeiss@idiv.de). 

#### ArcGIS

Take the dataset `Full_LUCAS2018-iDiv_dataset_02072021.csv` and define a new column "PA" that has "0" as values. This column will have 0 for non-protected, and 1 for protected sites. To fill the column, you need to download the shapefile of the Natura 2000 nature conservation areas [here](https://www.eea.europa.eu/data-and-maps/data/natura-12/natura-2000-spatial-data/natura-2000-shapefile-1). If you added this shapefile to your ArcGIS project, use the `Select by location` function in ArcGIS to filter for protected areas. Change the value of the "PA" column of the selected sites to "1". Finally, save the information about the points and their classification (i.e., ID and PA column) to `LUCAS_points_withPA.txt`.

#### R

Run the script `0_define_PA_nonPA.R`.

* Add PA column to the LUCAS data set
* Define functions, environmental variables, land cover types, and protected areas.
* Create empty R objects to be filled during the loop.
* Repeat the pairing of sites & the statistical test for their differences 1000 times.
* Merge the output of all the 1000 runs into one table per estimated parameter (i.e., p-value, effect size, and confidence intervalls). 
* Plot the estimated parameters subdivided into the three (four) land use cover types.

#### Reproduction instructions

To reproduce the analysis right from the beginning, you can to define the protected sites (e.g., sites within nature conservation areas) by your own (see [Instructions](#arcgis)).

# References

Associated publication:

[1] Zeiss et al. (2022): Challenges of and opportunities for protecting European soil biodiversity. Conservation Biology. https://doi.org/10.1111/cobi.13930

Related repository:

* JeMaNd_r/ManagementPlans
* FigShare DOI: 10.6084/m9.figshare.16698193
