#--------------------------------------#
#- Step Selection Functions --------####
#--- extract GIS data to steps --------# 
#---------- Jerod Merkle --------------#
#------------- Lab 5b ------------------#

library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(lubridate)
library(stringr)
library(parallel)
library(snow)
library(circular)
library(MASS)
library(MoveTools)


#set your working drive
setwd("~/UWyo/PhD project/YNP-ungulate-niche-proj/")

#---------------------------------------------#
# Load SSF data with availability sampled  ####
#---------------------------------------------#
# from lab SSF part a
load("SSFs_AvailSampled.RData")

#reduce data to 10 random individuals per species
spp <- unique(data$species)
spp_sample <- list()
for(i in 1:length(spp)){
  data_spp <- data %>% filter(species == spp[i])
  samp <- sample(unique(data_spp$id_yr_seas), 10)
  
  spp_sample[[i]] <- data %>% filter(id_yr_seas %in% samp)
}

data <- rbind(spp_sample[[1]], spp_sample[[2]], spp_sample[[3]],
                       spp_sample[[4]], spp_sample[[5]])

#--------------------------------------#
# reduce database to 1 step per day ####
#--------------------------------------#
# to reduce pseudoreplication (you don't have to do this, or you can do something different too)
# also doing this increases speed to fit the models later
# get id_yr_seas column fixed, so there are no NAs in it
data <- data %>% 
  arrange(id_yr_seas, date, case)

# grab one step per day, per ID, per year, per season...
# add some new columns
data <- data %>% 
  mutate(jul=yday(date),
         id_yr_seas_jul=paste(id_yr_seas, jul, sep="_"))

# create a reduced database of used points only that is randomized
tokeep <- data %>% 
  filter(case==1) %>% 
  dplyr::select(strata, id_yr_seas_jul) %>% 
  arrange(sample(nrow(.),nrow(.),replace=FALSE))  # randomize rows of tokeep

# grab the first step and its random steps per day
data <- data %>% 
  filter(strata %in% tokeep$strata[duplicated(tokeep$id_yr_seas_jul)==FALSE])

data$strata <- as.numeric(as.factor(data$strata))  # recreate ordered/numeric strata column
data$id_yr_seas_jul <- NULL   # remove the id_yr_seas column
data$jul <- NULL   # remove the id_yr_seas column
rm(tokeep)
head(data, 10)

n_distinct(data$strata)  # how many used steps are you left with? I have 1,819
# this gives you number of unique steps per id_yr_seas. I have about 122 each.
data %>% 
  group_by(id_yr_seas) %>% 
  dplyr::summarise(numb.strata=n_distinct(strata))


# -----------------------------------#
# Make data spatially aware again ####
# -----------------------------------#
head(data)
data <- data %>% 
  dplyr::rename(x_orig=x, #change x/y names, so we have x_orig and x_end of the step
                y_orig=y)

#turn into point dataframe and make spatially aware
data <- data %>% 
  arrange(id_yr_seas, date, case)
table(is.na(data$x_end))  #make sure this is all FALSE

# create a separate sf object for the source point (i.e., orig) and target point (i.e., end)
#for the end of steps
sf_end <- st_as_sf(data, coords=c("x_end","y_end"), dim="XY", crs=proj)
#for the start of steps (we'll need later)
sf_orig <- st_as_sf(data,coords=c("x_orig","y_orig"), dim="XY", crs=proj)

#---------------------------#
# bring in your GIS data ####
#---------------------------#
#bring in all your GIS data
elev <- raster("~/UWyo/PhD project/HOTR-proj/Data/GIS_data_clean/Elevation_meters_30m_gye.tif")
trasp <- raster("~/UWyo/PhD project/HOTR-proj/Data/GIS_data_clean/Aspect_TRASP_30m_gye.tif")
slope <- raster("~/UWyo/PhD project/HOTR-proj/Data/GIS_data_clean/Slope_degrees_30m_gye.tif")
tpi <- raster("~/UWyo/PhD project/HOTR-proj/Data/GIS_data_clean/TPI_unitless_30m_gye.tif")

#RAP tree cover and forb/grass biomass by year
#RAP data (forbs grass biomass and %cover by land cover type per year)
RAP_ann_biomass_file_list <- list.files(path = "~/UWyo/PhD project/HOTR-proj/Data/GIS_data_clean/", 
                                        pattern = glob2rx("*Biomass*AnnualForbsGrasses*"),
                                        full.names = TRUE)

RAP_per_biomass_file_list <- list.files(path = "~/UWyo/PhD project/HOTR-proj/Data/GIS_data_clean/", 
                                        pattern = glob2rx("*Biomass*PerennialForbsGrasses*"),
                                        full.names = TRUE)

RAP_treecov_file_list <- list.files(path = "~/UWyo/PhD project/HOTR-proj/Data/GIS_data_clean/", 
                                    pattern = glob2rx("*Cover_Trees*"),
                                    full.names = TRUE)

RAP_ann_biomass_list <- list()
for(i in 1:length(RAP_ann_biomass_file_list)){
  print(i)
  RAP_ann_biomass_list[[i]] <- raster(RAP_ann_biomass_file_list[[i]])
}

RAP_ann_biomass_names <- sub(".*RAP_gye_30m_reproj_", "", RAP_ann_biomass_file_list)
RAP_ann_biomass_names <- sub("*.tif", "", RAP_ann_biomass_names)
names(RAP_ann_biomass_list) <- RAP_ann_biomass_names

RAP_per_biomass_list <- list()
for(i in 1:length(RAP_per_biomass_file_list)){
  print(i)
  RAP_per_biomass_list[[i]] <- raster(RAP_per_biomass_file_list[[i]])
}

RAP_per_biomass_names <- sub(".*RAP_gye_30m_reproj_", "", RAP_per_biomass_file_list)
RAP_per_biomass_names <- sub("*.tif", "", RAP_per_biomass_names)
names(RAP_per_biomass_list) <- RAP_per_biomass_names

RAP_treecov_list <- list()
for(i in 1:length(RAP_treecov_file_list)){
  print(i)
  RAP_treecov_list[[i]] <- raster(RAP_treecov_file_list[[i]])
}

RAP_treecov_names <- sub(".*RAP_gye_30m_reproj_", "", RAP_treecov_file_list)
RAP_treecov_names <- sub("*.tif", "", RAP_treecov_names)
names(RAP_treecov_list) <- RAP_treecov_names

#add perennial and annual grass and forb biomass rasters by year
RAP_forbgrass_biomass <- list()

for(i in 1:length(RAP_per_biomass_list)){
  print(paste0("i: ", i))
  
  biomass_stack <- stack(RAP_ann_biomass_list[[i]], RAP_per_biomass_list[[i]])
  ann_per_biomass <-  calc(biomass_stack, sum)
  
  RAP_forbgrass_biomass[[i]] <- ann_per_biomass
}

names(RAP_forbgrass_biomass) <- c("forb_grass_biomass_2016", 
                                  "forb_grass_biomass_2017",
                                  "forb_grass_biomass_2018",
                                  "forb_grass_biomass_2019",
                                  "forb_grass_biomass_2020")

#duplicate 2020 rap biomass and tree cov as 6th list object in order to run extract raster on 2021 data
RAP_forbgrass_biomass[[6]] <- RAP_forbgrass_biomass[[5]]
RAP_treecov_list[[6]] <- RAP_treecov_list[[5]]


# ------------------------------------#
# relate target points to GIS data ####
# ------------------------------------#

# sf_end or target points first (all variables)
data$elev_target <- raster::extract(elev, st_transform(sf_end, crs=projection(elev)))
data$trasp_target <- raster::extract(trasp, st_transform(sf_end, crs=projection(trasp)))
data$slope_target <- raster::extract(slope, st_transform(sf_end, crs=projection(slope)))
data$tpi_target <- raster::extract(tpi, st_transform(sf_end, crs=projection(tpi)))

#extract rap data by year of gps dates
#split data into list based on years, then rbind data back together at end
yrs <-c(2017, 2018, 2019, 2020, 2021) #sampled data doesn't have 2016
data_yr <- list()

for(i in 1:length(yrs)){
  data_yr[[i]] <- data %>% filter(year(date) == yrs[[i]])
}


for(i in 1:length(data_yr)){
  data_yr[[i]]$forbgrass_biomass_target <- raster::extract(RAP_forbgrass_biomass[[i+1]], #account for missing 2016
                                           st_transform(sf_end %>% filter(year(date) == yrs[[i]]), 
                                                        crs=projection(RAP_forbgrass_biomass[[i+1]])))
  
  data_yr[[i]]$perc_treecov_target <- raster::extract(RAP_treecov_list[[i+1]], 
                                        st_transform(sf_end %>% filter(year(date) == yrs[[i]]),
                                                    crs=projection(RAP_treecov_list[[i+1]])))
}

#bind together each listed data per year as 1 data set
data <- rbind(data_yr[[1]], data_yr[[2]], data_yr[[3]],
                   data_yr[[4]], data_yr[[5]])



head(data)

#use Jerod's ndviExtract code to extract max IRG date, IRG values, annual NDVI to points
source("Z:/MODIS_NDVI/info/ndviEXTRACT.R")

data$maxIRGdate_target <- ndviEXTRACT(data, NDVImetric = "maxIRGdate", xname = "x_end",
                                   yname = "y_end", 
                                   NDVIfolder = "Z:/MODIS_NDVI/Bischof_calculations/",
                                   maxcpus = 11, datesname = "date", scaleIRG = TRUE,
                                   xyCRS = crs(sf_end))

data$IRGVals_target <- ndviEXTRACT(data, NDVImetric = "IRGVals", xname = "x_end",
                                      yname = "y_end", 
                                      NDVIfolder = "Z:/MODIS_NDVI/Bischof_calculations/",
                                      maxcpus = 11, datesname = "date", scaleIRG = TRUE,
                                      xyCRS = crs(sf_end))

data$csumNDVImax_target <- ndviEXTRACT(data, NDVImetric = "csumNDVImax", xname = "x_end",
                                yname = "y_end", 
                                NDVIfolder = "Z:/MODIS_NDVI/Bischof_calculations/",
                                maxcpus = 11, datesname = "date", scaleIRG = TRUE,
                                xyCRS = crs(sf_end))


data$try <- ndviEXTRACT(data, NDVImetric = "csumNDVImax", xname = "x_end",
                                       yname = "y_end", 
                                       NDVIfolder = "Z:/MODIS_NDVI/Bischof_calculations/",
                                       maxcpus = 11, datesname = "date", scaleIRG = TRUE,
                                       xyCRS = crs(sf_end))

#calculate absolute value of days from peak IRG 
data <- data %>%
  mutate(daysfromPeakIRG = abs(yday(date) - maxIRGdate_target))


#-----------------------------#
# relate lines to GIS data ####
#-----------------------------#
# first thing is to create sf lines dataframe with each row representing a step (used or available)

dataL <- lapply(1:nrow(data), function(i){
  return(st_linestring(rbind(as.matrix(data[i,c("x_orig","y_orig")]),
               as.matrix(data[i,c("x_end","y_end")])),dim="XY"))
})
dataL <- st_as_sfc(dataL, crs=proj)
dataL <- data.frame(strata=data$strata, case=data$case, geometry=dataL)
dataL <- st_as_sf(dataL, sf_column_name = "geometry")
dataL$date <- data$date
head(dataL)
nrow(data); nrow(dataL)  # these of course should be exactly the same!

# check to make sure each dataset has the same order of the rows (use strata and case to determine this)
head(data)
head(sf_end)
head(sf_orig)
head(dataL)


#take a look at the steps you have generated ####
smpl <- sample(unique(data$strata),1)   #grab a random strata to plot
plot(dataL$geometry[data$strata==smpl], main=paste("Strata:", smpl))
points(data$x_orig[data$strata==smpl], data$y_orig[data$strata==smpl], col="blue", pch=16)
points(data$x_end[data$strata==smpl], data$y_end[data$strata==smpl], col="orange", pch=16)
points(data$x_end[data$strata==smpl&data$case==1], data$y_end[data$strata==smpl&data$case==1], col="darkgreen", pch=16)
legend("topright", c("Source","Avail","Used"), pch=16, col=c("blue","orange","darkgreen"), bty="n")

# now do the same thing with mapview
smpl <- sample(unique(data$strata),1)   #grab a random strata to plot
mapview(dataL$geometry[data$strata==smpl], layer.name= smpl,
        map.types="Esri.WorldImagery")+
  mapview(sf_end[data$strata==smpl,], zcol="case", layer.name="Case")

rm(smpl)
  
#----------------------------------------#
# Now to relate the lines to GIS data ####
#----------------------------------------#

# you might consider trying your code on a small subset of lines first (make sure its works)
raster::extract(elev,st_transform(dataL[1:25,],crs=projection(elev)), fun=mean) #extracts the cell values that the line touches, then takes the mean of those values

# now to do all of your continuous data!
head(data)
beginCluster(type="SOCK")   #use multi core processing because it takes a while to assess what cells each line are touching

# Note: These take about 20 minutes per line / extraction (I have 11 cores, 30m res rasters, and 11k lines)
Sys.time()
data$elev_step <- raster::extract(elev, st_transform(dataL, crs=projection(elev)), fun=mean)
Sys.time()
data$trasp_step <- raster::extract(trasp, st_transform(dataL, crs=projection(trasp)), fun=mean)
Sys.time()
data$slope_step <- raster::extract(slope, st_transform(dataL, crs=projection(slope)), fun=mean)
data$tpi_step <- raster::extract(tpi, st_transform(dataL, crs=projection(tpi)), fun=mean)

#extract rap data by year of gps dates
#split data into list based on years, then rbind data back together at end
yrs <-c(2017, 2018, 2019, 2020, 2021) #sampled data doesn't have 2016
data_yr <- list()

for(i in 1:length(yrs)){
  data_yr[[i]] <- data %>% filter(year(date) == yrs[[i]])
}


for(i in 1:length(data_yr)){
  data_yr[[i]]$forbgrass_biomass_step <- raster::extract(RAP_forbgrass_biomass[[i+1]], #account for missing 2016
                                                           st_transform(dataL %>% filter(year(date) == yrs[[i]]), 
                                                                        crs=projection(RAP_forbgrass_biomass[[i+1]])), fun = mean)
  
  data_yr[[i]]$perc_treecov_step <- raster::extract(RAP_treecov_list[[i+1]], 
                                                      st_transform(dataL %>% filter(year(date) == yrs[[i]]),
                                                                   crs=projection(RAP_treecov_list[[i+1]])), fun = mean)
 
 
}

Sys.time()

#bind together each listed data per year as 1 data set
data <- rbind(data_yr[[1]], data_yr[[2]], data_yr[[3]],
              data_yr[[4]], data_yr[[5]])



endCluster()
head(data, 3)




# have to do some different things for descrete variables like landcover...
# basically, for each landcover variable, you need to create a seaprate raster of 
# whether or not it is that landcover (i.e., 1s and 0s). Then, when you extract
# to lines, you get the proportion of the line (i.e., step) that crosses a given landcover type
# my suggestion is to pick most important landcover variables here
# unique(landcov)   # this is nlcd landcov
# table(data$lc_target)   # you may need to bring in the legend to interpret
# table(data$landcov_target)   # you may need to bring in the legend to interpret
# barren <- landcov  #make empty landcov database called barron
# barren[] <- ifelse(values(landcov)==31, 1, 0)   #reclassify values so barren landcovers = 1 and all else = 0
# shrub <- landcov
# shrub[] <- ifelse(values(landcov)==52, 1, 0)
# herb <- landcov
# herb[] <- ifelse(values(landcov)%in% 71:74, 1, 0)
# plot(barren)
# plot(shrub)
# plot(herb)
# 
# beginCluster(type="SOCK")  
# data$barren_step <- raster::extract(barren, st_transform(dataL, crs=projection(barren)), fun=mean)
# data$shrub_step <- raster::extract(shrub, st_transform(dataL, crs=projection(shrub)), fun=mean)
# data$herb_step <- raster::extract(herb, st_transform(dataL, crs=projection(herb)), fun=mean)
# endCluster()
# head(data)
# 
# rm(elev, trasp, slope, tpi, Dist2Escape, landcov,
#    Dist2roads,treecov,iNDVI19,Cover_AnnualForbGrass,
#    Cover_BareGround,Cover_PerennialForbGrass,herb,shrub,barren)


# --------------------------------------------------------------#
# Calc whether a step crossed a linear feature, and how many ####
# --------------------------------------------------------------#
#load up your road layer
# roads <- st_read("./data/GIS", "roads_tiger_2016")
# crosses <- st_crosses(dataL, roads, sparse=FALSE)  # matrix of steps and the different roads and whether or not a cross occured
# dim(crosses)
# data$numb_roads_crossed <- apply(crosses, 1, sum)  # sum up across the different roads
# hist(data$numb_roads_crossed)
# table(data$numb_roads_crossed)
# table(data$numb_roads_crossed[data$case == 1])
# table(data$numb_roads_crossed[data$case == 0])
# rm(roads, crosses)

# now that you've extracted your GIS data, write out your environment! ####
# save your environment so we can load it back up another time
save.image("SSFs_VariablesExtracted.RData")
# load("SSFs_VariablesExtracted.RData")


