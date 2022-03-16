#--------------------------------------#
#- Step Selection Functions --------####
#---Availability sampling -------------# 
#---------- Jerod Merkle --------------#
#------------- Lab 5a -----------------#

library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(stringr)
library(parallel)
library(circular)
library(MASS)
 #library(devtools)
 #install_github("jmerkle1/MerkleLab-Code-Repository", subdir="MoveTools", force=TRUE)
library(MoveTools)


#set your working drive
setwd("~/UWyo/PhD project/YNP-ungulate-niche-proj/")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./Data/GPS data/Cleaned/allspp_cleanedGPSdata_greenup(mar_to_jun).rds")
head(data)


# do you have duplicates in your data
dup <- data %>% 
  st_drop_geometry() %>% 
  dplyr::select(id_yr_seas, date) %>% 
  duplicated()
table(dup)
data <- data %>% 
  filter(dup == FALSE)   #this is how to remove duplicates. But be careful. Pay attention to what you are doing here!
rm(dup)

#----------------------------------#
# calculate movement parameters ####
#----------------------------------#
data <- data %>% 
  arrange(id_yr_seas, date) #order database first
data$burst <- CalcBurst(data=data, id = TRUE, id_name="id_yr_seas", 
                        date_name="date", Tmax = 3600*3) # now I'm setting mine to just more than my fix rate, as I don't want to connect steps where missed fixes occurred
length(unique(data$burst))
# now for movement params
data <- CalcMovParams(data=data, id_name = "id_yr_seas", 
                      date_name = "date", burst=data$burst)

hist(data$dt/3600)
table(data$dt)   # check to make sure all your steps are the same or very similar length
table(data$dt/3600) # now in hrs   
# Note, because of small variations in when the collar actually takes a location,
# you may see a variety of values here. All your steps should be relatively close in dt.
# For example, if you have 3 hour data, you only want to use steps that 
# are between 2.8 and 3.2 hours apart. 
# This is where you say 'I don't want these steps connected for further development of the SSF'
# Anytime there is a FALSE in the StepFlag column, my 'step sampling' code knows NOT to include the step
#you will need to change the next two lines of code !!!!
sum(data$StepFlag)   # this is how many 'actual/usable' steps you currently have
data$StepFlag[data$dt < 3600*1.8] <- FALSE   # this says 'don't connect steps less than 1.8 hrs apart
data$StepFlag[data$dt > 3600*2.2] <- FALSE   # this says 'don't connect steps greater than 2.2 hrs apart
sum(data$StepFlag)   # this is how many 'actual/usable' steps you now have

table(data$dt[data$StepFlag == TRUE])

# are your turning angles and step lengths correlated?
head(data)
#loop by species
spp <- unique(data$species)

for(i in 1:length(spp)){
  data_spp <- data %>% filter(species == spp[i])
tarad <- rad(data_spp$rel.angle)
hist(tarad, main = spp[i])
# a refresher on what sin and cosine does to radians
#plot(seq(0,2*pi,length.out=100), sin(seq(0,2*pi,length.out=100)))  
#plot(seq(0,2*pi,length.out=100), cos(seq(0,2*pi,length.out=100)))
tasin <- sin(tarad)   # this is basically a metric of left and right turns
tacos <- cos(tarad)   # this is a metric of forward and backward turns (positive values are forward, negative are backward)
print(paste0("correl. right/left and step length ", spp[i], ": ",
       cor(tasin,data_spp$dist, use="pairwise.complete.obs", method="pearson")))  #this is your pearson correlation for step length and going right versus left
print(paste0("correl. back/forward and step length ", spp[i], ": ",
       cor(tacos,data_spp$dist, use="pairwise.complete.obs", method="pearson")))   #this is your pearson correlation for step length and going forward versus backward
rm(tarad,tasin,tacos)
}

#--------------------------#
# sampling random steps ####
#--------------------------#

# take out of sf, and keep simply as a dataframe
proj <- st_crs(data)   # grab the projection for later
data <- data %>% 
  mutate(x = st_coordinates(.)[,1], # add x and y as columns
         y = st_coordinates(.)[,2]) %>% 
  st_drop_geometry()  #need to remove the geometry column
head(data)
class(data)


#------------------------------------------#
# sample based on empirical distribution ####
#------------------------------------------#
#loop by species
data_list <- list()
spp <- unique(data$species)

for(i in 1:length(spp)){
data_spp <- data %>% filter(species == spp[[i]]) %>% 
  filter(dist > 100) #only take steps > 100 m

data_list[[i]] <- DrawRandSteps_emperical(data=data_spp, nr=5,   # nr = number of random steps 
                                simult=FALSE,        # should it take a step length and a turning angle simultaneously?
                                id_name="id_yr_seas", date_name="date", x_name="x", y_name="y",   #what are the names of your id and date columns?
                                withinID=TRUE,        # should sample from same individual or if FALSE all individuals in dataset?
                                distMax = Inf, uniform.angles=FALSE)

}

names(data_list) <- spp

head(data_list[[1]], 20)
str(data_list[[1]])

for(i in 1:length(data_list)){
hist(data_list[[i]]$dist[data_list[[i]]$case == 1], 
     xlab = "distance, case = 1", main = spp[i], breaks = 20)
hist(data_list[[i]]$dist[data_list[[i]]$case == 0], 
     xlab = "distance, case = 0", main = spp[i], breaks = 20)
hist(data_list[[i]]$rel.angle[data_list[[i]]$case==1], 
     xlab = "turn angle, case = 1", main = spp[i])
hist(data_list[[i]]$rel.angle[data_list[[i]]$case==0],
     xlab = "turn angle, case = 0", main = spp[i])
}


#combine data list to one dataframe
data <- rbind(data_list[[1]], data_list[[2]], data_list[[3]], data_list[[4]],
              data_list[[5]])


###JUST used empirical above, did not run below code###
# ---------------------------------------------#
# sampling based on parametric distribution ####
# ---------------------------------------------#
# if you want want... 

# if you want to use weibull or gamma distribution for step lengths, 
# need to convert StepFlag to FALSE when dist or speed is 0
for(i in 1:length(data_list)){
table(data_list[[i]]$dist[data_list[[i]]$StepFlag == TRUE] == 0)   # how many of these do you have?  If more than a few, need to think about this and discuss
data_list[[i]]$StepFlag[data_list[[i]]$dist <= 0] <- FALSE  # remove those or add some small value to the step length
table(data_list[[i]]$speed[data_list[[i]]$StepFlag == TRUE] == 0) # just to check


# figure out which distribution fits your step length data_list[[i]] best
fit_weib <- fitdistr(data_list[[i]]$dist[data_list[[i]]$StepFlag == TRUE],densfun="weibull")
fit_gamma <- fitdistr(data_list[[i]]$dist[data_list[[i]]$StepFlag == TRUE],densfun="gamma")
fit_exp <- fitdistr(data_list[[i]]$dist[data_list[[i]]$StepFlag == TRUE],densfun="exponential")

fit_weib$n; fit_gamma$n; fit_exp$n  # these need to all be the exact same to compare AIC
aic_tbl <- AIC(fit_weib, fit_gamma, fit_exp)  
aic_tbl[order(aic_tbl$AIC),]   # look for lowest AIC
rm(fit_exp,fit_gamma,fit_weib,aic_tbl)
}

data <- DrawRandSteps_parametric(data=data, nr=5,   # nr = number of random steps
                                  step_distr = "weibull", # Can choose weibull, gamma, or exponential (for step length/speed)
                                  ta_distr="vonmeses", # Can choose wrappedcauchy or vonmeses (for turning angles)
                                  speed=FALSE,        # should we use speed because your steps are each not nearly the same dt?
                                  withinID=FALSE,        # should sample from same individual or other individuals in dataset?
                                  id_name="id_yr_seas", date_name="date", x_name="x", y_name="y",   #what are the names of your id and date columns?
                                  uniform.angles=FALSE, distMax = Inf)
head(data,12)

hist(data$dist[data$case == 1])
hist(data$dist[data$case == 0])
hist(data$rel.angle[data$case==1])
hist(data$rel.angle[data$case==0])



# now that you've sampled your available steps, write out your environment! ####
# save your environment so we can load it back up another time ####
save.image("SSFs_AvailSampled.RData")
# load("SSFs_AvailSampled.RData")



