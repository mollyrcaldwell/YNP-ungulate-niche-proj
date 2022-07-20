#---------------------------------#
#- State Space and SSFs -------####
#--- -----------------------------# 
#---------- 4/6/2022 ---------#

library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(stringr)
library(parallel)
library(circular)
library(MASS)
# library(devtools)
# install_github("jmerkle1/MerkleLab-Code-Repository", subdir="MoveTools", force=TRUE)
library(MoveTools)
library(moveHMM)


#set your working drive
setwd("~/UWyo/PhD project/YNP-ungulate-niche-proj/")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./data/GPS data/Cleaned/allspp_cleanedGPSdata_greenup(mar_to_jun).rds")
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


# Remove burstS for which there was only 1 observation
# this helps because there can't really be a sequence of steps when there is only 1 point in a burst
burts_1pt <- data %>%
  st_drop_geometry() %>%
  count(burst) %>%
  filter(n == 1) %>%
  pull(burst)

data <- data %>%
  filter(!burst %in% burts_1pt)
rm(burts_1pt)

#list data by species
data_spp <- split(data, f = data$species)

#----------------------#
# State space model ####
#----------------------#
# from Morales et al. 2004 and Beyer et al. 2013

# Organize your data ...

## first pull out only the ID, date, and coordinates. you MUST name them exactly these things
## to use the functions in moveHMM
df_spp <- list()

for(i in 1:length(data_spp)){
df_spp[[i]] <- data.frame(ID=data_spp[[i]]$id_yr_seas,
                 Time=data_spp[[i]]$date,
                 X=st_coordinates(data_spp[[i]])[,1],
                 Y=st_coordinates(data_spp[[i]])[,2])
}

names(df_spp) <- names(data_spp)

## this function preps the data to run a state space model
## put in the dataframe you just made, the coordinate type and the name of the coordinate columns
data.hmm_spp <- list()

for(i in 1:length(df_spp)){
data.hmm_spp[[i]] <- prepData(df_spp[[i]], type="UTM", coordNames = c("X", "Y"))
}

names(data.hmm_spp) <- names(data_spp)

## this next step will first plot the movement of each of your individuals (colored by ID)
## and will then plot various movement parameters for each individual animal
# plot(data.hmm, compact = TRUE)  

summary(data.hmm_spp[[1]]) ##gives you summary information

## if you have steps with step length zero, you have to account for that in the models
## by including a zero-mass starting parameter because the strictly positive distributions
## (i.e., gamma) are inadequate to use. You can either account for this by adding in a zero-mass
## parameter, or you can add a tiny amount of positive variation to the steps with length 0
## we are going to do that for similicity sake. For more information about how to deal with zero-
## inflation, you can read about it in the moveHMM vignette.

for(i in 1:length(data.hmm_spp)){
## first check to see if you have any steps with length 0
# table(data.hmm$step==0) ## for my data, there are 22 steps with step length 0

## add in some very small variation for those 22 steps
data.hmm_spp[[i]]$step <- ifelse(data.hmm_spp[[i]]$step==0, 0+runif(1, 0, .001), 
                                 data.hmm_spp[[i]]$step)

## after you add in the variation, these should all be FALSE, if they are not, something is wrong
# table(data.hmm$step==0)
}


## before we run models with different states, we need to set some starting parameters for the models to start from
## we can use the mean and distribution of our actual data to decide on these
## looking at the distribution of steps and turning angles will also let us decide what
## distribution to use
step_angle_summ <- data.frame(step_mean = "NA", step_sd = "NA", angle_mean = "NA",
                              angle_sd = "NA", species = "NA")

for(i in 1:length(data.hmm_spp)){
print(hist(data.hmm_spp[[i]]$step, main = names(data.hmm_spp)[i]))
print(hist(data.hmm_spp[[i]]$angle, main = names(data.hmm_spp)[i]))

step_mean <- mean(data.hmm_spp[[i]]$step, na.rm=T)
step_sd <- sd(data.hmm_spp[[i]]$step, na.rm=T)

angle_mean <- mean(data.hmm_spp[[i]]$angle, na.rm=T)
angle_sd <- sd(data.hmm_spp[[i]]$angle, na.rm=T)

species <- names(data.hmm_spp)[i]

df_temp <- cbind(step_mean, step_sd, angle_mean, angle_sd, species)

step_angle_summ <- rbind(step_angle_summ, df_temp)
}

step_angle_summ <- step_angle_summ %>% filter(species != "NA") %>%
  mutate(step_mean = as.numeric(step_mean)) %>%
  mutate(step_sd = as.numeric(step_sd)) %>%
  mutate(angle_mean = as.numeric(angle_mean)) %>%
  mutate(angle_sd = as.numeric(angle_sd))

#-----------------------#
# DOUBLE STATE MODEL ####
#-----------------------#

## give some starting parameter values, because we are trying to identify 2 states now, you need to put in
## the starting parameters that you think are close for each state.

# identify cores (use 1 less than you have)
no_cores <- detectCores()-1
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("data.hmm_spp", "step_angle_summ"))

doubleStateModel_list <- list()
for(i in 1:length(data.hmm_spp)){
  print(names(data.hmm_spp)[i])
mu0 <- c((step_angle_summ$step_mean[i] - step_angle_summ$step_mean[i]/2),
         (step_angle_summ$step_mean[i] + step_angle_summ$step_mean[i]/2)) # step mean (need one for each state, so should have 2)
sigma0 <- c((step_angle_summ$step_sd[i] - step_angle_summ$step_sd[i]/2),
            (step_angle_summ$step_sd[i] + step_angle_summ$step_sd[i]/2)) # step SD 
stepPar0 <- c(mu0,sigma0) ## combine the mean and SD into one object for the model

angleMean0 <- c((step_angle_summ$angle_mean[i] - step_angle_summ$angle_mean[i]/2),
                (step_angle_summ$angle_mean[i] + step_angle_summ$angle_mean[i]/2)) # angle mean
kappa0 <- c((step_angle_summ$angle_sd[i] - step_angle_summ$angle_sd[i]/2),
            (step_angle_summ$angle_sd[i] + step_angle_summ$angle_sd[i]/2)) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

## run the double state model
doubleStateModel_list[[i]] <- fitHMM(data.hmm_spp[[i]], nbStates = 2, 
                           stepPar0 = stepPar0,
                           anglePar0 = anglePar0, verbose=2,
                           stepDist = "gamma", angleDist = "vm")

}

stopCluster(clust)   # you must stop the parallelization framework

names(doubleStateModel_list) <- names(data_spp)

##can look at model output
doubleStateModel_list[[1]]

## can look at the CI for the model
CI(doubleStateModel_list[[1]])

## then plot the model - this returns a plot of the distribution of step length and and a plot of the turning angle for
## the two states determined by the model; it then will give you an output for each indivdiual animal, with the 
## states colored in orange and blue
plot(doubleStateModel_list[[1]])



## this gives you a series of plots for each individual animal. The top plot shows you if they are in state 1 or 2,
## and the bottom two plots gives you the probability that the animal is in state 1 (middle plot) or state 2 (bottom plot)
#plotStates(doubleStateModel_list[[1]])
rm(angleMean0, anglePar0, kappa0, mu0, sigma0, stepPar0, df, data.hmm)

#--------------------------#
# add most likely state ####
#--------------------------#

# State 1 is the encamped state (short step, large angles)
# State 2 is traveling (long steps, small turning angles)

for(i in 1:length(doubleStateModel_list)){
  print(i)
# Get the state
state <- viterbi(doubleStateModel_list[[i]])

# Get the probability of being in each state
sp <- stateProbs(doubleStateModel_list[[i]])

# Add states to original dataframe
data_spp[[i]]$state <- state
data_spp[[i]]$prob_state1 <- sp[, 1]
data_spp[[i]]$prob_state2 <- sp[, 2]

# Add columns with step length and angle to data
data_spp[[i]]$hmm.dist <- doubleStateModel_list[[i]]$data$step
data_spp[[i]]$hmm.angle <- doubleStateModel_list[[i]]$data$angle

}

rm(sp, state)

table(data_spp[[5]]$state)


#-----------------------------------------#
# organize to identify movement bursts ####
#-----------------------------------------#

for(i in 1:length(data_spp)){
  
print(i)
# Add a new column with Timestamp_EndOfStep for each consecutive steps of each burst
data_spp[[i]] <- data_spp[[i]] %>%
  group_by(burst) %>%
  mutate(date_EndOfStep=c(date[-1], NA)) %>%
  ungroup() %>%
  as.data.frame() %>%   # take out of tibble
  st_set_geometry("geometry")  # remake into sf object

# Loop across all bursts to identify first locations when animals are changing states
data_spp[[i]] <- data_spp[[i]] %>%
  group_by(burst) %>%
  # First location of state 2 (animal is starting to travel) will lead to "1" (from 2-1)
  # First location of state 1 (animal is switching to encamped state) will lead to "-1" (from 1-2)
  mutate(transition=c(NA,   # Add NA for first location, because we don't know if it's transitioning or not!
                      na.omit(lead(state)-state))) %>% # Subtract each line by the line before (which will give us a vector of length(state)-1)
  ungroup() %>%
  as.data.frame() %>%   # take out of tibble
  st_set_geometry("geometry")  # remake into sf object

# table(data_spp[[i]]$transition)
# 
# head(data_spp[[i]])

# Add a column to ID locations that are transitions, and transition into which state
data_spp[[i]] <- data_spp[[i]] %>%
  mutate(first_pt_encamped=ifelse(data_spp[[i]]$transition == -1, 1, 0),
         first_pt_traveling=ifelse(data_spp[[i]]$transition == 1, 1, 0))

# Add a column to identify last location of a burst, and if it's the last of a sequence of the encamped or traveling state
## Last location from a burst is the end of an encamped step
data_spp[[i]] <- data_spp[[i]] %>%
  group_by(burst) %>%
  mutate(last_pt_encamped=lag(ifelse(is.na(lead(dist)) == TRUE & state == 1, 1, 0))) %>% 
  ungroup() %>%
  as.data.frame() %>%   # take out of tibble
  st_set_geometry("geometry")  # remake into sf object

## Last location from a burst is the end of a traveling step
data_spp[[i]] <- data_spp[[i]] %>%
  group_by(burst) %>%
  mutate(last_pt_traveling=lag(ifelse(is.na(lead(dist)) == TRUE & state == 2, 1, 0))) %>% 
  ungroup() %>%
  as.data.frame() %>%   # take out of tibble
  st_set_geometry("geometry")  # remake into sf object

## Add column with the actual behavior, with NAs for the last step of bursts, 
# because we don't have step length nor turning angle of following step.
data_spp[[i]] <- data_spp[[i]] %>%
  mutate(behavior=ifelse(is.na(dist) == TRUE, NA,
                         ifelse(state == 1, "encamped", "traveling")))
}

# number of points actually attributed to a state and also part of a burst
table(data_spp[[3]]$behavior)
table(is.na(data_spp[[3]]$behavior))   # 100 NAs
length(unique(data_spp[[3]]$burst)) # the above should be the same number as this


#----------------------------------#
# keep only traveling movements ####
#----------------------------------#

travelState_list <- list()

for(i in 1:length(data_spp)){
# Keep only traveling steps (these are basically the origin locations of traveling steps):
travelSteps <- data_spp[[i]] %>%
  filter(behavior == "traveling" & is.na(behavior) == FALSE)

# End locations that are the last one of a burst, so for which behavior = NA, but the step before was in the travel state:
lastStepOfBurstTravel <- data_spp[[i]] %>%
  filter(last_pt_traveling == 1 & is.na(last_pt_traveling) == FALSE)

# Keep also locations that are the beginning of a first encamped step, right after a transition from traveling to encamped
# (end locations of traveling steps that are the last of a travel burst):
firstEncampedLocation <- data_spp[[i]] %>%
  filter(first_pt_encamped == 1 & !is.na(first_pt_encamped))

# Bind them together and order data_spp[[i]]
travelState <- rbind(travelSteps, lastStepOfBurstTravel, firstEncampedLocation) 
travelState_list[[i]] <- travelState %>%   #order data_spp[[i]] by id and by date
  dplyr::arrange(id_yr_seas, date)
}

rm(travelSteps, lastStepOfBurstTravel, firstEncampedLocation)
head(travelState_list[[1]])


#########################################################
# 5) IDENTIFY TRAVEL BURSTS #############################
#########################################################

for(i in 1:length(travelState_list)){
  print(i)
# Remove old burst column
travelState_list[[i]] <- travelState_list[[i]] %>%
  dplyr::select(!burst)

# Order travelState_list[[i]]
travelState_list[[i]] <- travelState_list[[i]] %>%   #order data by id and by date
  arrange(id_yr_seas, date)

# Change travelState_list[[i]] out of sf
proj <- st_crs(travelState_list[[i]])   # grab the projection for later
travelState_list[[i]] <- travelState_list[[i]] %>%
  mutate(x=st_coordinates(.)[,1],
         y=st_coordinates(.)[,2]) %>%
  st_drop_geometry()

# head(travelState_list[[i]])
# colnames(travelState_list[[i]])

# NAs in first_pt_traveling column mean that it's the first step of a burst - so we can replace them by 1s
# Doesn't really matter if it's first travel step in reality or not, but we just need to count it as start of a new burst
# So we'll create a new column (first_travel_state4TravelBurst)

# table(is.na(travelState_list[[i]]$first_pt_traveling))

travelState_list[[i]] <- travelState_list[[i]] %>%
  mutate(first_travel_state4TravelBurst=ifelse(is.na(first_pt_traveling) == TRUE, 1, first_pt_traveling))
# table(is.na(travelState_list[[i]]$first_travel_state4TravelBurst)) # There should be no NAs now!

# Add new column with bursts, to identify every bout of travel
# We'll use the transition from traveling to encamped to separate bursts
# i.e., when cumsum() of first_pt_traveling = 1

travelState_list[[i]] <- travelState_list[[i]] %>%
  mutate(SplitID=c(0, diff(as.numeric(as.factor(id_yr_seas)))))
table(travelState_list[[i]]$SplitID) # this should be the number of your ids or id_yr_seas'

# create TravelBurst, which denotes separate sequences of points that are in travel state
travelState_list[[i]] <- travelState_list[[i]] %>%
  mutate(TravelBurst=cumsum(ifelse(SplitID != 0, 1,
                                   ifelse(first_travel_state4TravelBurst >= 1, 1, 0))))

# Make sure there are no NAs
# table(is.na(travelState_list[[i]]$TravelBurst))
# 
# head(travelState_list[[i]], 15)

# How many unique TravelBurst?
# length(unique(travelState_list[[i]]$TravelBurst)) # 1614 unique TravelBursts
# 
# # How many TravelBursts are at the end of a burst?
# sum(is.na(travelState_list[[i]]$Timestamp_EndOfStep == TRUE)) # I have 0, but might need the code below to deal with some

Bursts2Remove <- unique(travelState_list[[i]][is.na(travelState_list[[i]]$Timestamp_EndOfStep) == TRUE, c("TravelBurst")])
# length(Bursts2Remove) 

# Remove them
travelState_list[[i]] <- travelState_list[[i]] %>%
  filter(!TravelBurst %in% Bursts2Remove)
length(unique(travelState_list[[i]]$TravelBurst)) # 3738-114 = 3624, perfect!
rm(Bursts2Remove)

# remove unused columns
travelState_list[[i]]$dist <- NULL
travelState_list[[i]]$dt <- NULL
travelState_list[[i]]$speed <- NULL
travelState_list[[i]]$abs.angle <- NULL
travelState_list[[i]]$rel.angle <- NULL
travelState_list[[i]]$StepFlag <- NULL
travelState_list[[i]]$splitID <- NULL

}

names(travelState_list) <- names(data_spp)
#----------------------------------------------------#
# summary information about each traveling burst #####
#----------------------------------------------------#
# Create new dataframe where each step represents the beginning and end 
# of a traveling sequences (ie.g,. burst) of poitns 
travel_summ_df <- data.frame(mean_steps = "NA", sd_steps = "NA", med_steps = "NA",
                             species = "NA")
for(i in 1:length(travelState_list)){
# First, order dataframe
travelState_list[[i]] <- travelState_list[[i]] %>%   #order data by id and by date
  arrange(id_yr_seas, date)

UniqueBursts <- travelState_list[[i]] %>%
  group_by(TravelBurst) %>%
  summarise(TravelBurst = first(TravelBurst), # Unique burst ID
            id_yr_seas = first(id_yr_seas),
            date_orig = first(date), # Timestamp of the first step of the burst (first row, 6th column)
            date_end = last(date), # Timestamp at the end of the last step of the burst
            x_orig = first(x), # x_orig from that first step of the burst
            y_orig = first(y), # y_orig from that first step of the burst
            x_end = last(x), # x_end from that last step of the burst
            y_end = last(y), # y_end from that last step of the burst
            NbrTravelSteps = n()-1) %>% # The number of consecutive traveling steps in that burst (-1 because that last location is the start of an encamped step))
  as.data.frame()   # take out of tibble


# head(UniqueBursts)
# str(UniqueBursts)

# length(unique(travelState_list[[i]]$TravelBurst)) # 1614 unique TravelBursts
# nrow(UniqueBursts)  # same number as above

#table(is.na(UniqueBursts)) # No NAs at all!


# How long do traveling bursts last?
# table(UniqueBursts$NbrTravelSteps)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  17  20  23  28 
# 124 238 308 351 320 185  34   7   4   3   2   8  12   9   5   1   1   1   1
mean_steps <- mean(UniqueBursts$NbrTravelSteps)
sd_steps <- sd(UniqueBursts$NbrTravelSteps)
# mean = 4.065056
med_steps <- median(UniqueBursts$NbrTravelSteps)
# median = 4
rm(UniqueBursts)
species <- names(travelState_list)[i]

df <- cbind(mean_steps, sd_steps, med_steps, species)

travel_summ_df <- rbind(travel_summ_df,df)
}

travel_summ_df <- travel_summ_df %>% filter(species != "NA")

write.table(travel_summ_df, file ="./Code output/travel_steps_summ_state space ssf.txt",
            sep = ",", quote = FALSE, row.names = FALSE)

# -------------------------------------
# Create a similar dataframe, but more prepared for the moveparams function 
# -------------------------------------
data.orig_list <- data_spp

for(i in 1:length(data_spp)){

bursts.df.head <- travelState_list[[i]] %>%
  group_by(TravelBurst) %>%
  slice_head(n=1) %>%
  ungroup()

bursts.df.tail <- travelState_list[[i]] %>%
  group_by(TravelBurst) %>%
  slice_tail(n=1) %>%
  ungroup()

data_spp[[i]] <- rbind(bursts.df.head, bursts.df.tail)
data_spp[[i]] <- data_spp[[i]] %>%   #order data by id and by date
  arrange(id_yr_seas, date) %>% 
  as.data.frame()
rm(bursts.df.head, bursts.df.tail)
}

# Look at final dataframe
head(data_spp[[1]])
colnames(data_spp[[1]])

#-------------#
# plotting ####
#-------------#


# some plotting to see what you did!

# Make sf dataframe
data <- st_as_sf(data, coords = c("x","y"), dim = "XY", crs = proj)


burst2plot <- sample(unique(data$TravelBurst),1)
hrs_around_burst <- 0   # how many hours before and after the burst2plot to plot

dat.tmp <- data[data$TravelBurst == burst2plot,]
dat.orig.tmp <- data.orig %>% 
  filter(id_yr_seas == dat.tmp$id_yr_seas[1] & 
           date >= dat.tmp$date[1]-hrs_around_burst*3600 &
           date <= dat.tmp$date[2]+hrs_around_burst*3600)

mapview(dat.orig.tmp, zcol="behavior", layer.name="origional points") +
  mapview(st_cast(st_combine(dat.orig.tmp), "LINESTRING"), color="grey", layer.name="origional line") +
  mapview(dat.tmp, layer.name="new burst points") +
  mapview(st_cast(st_combine(dat.tmp), "LINESTRING"), layer.name="new burst line")
  

rm(burst2plot, hrs_around_burst, dat.tmp, dat.orig.tmp)


#--------------------------#
# sampling random steps ####
#--------------------------#

# No need to calculate bursts again! We already have unique TravelBurst!
length(unique(data$TravelBurst)) # 1615 unique TravelBurst

table(table(data$TravelBurst))  # These should all be 2



for(i in 1:length(data_spp)){
   data_spp[[i]] <- st_as_sf(data_spp[[i]], coords = c("x", "y"))
# Calculate movement parameters
data_spp[[i]] <- CalcMovParams(data = data_spp[[i]], id_name = "id_yr_seas", 
                      date_name = "date", burst = data_spp[[i]]$TravelBurst)


# This is where you say 'I don't want these steps connected for further development of the SSF'
# Note from JAM: Anytime there is a FALSE in the StepFlag column, my 'step sampling' code knows NOT to include the step
#sum(data_spp[[i]]$StepFlag) # How many 'actual/usable' steps you currently have - there shouldn't be any! 

# Trick data_spp[[i]]frame so we can use Jerod's function to draw random steps
## Change StepFlag column for TRUE and FALSE
data_spp[[i]] <- data_spp[[i]] %>%
  mutate(StepFlag= ifelse(is.na(dist)==FALSE, TRUE, FALSE))
# head(data_spp[[i]])
# sum(data_spp[[i]]$StepFlag) # Should be back to your number of travelbursts
# length(unique(data_spp[[i]]$TravelBurst))


# hist(data_spp[[i]]$dist[data_spp[[i]]$StepFlag == TRUE]/1000)
# mean(data_spp[[i]]$dist[data_spp[[i]]$StepFlag == TRUE]/1000)  # mean 0.9176365
# median(data_spp[[i]]$dist[data_spp[[i]]$StepFlag == TRUE]/1000) # median 0.580640
# quantile(data_spp[[i]]$dist[data_spp[[i]]$StepFlag == TRUE]/1000, probs=0.95) 

# If you want to use weibull or gamma distribution for step lengths, 
# need to convert StepFlag to FALSE when dist or speed is 0
# table(data_spp[[i]]$speed[data_spp[[i]]$StepFlag == TRUE] <= 0)   # all false
# data_spp[[i]]$StepFlag[data_spp[[i]]$speed <= 0] <- FALSE  # Remove those if needed How can there be steps of length 0 with the HMM process!!!!! Loop and come back, i guess.
# table(data_spp[[i]]$speed[data_spp[[i]]$StepFlag == TRUE] <= 0) # Just to check again
# table(data_spp[[i]]$dist[data_spp[[i]]$StepFlag == TRUE] <= 0) # Just to check again

# Since they weren't calculated because there aren't 'previous steps' to determine it, put some random turning angles
data_spp[[i]]$rel.angle <- data_spp[[i]]$abs.angle

# # take out of sf, and keep simply as a data_spp[[i]]frame
data_spp[[i]] <- data_spp[[i]] %>%
  mutate(x = st_coordinates(.)[,1], # add x and y as columns
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()  #need to remove the geometry column


}
 # ---------------------------------------------#
# sampling based on parametric distribution ####
# ---------------------------------------------#

#create variation for 0 speeds 
for(i in 1:length(data_spp)){
data_spp[[i]]$dist<- ifelse(data_spp[[i]]$dist==0, 0+runif(1, 0, .001), 
                             data_spp[[i]]$dist)

data_spp[[i]]$speed<- ifelse(data_spp[[i]]$speed==0, 0+runif(1, 0, .001), 
                            data_spp[[i]]$speed)
}


data_spp_sampled <- list()
for(i in 1:length(data_spp)){
  print(i)
data_spp_sampled[[i]] <- DrawRandSteps_parametric(data=data_spp[[i]], nr=5,   # nr = number of random steps
                                  step_distr = "weibull", # Can choose weibull, gamma, or exponential (for step length/speed)
                                  ta_distr="vonmeses", # Can choose wrappedcauchy or vonmeses (for turning angles)
                                  speed=T,        # Check doing tru and doing false, could be big differences
                                  withinID=FALSE,        # should sample from same individual or other individuals in dataset?
                                  id_name="id_yr_seas", date_name="date", x_name="x", y_name="y",   #what are the names of your id and date columns?
                                  uniform.angles=TRUE, distMax = Inf)

}

head(data,12)

hist(data$dist[data$case == 1])
hist(data$dist[data$case == 0])
hist(data$rel.angle[data$case==1])
hist(data$rel.angle[data$case==0])



# now that you've sampled your available steps, write out your environment! ####
# save your environment so we can load it back up another time ####
save.image("StateSpaceSSFs_AvailSampled.RData")
# load("StateSpaceSSFs_AvailSampled.RData")

# you will now need to extract some GIS data and fit models using your code from labs 5 and 6. 

