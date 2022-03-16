#--------------------------------------#
#- Step Selection Functions --------####
#--- Initial analyses -----------------# 
#---------- Jerod Merkle --------------#
#------------- Lab 5c ------------------#

library(dummies)
library(survival)
library(car)
library(raster)
library(sf)
library(tidyverse)
library(ggplot2)
library(ggeffects)
library(lme4)
library(itsadug)
library(visreg)
library(MoveTools)

#set your working drive
setwd("~/UWyo/PhD project/YNP-ungulate-niche-proj/")

#--------------------------------------------#
# Load SSF data with variables extracted  ####
#--------------------------------------------#
# from lab RSF part b
load("SSFs_VariablesExtracted.RData")
head(data)

#---------------------------------------#
#--- Get ready for Analysis ---------####
#---------------------------------------#

# some simple scaling so your variables are on similar scales
hist(data$dist)
data$dist <- data$dist/1000   #get your distance between source and target point variable into km
hist(data$dist)

data$elev_target <- data$elev_target/1000   # get into KM

data$elev_step <- data$elev_step/1000   # get into KM

data$rank_forbgrassbiomass_target <- (data$forbgrass_biomass_target/max(data$forbgrass_biomass_target))
data$rank_perctreecov_target <- data$perc_treecov_target/100
data$rank_forbgrassbiomass_step <- lapply(data$forbgrass_biomass_step, function(value) value/max(unlist(data$forbgrass_biomass_step)))
data$rank_perctreecov_step <- lapply(data$perc_treecov_step, function(value) value/100)
data$rank_csumNDVImax_target <- data$csumNDVImax_target/max(data$csumNDVImax_target, na.rm = T)
data$rank_IRGVals_target <- data$IRGVals_target/max(data$IRGVals_target, na.rm = T)

# # ---------------------------------------#
# # some work on the landcover variable ####
# # ---------------------------------------#
# #fix lc variable, so we have actual variable names for landcover
# # landcover variables are often numeric, so want to merge with their meaning
# legend <- read.csv("./data/GIS/nlcd_legend.csv") #bring in landcov legend
# 
# legend <- legend %>% 
#   dplyr::rename(lc_target=value,  # rename value column to match your landcover column name in data
#          landcov_target=class_general)  # rename the class column to the name you want to see in data
# 
# # note that this is ONLY for the target point
# data <- data %>%  # add a new column with actual landcover values 
#   left_join(legend[,c("landcov_target","lc_target")])
# 
# rm(legend)
# head(data)

# some basic selection ratios of your landcover (categorical variable)
# this will help you 1) get a feel for what is selected vs. avoided
# and 2) help which choosing a reference category (for the future)
# round(table(data$landcov_target[data$case==1])/sum(data$case==1),3) # used pay attention here to very small number or 0s - if you find some, those classes should be identified as reference category
# round(table(data$landcov_target[data$case==0])/sum(data$case==0),3)  # available


# Jerod's thoughts on how to choose a reference category:
# 1. Should be a variable you do not care about too much
# 2. should not be a rare landcover on the landscape. Should have at least 5 or 10% in availability
# 3. Best to choose one where the use is pretty similar to the availability (i.e., they use it close to its availability)
# 4. Sometimes its OK to merge a few categories that you don't care about into the reference (e.g,. water, rock, ice, etc.)

# create dummy variables for landcover variable (don't forget, need to drop one in analysis)
# table(is.na(data$landcov_target))  #this ideally should be ALL FALSE. If not, then you should remove the lines with NAs
# data$landcover_target <- data$landcov_target
# dum <- data %>% 
#   dplyr::select(landcov_target) %>% 
#   dummy.data.frame(names="landcov_target") #create dummy variables
# 
# head(dum)
# length(unique(data$landcov_target)) == ncol(dum)   # this should be TRUE
# # now cbind up the new columns
# data <- cbind(data, dum)
# head(data)
# rm(dum)


#--------------------------------------#
# Check Correlation among variables ####
#--------------------------------------#
#start by checking correlation
vars_target <- c("elev_target","trasp_target","slope_target","tpi_target",   # for the target point variables
               "rank_perctreecov_target", "rank_forbgrassbiomass_target",
               "rank_IRGVals_target", "rank_csumNDVImax_target")
correl <- data %>% 
  dplyr::select(all_of(vars_target)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)
ifelse(abs(correl)>.5, correl, NA)  # a cleaner way to look at it

# elev and NDVI max correl, and forbgrass biomass and perc tree cov negatively correlated
#tpi and trasp correlated, left out tpi from models

vars_step <- c(  "elev_step","trasp_step","slope_step",  # for the step variables
                 "rank_perctreecov_step", "rank_forbgrassbiomass_step")
correl <- data %>% 
  dplyr::select(all_of(vars_step)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)
ifelse(abs(correl)>.5, correl, NA)  # a cleaner way to look at it

# I cannot have slope_step and Dist2Escape_step in same model. I'll go with slope_step
# I cannot have PerennialForbGrass_step and treecov_step in same model. I'll go with PerennialForbGrass_step

# update variable objects less the ones you deemed too correlated with others
vars_target <- c("elev_target","trasp_target","slope_target",   # for the target point variables
                 "rank_perctreecov_target", "rank_forbgrassbiomass_target",
                 "rank_IRGVals_target", "rank_csumNDVImax_target")
vars_step <- c(  "elev_step","trasp_step","slope_step",   # for the step variables
                 "rank_perctreecov_step", "rank_forbgrassbiomass_step")

# now concatinate all the variables together so you have all the vars you'll need in one object
vars_all <- c(vars_target, vars_step)


# ----------------------------------#
# review variables distributions ####
# ----------------------------------#
# data$case_fact <- ifelse(data$case == 1, "Used","Available")
# for(i in 1:length(vars_all)){   #once you run this, use the back button on your plots tab to take a look at the distributions
#   # Use semi-transparent fill
#   print(ggplot(data, aes_string(x=vars_all[i], fill="case_fact")) +
#     geom_density(alpha=0.4))
# }  #push backwards on the graph window to have a look


library(ggridges)

data$case_fact <- ifelse(data$case == 1, "Used","Available")
for(i in 1:length(vars_all)){   #once you run this, use the back button on your plots tab to take a look at the distributions
  # Use semi-transparent fill
  print(ggplot(data, aes_string(x=vars_all[i], y = "species", fill="case_fact")) +
          geom_density_ridges(alpha=0.4))
}  


# -----------------------------------#
# remove strata with missing data ####
# -----------------------------------#
#view data table of nas
table(is.na(data))
apply(data, 2, function(x){table(is.na(x))})
#remove lines with even 1 NA. You want this, so you have same sample size for both target and step models
nrow(data)
data <- data[apply(data[,c("case",vars_all)],1,anyNA)==FALSE,]
nrow(data)
rm(vars_all)

#remove strata where there is no used case
tokeep <- data %>% 
  group_by(strata) %>% 
  dplyr::summarise(sums=sum(case)) %>% 
  filter(sums == 1)

head(tokeep)

data <- data %>% 
  filter(strata %in% tokeep$strata)
rm(tokeep)

#how many controls do we have for each case (OK to be missing a few here and there, but shouldn't have too many with just a few)
tbl <- data %>% 
  group_by(strata) %>% 
  dplyr::summarize(lengths=length(case))
table(tbl$lengths)

# how to keep only the strata with > X number of case/controls
nrow(data)
min_numb_in_strata <- 5 # what the is is minimum number of case/controls you want per strata?
data <- data %>% 
  filter(strata %in% tbl$strata[tbl$lengths >= min_numb_in_strata])
nrow(data)
table(as.numeric(tapply(data$case,data$strata,length)))   # this provides the distribution of the number of case/controls in each strata
rm(tbl, min_numb_in_strata)

#------------------------------#
# Fit models! --------------####
#------------------------------#
data <- data %>% 
  arrange(id_yr_seas, strata, case)  # this should be ordering by id, then by time

#paramterize with conditional logistic regression - specifying cluster will envoke robust SE estimates using General estimating equations
# Add dist and log(dist) as per Forrester et al. 2009 (Ecology) and Avgar et al. 2016 (MEE)

#point (target models, vegetation quality)
#loop by species
mpoint_list <- list()
spp <- unique(data$species)


for(i in 1:length(spp)){
  data_spp <- data %>% filter(species == spp[i])
  
mpoint <- clogit(case~dist+log(dist)+elev_target+trasp_target+slope_target+
                   rank_perctreecov_target + rank_forbgrassbiomass_target +
                   rank_IRGVals_target + rank_csumNDVImax_target +
                   strata(strata)+ cluster(id_yr_seas), # cluster term says that you assume steps within each cluster are correlated with each other, but steps among the clusters are independent
                 x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
                 method = "efron",data=data_spp)

mpoint_list[[i]] <- mpoint
}

names(mpoint_list) <- spp

summary(mpoint_list[[1]])
confint(mpoint_list[[1]])

#create table of coeff and summary stats, save as txt file
for(i in 1:length(mpoint_list)){
  summ_stat <- as.data.frame(round(coef(summary(mpoint_list[[i]])), digits = 2))
  summ_stat <- tibble::rownames_to_column(summ_stat, "Variable")
  write.table(summ_stat, file = paste0("./Code output/ssf_summ_point_vegqual_movelab5_", spp[[i]], ".txt"),
              sep = ",", quote = FALSE, row.names = FALSE)
}

  ##plot coefficients
  library(sjPlot)
  for(i in 1:length(mpoint_list)){
   print(sjPlot::plot_model(mpoint_list[[i]], 
                     title = paste0(spp[i], " endpoints model coeff.")))
  }

# test variance inflation factors (they should all be below 3-4 or so)
vif(glm(case~dist+log(dist)+elev_target+trasp_target+slope_target+
          rank_perctreecov_target + rank_forbgrassbiomass_target, data=data))

#point (target models, NO vegetation quality)
#loop by species
mpoint_list_nv <- list()
spp <- unique(data$species)


for(i in 1:length(spp)){
  data_spp <- data %>% filter(species == spp[i])
  
  mpoint <- clogit(case~dist+log(dist)+elev_target+trasp_target+slope_target+
                     rank_perctreecov_target + rank_forbgrassbiomass_target +
                     strata(strata)+ cluster(id_yr_seas), # cluster term says that you assume steps within each cluster are correlated with each other, but steps among the clusters are independent
                   x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
                   method = "efron",data=data_spp)
  
  mpoint_list_nv[[i]] <- mpoint
}

names(mpoint_list_nv) <- spp

summary(mpoint_list_nv[[1]])
confint(mpoint_list_nv[[1]])



##plot coefficients
library(sjPlot)
for(i in 1:length(mpoint_list_nv)){
  print(sjPlot::plot_model(mpoint_list_nv[[i]], 
                           title = paste0(spp[i], " points model coeff.")))
}

# test variance inflation factors (they should all be below 3-4 or so)
vif(glm(case~dist+log(dist)+elev_target+trasp_target+slope_target+
          rank_perctreecov_target + rank_forbgrassbiomass_target, data=data))


#create table of coeff and summary stats, save as txt file
for(i in 1:length(mpoint_list_nv)){
  summ_stat <- as.data.frame(round(coef(summary(mpoint_list_nv[[i]])), digits = 2))
  summ_stat <- tibble::rownames_to_column(summ_stat, "Variable")
  write.table(summ_stat, file = paste0("./Code output/ssf_summ_point_movelab5_", spp[[i]], ".txt"),
              sep = ",", quote = FALSE, row.names = FALSE)
}



#step models
mstep_list <- list()

for(i in 1:length(spp)){
  data_spp <- data %>% filter(species == spp[i])
mstep <- clogit(case~dist+log(dist)+elev_step+trasp_step+slope_step+
                  rank_perctreecov_step + rank_forbgrassbiomass_step +
                  strata(strata)+cluster(id_yr_seas),
                x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
                method = "efron",data=data_spp)

mstep_list[[i]] <- mstep
}

names(mstep_list) <- spp
 
summary(mstep_list[[1]])

##plot coefficients
library(sjPlot)
for(i in 1:length(mstep_list)){
  print(sjPlot::plot_model(mstep_list[[i]], 
                           title = paste0(spp[i], " step length model coeff.")))
}
Q
# test variance inflation factors
vif(glm(case~dist+log(dist)+elev_step+trasp_step+slope_step+
          tpi_step+rank_perctreecov_step + rank_forbgrassbiomass_step,
        data=data))


round(summary(mpoint)$coefficients,3)
round(summary(mstep)$coefficients,3)

#create table of coeff and summary stats, save as txt file
for(i in 1:length(mstep_list)){
  summ_stat <- as.data.frame(round(coef(summary(mstep_list[[i]])), digits = 2))
  summ_stat <- tibble::rownames_to_column(summ_stat, "Variable")
  write.table(summ_stat, file = paste0("./Code output/ssf_summ_step_movelab5_", spp[[i]], ".txt"),
              sep = ",", quote = FALSE, row.names = FALSE)
}


#test QIC between the three models
qic_list <- list()
for(i in 1:length(spp)){
qic <- cbind(model=c("point- forage quality", "point", "step"),rbind(CalcQIC(mpoint_list[[i]], details=TRUE), 
                                    CalcQIC(mpoint_list_nv[[i]], details=TRUE),   
                                    CalcQIC(mstep_list[[i]], details=TRUE)))

qic_list[[i]] <- qic
}

names(qic_list) <- spp

#save as txt file
for(i in 1:length(qic_list)){
write.table(qic_list[[i]], file = paste0("./Code output/ssf_qic_movelab5_", spp[[i]], ".txt"),
            sep = ",", quote = FALSE, row.names = FALSE)
}

# QICR is what you want to minimize. It's like AIC, but for models with a cluster term specified
# Check and make sure n and nevent are the exact same between the different models you are comparing!

