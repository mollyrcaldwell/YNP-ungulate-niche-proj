#--------------------------------------------#
#- Resource Selection Functions ----------####
#------------- Lab 4 ------------------------#
#-- Part d - Model Selection and validation -#
#---------- Jerod Merkle --------------------#

library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(stringr)
library(move)
library(dummies)
library(lme4)
library(parallel)
library(snow)
library(ggplot2)
library(MuMIn)
library(MASS)
library(car)
library(pROC)
library(visreg)
library(glmmTMB)

#library(devtools)
#install_github("jmerkle1/MerkleLab-Code-Repository", subdir="MoveTools", force=TRUE)
library(MoveTools)

#set your working drive
setwd("~/UWyo/PhD project/YNP-ungulate-niche-proj/")

#--------------------------------------------#
# Load RSF data with variables extracted  ####
#--------------------------------------------#
# from lab RSF part b
load("RSFs_VariablesExtracted.RData")

# -----------------------------------#
# dealing with landcover variable ####
# -----------------------------------#

# some basic selection ratios of your landcover (categorical variable)
# this will help you 1) get a feel for what is selected vs. avoided
# and 2) help which choosing a reference category (for the future)
# round(table(data_pop$landcov[data_pop$used==1])/sum(data_pop$used==1),3) # used pay attention here to very small number or 0s - if you find some, those classes should be identified as reference category
# round(table(data_pop$landcov[data_pop$used==0])/sum(data_pop$used==0),3)  # available

# Jerod's thoughts on how to choose a reference category:
# 1. Should be a variable you do not care about too much
# 2. Should not be a rare landcover on the landscape. Should have at least ~ 5 or 10% in availability
# 3. Best to choose one where the use is pretty similar to the availability (i.e., they use it close to its availability). This helps with interpreting the coefficients of the other land cover types
# 4. Sometimes it's OK to merge a few categories that you don't care about into the reference (e.g,. water, rock, ice, etc.)

# Once you've made some decisions, here is some code to organize your reference category
# Doing the following is also necessary to use some of the automated model selection methods

# create a new landcover factor column, with the correct landcover as reference
# data_pop$landcov2 <- as.character(data_pop$landcov)  # note that I am changing it to a character vector
# table(data_pop$landcov2)
# # reclassify one at a time (renaming the ones you want in the reference category to reference)
# data_pop$landcov2[data_pop$landcov2=="developed"] <- "reference"
# data_pop$landcov2[data_pop$landcov2=="water"] <- "reference"
# data_pop$landcov2[data_pop$landcov2=="wetlands"] <- "reference"
# table(data_pop$landcov2)
# 
# # now specify your new lanccover variable as a factor
# data_pop$landcov2 <- as.factor(data_pop$landcov2)
# levels(data_pop$landcov2)   #these are the current levels. they are in alphabetical order, but you want reference to be the first one!
# data_pop$landcov2 <- relevel(data_pop$landcov2, ref=3) #the value of ref is where in the string of levels that your reference is located. Mine is 3rd in order
# levels(data_pop$landcov2)   #check and make sure it is right (your reference category must be the first level)

# --------------------------------#
# Model Selection -  prep data ####
# --------------------------------#

#convert forbgrass biomass to percentage
data_comm <- data_comm %>%
  mutate(perc_forbgrass_biomass = (forbgrass_biomass/max(forbgrass_biomass))*100)

# deal with NAs...
# How many NAs are in each column?
apply(data_comm, 2, function(x){table(is.na(x))})

#remove lines where there is an at least 1 NA in your columns
variables <- c("elev","trasp","slope", "forbgrass_biomass", "perc_treecov")
nrow(data_comm)
data_comm <- data_comm[apply(data_comm[,variables],1,anyNA)==FALSE,]
nrow(data_comm)   # Hopefully you aren't deleting too many. 

# I suggest reducing your database so you can go through the motions of this section a bit more quickly
#data_comm_red <- sample_n(data_comm, 4000)  # Randomly sample to about 4k points. Models seem to fit pretty quickly with this.
#sample 3 individuals per species
spp <- unique(data_comm$species)
spp_sample <- list()
for(i in 1:length(spp)){
  data_spp <- data_comm %>% filter(species == spp[i])
  samp <- sample(unique(data_spp$id_yr_seas), 2)
  
  spp_sample[[i]] <- data_comm %>% filter(id_yr_seas %in% samp)
}
  
data_comm_red <- rbind(spp_sample[[1]], spp_sample[[2]], spp_sample[[3]],
                       spp_sample[[4]], spp_sample[[5]])
  
#-------------------------------------#
# Manual forward stepwise approach ####
#-------------------------------------#

# build a full model here with all of your variables including your new landcover2 factor column 
# OK to include ones that are correlated with each other...
# Note that you will get convergence issues because of so many correlated variables!!!!! That's OK.
# but you should not interpret the coefficients of this FULL model. it is just to have all the variables in one place

full_comm <- glmer(used ~ elev + trasp + slope + perc_treecov + perc_forbgrass_biomass +
                    (1|species/id_yr_seas), # For this exercise, only use random intercept (things'll go faster)
                  family = binomial(link = "logit"), data = data_comm_red, na.action="na.fail")
summary(full_comm)


full_comm_ri <- glmmTMB(used ~ elev + trasp + slope + perc_treecov + perc_forbgrass_biomass +
                     (1|species/id_yr_seas) + (0+elev|species/id_yr_seas) + 
                       (0+trasp|species/id_yr_seas) + (0+slope|species/id_yr_seas) +  
                     (0+perc_treecov|species/id_yr_seas) +
                     (0+perc_forbgrass_biomass|species/id_yr_seas), 
                   family = binomial(), data = data_comm_red, na.action="na.fail")
 
summary(full_comm_ri)

# The stepwise model selection starts here
# now, fit a null model, with only an intercept
# null_comm <- glmer(used ~ 1 + (1|species/id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# 
# #try adding each variable in the full model, one by one to the null model
# addterm(null_comm, full_comm, sorted=TRUE, test="Chisq")
# 
# #if there is a variable with support, add it to the next model
# m1_comm <- glmer(used ~ elev + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# #check summary, look at the sign and significance of the new variable
# summary(m1_comm)
# 
# #now check to see if any new terms provide any more empirical support
# addterm(m1_comm, full_comm, test="Chisq", sorted=TRUE)
# #add in the best supported new term
# m2_comm <- glmer(used ~ landcov2 + slope + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# 
# # check vifs with new term (make sure the new variable isn't correlated with another variable). Don't want VIFs > 4... 
# # If it is correlated, take it out and move on to the next most important variable in the addterm() result
# vif(glm(used ~ landcov2 + slope, family = binomial(link = "logit"), data = data_comm_red))
# # if VIFs are OK, check the parameters, have any of the signs flipped from the previous model? Are they similar? Are they all still significant?
# summary(m2_comm)
# 
# # add a new term
# addterm(m2_comm, full_comm, test="Chisq", sorted=TRUE)
# m3_comm <- glmer(used ~ landcov2 + slope + PerennialForbGrass + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# vif(glm(used ~ landcov2 + slope + PerennialForbGrass, family = binomial(link = "logit"), data = data_comm_red))
# summary(m3_comm)
# 
# addterm(m3_comm, full_comm, test="Chisq", sorted=TRUE)
# m4_comm <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi, family = binomial(link = "logit"), data = data_comm_red))
# summary(m4_comm)
# 
# addterm(m4_comm, full_comm, test="Chisq", sorted=TRUE)
# m5_comm <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + Dist2roads + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + Dist2roads, family = binomial(link = "logit"), data = data_comm_red))
# summary(m5_comm)
# 
# addterm(m5_comm, full_comm, test="Chisq", sorted=TRUE)
# m6_comm <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + Dist2roads + trasp +  + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + Dist2roads + trasp , family = binomial(link = "logit"), data = data_comm_red))
# summary(m6_comm)
# 
# addterm(m6_comm, full_comm, test="Chisq", sorted=TRUE)
# m7_comm <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev, family = binomial(link = "logit"), data = data_comm_red))
# summary(m7_comm)
# 
# addterm(m7_comm, full_comm, test="Chisq", sorted=TRUE)
# m8_comm <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + treecov + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + treecov, family = binomial(link = "logit"), data = data_comm_red))
# # my treecov is pretty correlated with landcov2, gonna use next variable in line with addterm()
# m8_comm <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround, family = binomial(link = "logit"), data = data_comm_red))
# summary(m8_comm)
# 
# addterm(m8_comm, full_comm, test="Chisq", sorted=TRUE)
# m9_comm <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + iNDVI19 + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + iNDVI19, family = binomial(link = "logit"), data = data_comm_red))
# summary(m9_comm)
# 
# addterm(m9_comm, full_comm, test="Chisq", sorted=TRUE)
# m10_comm <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + (1|id_yr_seas), family = binomial(link = "logit"), data = data_comm_red)
# vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + iNDVI19 + AnnualForbGrass, family = binomial(link = "logit"), data = data_comm_red))
# summary(m10_comm)
# addterm(m10_comm, full_comm, test="Chisq", sorted=TRUE)
# 
# # stop here, last variable does not improve AIC all that much... 
# # Or, Likelihood ratio test not significant and/or other variables created too high of vifs
# 
# best_comm <- m10_comm    # this is the best model based on step-wise model selection
# summary(best_comm)

# ----------------------------------------#
# Hypothesis-based AIC model selection ####
# ----------------------------------------#

# True AIC model selection is about testing relative empirical support for competing hypotheses formulated as models
# For example, I could assess empirical support for models containing only topographic variables vs. topo+habitat vs. topo+habitat+human

topo <- glmmTMB(used ~ elev + trasp + slope +
                  (1|species/id_yr_seas) + (0+elev|species/id_yr_seas) + 
                  (0+trasp|species/id_yr_seas) + (0+slope|species/id_yr_seas), 
                family = binomial(), data = data_comm_red, na.action="na.fail")

# topo_habitat <- glmmTMB(used ~ elev + trasp + slope + perc_treecov + perc_forbgrass_biomass +
#                           (1|species/id_yr_seas) + (0+elev|species/id_yr_seas) + 
#                           (0+trasp|species/id_yr_seas) + (0+slope|species/id_yr_seas) +  
#                           (0+perc_treecov|species/id_yr_seas) +
#                           (0+perc_forbgrass_biomass|species/id_yr_seas), 
#                         family = binomial(), data = data_comm_red, na.action="na.fail")

topo_habitat <- full_comm_ri

null <- glmmTMB(used ~ 1 + (1|species/id_yr_seas), 
                family = binomial(), data = data_comm_red, na.action="na.fail")


AIC_tble <- model.sel(topo, topo_habitat, null, rank="AIC")  # make AIC table for these three models
AIC_tble

# --------------------------#
# Dredge function - uses ####
# --------------------------#
?dredge

# A note on the dredge function. It is not wise to use dredge function to try all combinations of your variables.
# This is dredging and will increase the chances of a type I error, where you find some relationship but it is not there.
# The idea is simply that trying so many combinations of variables the chances increase that you 'find' a relationship due to random chance.


# That being said, the dredge function has many arguments that can help you with your model selection needs.
# let's have a look:

#start by creating a correlation matrix
variables <- c("elev","slope","trasp","perc_treecov", "perc_forbgrass_biomass")
correl <- data_comm_red %>% 
  dplyr::select(all_of(variables)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)
ifelse(abs(correl)>.5, correl, NA)  # a cleaner way to look at it

# elev and forb grass biomass correlated, tree cover and forb grass biomass correlated

# build a matrix that provides information on which variables should and should not be included in the same model
# smat <- abs(correl) <= .5    #lets use 0.5 as the cut-off for correlation, where variables that are correlated larger than this number CANNOT be included in same model
# smat[!lower.tri(smat)] <- NA #Make everything in the matrix NA except (!lower.tri)the lower triangle
# smat   # when FALSE, the two variables cannot be in same model
# i <- as.vector(smat == FALSE & !is.na(smat))
# sexpr <- parse(text = paste("!(", paste("(",variables[col(smat)[i]], " && ",variables[row(smat)[i]], ")",sep = "", collapse = " || "), ")"))
# sexpr   # this tells which variables cannot be in same model. Check help with dredge to see all the different ways you can specify things here
# rm(smat, i, correl)


# build a full model here with all of your variables including your new landcover2 factor column 
# OK to include ones that are correlated with each other...
# Note that you will get convergence issues because of so many correlated variables!!!!! That's OK.
# but you should not interpret the coefficients of this FULL model. it is just to have all the variables in one place

full_comm <- topo_habitat
summary(full_comm)


# dredge tests combinations of variables from the variables in a given FULL model, 
# based on your smat, and how many variables you want to keep in each model.
# Note that dredge can also help you with interaction terms and polynomial terms too (not shown below, though)
# we will run this using parrallel processing

#next three lines set up the parralell processing
cores <- detectCores()-1
clust <- makeCluster(getOption("cl.cores", cores), type = "SOCK")
invisible(clusterEvalQ(clust, library(glmmTMB)))   # you need to export lme4 to each core
clusterExport(clust, "data_comm_red")    # you need to export your database to each core

dredge_comm <- pdredge(full_comm, beta="none", rank="AIC", cluster=clust,
                      #subset=sexpr, #this is where you specify which variables can and cannot be in the same model together
                      #fixed=c("slope", "elev"), # perhaps you want to specify variables that MUST be in each model?
                      m.lim=c(4,NA)) # tried to speed up the method by setting m.min to 10 (i.e., it must have 10 variables in each model). The second number is the upper level. Make this number the max variables that you have. NA means max.
stopCluster(clust)   # end parrallel processing
head(dredge_comm)    #these are the top models
nrow(dredge_comm)     # this is how many models it tested
dredge_comm    #these are the top models

top_mod <- glmmTMB(used ~ elev + trasp + perc_treecov + perc_forbgrass_biomass +
                     (1|species/id_yr_seas) + (0+elev|species/id_yr_seas) + 
                     (0+trasp|species/id_yr_seas) +  
                     (0+perc_treecov|species/id_yr_seas) +
                     (0+perc_forbgrass_biomass|species/id_yr_seas), 
                   family = binomial(), data = data_comm_red, na.action="na.fail")
summary(top_mod)

# -----------------------------------------#
# Plotting predicted responses, quickly ####
# -----------------------------------------#

# can fairly quickly plot responses by variable with visreg package
visreg(top_mod, xvar = "elev", type="conditional", scale="linear", 
       trans=exp, ylab="Relative Probability of Use")   # plot on linear scale (before transforming back to binomial world)

visreg(top_mod, xvar = "trasp", type="conditional", scale="linear", 
       trans=exp, ylab="Relative Probability of Use")   # plot on linear scale (before transforming back to binomial world)

visreg(top_mod, xvar = "perc_treecov", type="conditional", scale="linear", 
       trans=exp, ylab="Relative Probability of Use")   # plot on linear scale (before transforming back to binomial world)

visreg(top_mod, xvar = "perc_forbgrass_biomass", type="conditional", scale="linear", 
       trans=exp, ylab="Relative Probability of Use")   # plot on linear scale (before transforming back to binomial world)

visreg(top_mod, xvar = "elev", by="trasp", type="conditional", scale="linear", 
       trans=exp, ylab="Relative Probability of Use")   # plot on linear scale (before transforming back to binomial world)


#------------------------#
# calculate R squared ####
#------------------------#
#calculation of R^2 for mixed models based on Nakagawa and Schielzeth 2013, adn Johnson 2014.
r.squaredGLMM(topo_habitat)

# do AUC/ROC for logistic model (Note: because there is contamination in an RSF, Boyce recommends NOT doing AUC, but using non-parametric k-folds)
boxplot(predict(top_mod, type="response")~topo_habitat[["frame"]]$used, xlab="used(1) versus available(0)",
        ylab="predicted probability of use")
rocauc <- roc(topo_habitat[["frame"]]$used, predict(topo_habitat, type="response"), ci=TRUE)
rocauc   #if it is near 0.5, not a good model!
plot(rocauc)
# sensitivity is true-positive rate
# 1-specificity is the false positive rate


#----------------------------#
# kfolds cross validation ####
#----------------------------#

#Initial code to understand out how predicted values work

predvals_logit <- predict(top_mod, type="response")   #get the predicted values on the logistic scale
hist(predvals_logit) # so this is the predicted values on the logit scale, i.e., exp(model)/(1+exp(model))
hist(logit(predvals_logit))  #use the logit function to get linear predictions

predvals_linear <- predict(top_mod, type="link")   #another way to get linear predictions
hist(predvals_linear)

hist(exp(predvals_linear)/(1+exp(predvals_linear)))   #now we can get the predvals_logit by doing a logistic transformation of the linear predictors

hist(exp(predvals_linear))    # the exponent of the linear predictors
RSF_predictions <- exp(predvals_linear)   #this is the regular RSF formation (although we did not drop the intercept, but not a big deal because everything is relative)
RSF_predictions <- (RSF_predictions-min(RSF_predictions))/(max(RSF_predictions)-min(RSF_predictions))   #scale bewteen 0 and 1, based on Manley 
hist(RSF_predictions)


# FIRST, do validation without doing k-folds partitioning (helps understand what is going on)
# get 10 bins of available data using quantiles (i.e., equal number of points per bin)
quant <- quantile(RSF_predictions[top_mod[["frame"]]$used==0], 
                  probs = seq(from = 0,to = 1, length.out = 11))
quant[1] <- 0    #make sure the lower quantile is 0
quant[length(quant)] <- Inf    #make sure the upper quantil is Inf
quant
#this shows that there are an equal number of random points in each bin
table(factor(findInterval(RSF_predictions[top_mod[["frame"]]$used==0],quant), levels = 1:10))

#break up predictions into the bins you specified with your random points
# so we know how many 1s are in each bin
int <- factor(findInterval(RSF_predictions[top_mod[["frame"]]$used==1],quant), levels = 1:10)   
table(int)

#plot it
par(mai=c(.4,.4,.02,.02))
plot(1:10, as.numeric(table(int)), ylab="Number of used values in each bin",
     ylim=c(1,max(as.numeric(table(int)))), xlab="RSF score bin", 
     type="b", cex.axis=0.6,tcl=-.1, mgp=c(1,0,0), cex.lab=.8)
legend("topleft", paste("Spearman Rank Correlation = ", round(cor(1:10, table(int), 
                                  method = "spearman"),3), sep=""), bty="n", cex=.7)


# SECOND, run the cross validation by splitting up your data into testing and training partitions
#don't use on the nested species- doesn't make a lot of sense
cv_comm <- CalcKfoldCV(data=data_comm_red,
                      formula="used ~ elev + trasp + perc_treecov + perc_forbgrass_biomass +
                     (1|species/id_yr_seas) + (0+elev|species/id_yr_seas) + 
                     (0+trasp|species/id_yr_seas) +  
                     (0+perc_treecov|species/id_yr_seas) +
                     (0+perc_forbgrass_biomass|species/id_yr_seas)",
                      modType="glmmTMB",
                      binVar="species",
                      k=5,
                      modFam="binomial",
                      resp="used")
cv_comm

mean(cv_comm)  # see Boyce et al. eval RSFs for info on how to interpret
sd(cv_comm)

# note that often folks repeat this CV analysis a number of times to obtain variation inherent in sampling of different folds. 



#graph & summary  topographic/habitat model results
summ_th <- summary(topo_habitat)
confint(topo_habitat)
coef <- coefficients(topo_habitat)  # have a look at the coefficients for each id_yr_seas

coef_df <- coef[["cond"]][["species"]]

sd_df <- data.frame(variable = c("elev", "trasp", "slope", "perc_treecov", "perc_forbgrass_biomass"),
                    sd = c(3.4409e-05, 5.5962e-03, 6.6620e-02, 5.1124e-02, 9.0023e-02))

coef_df <- coef_df %>%
  mutate(CI_elev = 1.96 * (sd_df$sd[1]/sqrt(2))) %>%
  mutate(CI_trasp = 1.96 * (sd_df$sd[2]/sqrt(2))) %>%
  mutate(CI_slope = 1.96 * (sd_df$sd[3]/sqrt(2))) %>%
  mutate(CI_perc_treecov = 1.96 * (sd_df$sd[4]/sqrt(2))) %>%
  mutate(CI_forbgrass_biomass = 1.96 * (sd_df$sd[5]/sqrt(2))) %>%
  mutate(species = row.names(coef_df)) %>%
  dplyr::select(-CI)

coef_df_l <- pivot_longer(coef_df[,c(1:6, 12)], cols = c("(Intercept)", "elev", "trasp", "slope", 
                                            "perc_treecov", "perc_forbgrass_biomass"),
                          names_to = "Variable") %>%
  mutate(Variable = if_else(Variable == "perc_forbgrass_biomass", "forbgrass_biomass", Variable)) %>%
  rename("Coefficient" = value)

ci_df_l <- pivot_longer(coef_df[,c(7:12)], cols = c("CI_elev", "CI_trasp", "CI_slope", "CI_perc_treecov", 
                                                    "CI_forbgrass_biomass"))

ci_df_l <- ci_df_l %>%
  #mutate(Variable = str_replace(name, "CI_", "")) %>%
  #dplyr::select(-name) %>%
  rename("CI" = value)

graph_tbl <- left_join(coef_df_l, ci_df_l, by = c("species", "Variable")) %>%
  filter(Variable != "(Intercept)")

coef_tbl <- coef_df %>%
  mutate(elev = paste0(round(elev, 3), " +/- ", round(CI_elev,3))) %>%
  mutate(trasp = paste0(round(trasp, 3), " +/- ", round(CI_trasp,3))) %>%
  mutate(slope = paste0(round(slope, 3), " +/- ", round(CI_slope,3))) %>%
  mutate(perc_treecov = paste0(round(perc_treecov, 3), " +/- ", round(CI_perc_treecov,3))) %>%
  mutate(perc_forbgrass_biomass = paste0(round(perc_forbgrass_biomass, 3), " +/- ", round(CI_forbgrass_biomass,3))) %>%
  dplyr::select(`(Intercept)`, elev, trasp, slope, perc_treecov, perc_forbgrass_biomass) 

row.names(coef_tbl) <- NULL

write.table(coef_tbl, "./Code output/communityRSF_topohabitat_coefCI.txt",
            sep = ",", quote = FALSE, row.names = FALSE)




#graph
graph_tbl$Variable <- factor(graph_tbl$Variable, levels = c("elev", "trasp",
                                "slope", "perc_treecov", "forbgrass_biomass"))

ggplot(graph_tbl, aes(x = Variable, y = Coefficient)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_point(aes(x = Variable, 
                 y = Coefficient)) + 
  geom_linerange(aes(x = Variable, 
                     ymin = (Coefficient - CI),
                     ymax = (Coefficient + CI)),
                 lwd = 1) +
    facet_wrap(~species) +
  coord_flip()
