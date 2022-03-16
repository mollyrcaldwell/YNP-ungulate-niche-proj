##############################################################
# Function to extract spatially and temporally explicit data 
# from a dated raster stacks created from ndviCALC 
# Returns a vector with length of line data with the mean line/date values of NDVImetric
####### 7 March 2022 ###########################################


##########FUNCTION STARTS HERE ######################
ndviEXTRACT_lines <- function(line_data=dataL, NDVImetric="fittedVals", NDVIfolder="E:/ndvi_bisc_calcs",
                        maxcpus=4, geom_name= "geometry", datesname="date", scaleIRG=TRUE,
                        xyCRS=CRS("+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs+towgs84=0,0,0")){
  
  #some checks
  if(inherits(line_data, "data.frame") == FALSE) stop("line_data is not a dataframe")
  if(any(c(geom_name) %in% colnames(line_data) == FALSE) == TRUE) stop("you do not have a geometry column")
   if(inherits(line_data$datesname, "POSIXct") != TRUE) 
    stop("line_data[,datesname] is not POSIXct")
  # if(any(is.na(line_data$datesname) == TRUE)) 
    # stop("You have NAs in your date column")
  drs <- dir(NDVIfolder)
  if(NDVImetric %in% c("maxNDVIdate","maxIRGdate", "springStart","springEnd","csumNDVImax","maxBrownDowndate",
                       "fittedVals","IRGVals","xmidS_SD", "springLength","sumIRG","SpringScale")==FALSE){
    stop("The NDVImetric must only be maxNDVIdate, maxIRGdate, springStart, springEnd, csumNDVImax, fittedVals, IRGVals, maxBrownDowndate, xmidS_SD,SpringScale, springLength, or sumIRG.")
  }
  if(NDVImetric %in% c("maxNDVIdate","maxIRGdate", "springStart","springEnd","csumNDVImax","springLength")){
    if(length(grep(NDVImetric, drs)) == 0)
      stop("You have an NDVImetric name that does not correspond to files in your NDVIfolder")
  }else{
    if(length(grep("dlcparams", drs)) == 0)
      stop("You don't have dlcparams in your NDVIfolder")
  }
  if(all(c("raster","rgdal","snowfall") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: raster, rgdal, and snowfall")
  require("raster")
  require("rgdal")
  require("snowfall")
  require("lubridate")
  
  tz <- attr(data$datesname,"tzone")
  
  data <- line_data
  data$year <- as.numeric(year(data$datesname))
  data$jul <- as.numeric(yday(data$datesname))
  data$jul[data$jul > 365] <- 365
  data$unique <- 1:nrow(data)
  u <- unique(data$year)
  
  cpus <- ifelse(length(u)>maxcpus, maxcpus, length(u))
  
  if(cpus > 1){   #this code in this if statement replicates the code in the else part, just no paralell processing
    sfInit(parallel = T, cpus = cpus)   #must change the cpus
    sfExport("data", "u", "geom_name", "NDVImetric", "NDVIfolder", "xyCRS","scaleIRG")
    sfLibrary(raster)
    sfLibrary(rgdal)
    vals <- do.call(rbind, sfClusterApplyLB(1:length(u), function(i){
      temp <- data[data$year==u[i],]
      temp <- temp[is.na(temp[,geom_name])==FALSE,]
      coordinates(temp) <- geom_name
      proj4string(temp) <- xyCRS
      
      if(NDVImetric %in% c("maxNDVIdate","maxIRGdate","springStart","springEnd","csumNDVImax","springLength")){
        stk <- stack(paste(NDVIfolder, "/", NDVImetric, ".grd", sep="")) 
        stk <- stk[[paste("X",u[i],sep="")]]
        vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
        toreturn <- data.frame(unique=temp$unique, vals=vals)
        names(toreturn) <- c("unique", NDVImetric)
        return(toreturn)
      }else{
        if(NDVImetric %in% c("xmidS_SD","SpringScale","maxBrownDowndate")){
          if(NDVImetric == "xmidS_SD"){
            stk <- stack(paste(NDVIfolder, "/dlcparams", ".grd", sep="")) 
            stk <- stk[[grep(paste("xmidS_SD_", u[i], sep=""),names(stk))]]
            vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }
          if(NDVImetric == "SpringScale"){
            stk <- stack(paste(NDVIfolder, "/dlcparams", ".grd", sep="")) 
            stk <- stk[[grep(paste("scalS_", u[i], sep=""),names(stk))]]
            vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }
          if(NDVImetric == "maxBrownDowndate"){
            stk <- stack(paste(NDVIfolder, "/dlcparams", ".grd", sep="")) 
            stk <- stk[[grep(paste("xmidA_", u[i], sep=""),names(stk))]]
            vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }
        }else{
          stk <- stack(paste(NDVIfolder, "/dlcparams", ".grd", sep="")) 
          stk <- stk[[grep(u[i],names(stk))]]
          stk <- stk[[1:4]] #leave out the SDs
          vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
          # vals <- extract(stk, 30:32)
          time=1:365
          if(NDVImetric == c("fittedVals")){
            vals <- do.call(rbind, lapply(1:nrow(vals), function(e){
              (1/(1+exp((vals[e,1]-time)/vals[e,2])))-(1/(1+exp((vals[e,3]-time)/vals[e,4])))
            }))
            vals <- do.call(c, lapply(1:nrow(vals), function(e){
              return(vals[e,temp$jul[e]])
            }))
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }else{
            vals <- do.call(rbind, lapply(1:nrow(vals), function(e){
              temp <- (1/(1+exp((vals[e,1]-time)/vals[e,2])))
              temp <- c((diff(temp)/diff(time)),NA)
              temp[temp < 0] <- 0
              if(scaleIRG==TRUE){
                if(all(is.na(temp))==FALSE){
                  temp <- (temp-min(temp, na.rm=TRUE))/(max(temp, na.rm=TRUE)-min(temp, na.rm=TRUE))
                }
              }
              if(NDVImetric=="sumIRG"){
                temp[is.na(temp)==TRUE] <- 0
                return(cumsum(temp))
              }else{
                return(temp)
              }
            }))
            vals <- do.call(c, lapply(1:nrow(vals), function(e){
              return(vals[e,temp$jul[e]])
            }))
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }
        }
      }
    }))
    sfStop()
  }else{
    vals <- do.call(rbind, lapply(1:length(u), function(i){
      temp <- data[data$year==u[i],]
      temp <- temp[is.na(temp[,xname])==FALSE,]
      coordinates(temp) <- c(xname, yname)
      proj4string(temp) <- xyCRS
      
      if(NDVImetric %in% c("maxNDVIdate","maxIRGdate","springStart","springEnd","csumNDVImax","springLength")){
        stk <- stack(paste(NDVIfolder, "/", NDVImetric, ".grd", sep="")) 
        stk <- stk[[paste("X",u[i],sep="")]]
        vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
        toreturn <- data.frame(unique=temp$unique, vals=vals)
        names(toreturn) <- c("unique", NDVImetric)
        return(toreturn)
      }else{
        if(NDVImetric %in% c("xmidS_SD","SpringScale","maxBrownDowndate")){
          if(NDVImetric == "xmidS_SD"){
            stk <- stack(paste(NDVIfolder, "/dlcparams", ".grd", sep="")) 
            stk <- stk[[grep(paste("xmidS_SD_", u[i], sep=""),names(stk))]]
            vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }
          if(NDVImetric == "SpringScale"){
            stk <- stack(paste(NDVIfolder, "/dlcparams", ".grd", sep="")) 
            stk <- stk[[grep(paste("scalS_", u[i], sep=""),names(stk))]]
            vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }
          if(NDVImetric == "maxBrownDowndate"){
            stk <- stack(paste(NDVIfolder, "/dlcparams", ".grd", sep="")) 
            stk <- stk[[grep(paste("xmidA_", u[i], sep=""),names(stk))]]
            vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }
        }else{
          stk <- stack(paste(NDVIfolder, "/dlcparams", ".grd", sep="")) 
          stk <- stk[[grep(u[i],names(stk))]]
          stk <- stk[[1:4]] #leave out the SDs
          vals <- extract(stk, spTransform(temp, CRS(projection(stk))), fun = mean)
          # vals <- extract(stk, 30:32)
          time=1:365
          if(NDVImetric == c("fittedVals")){
            vals <- do.call(rbind, lapply(1:nrow(vals), function(e){
              (1/(1+exp((vals[e,1]-time)/vals[e,2])))-(1/(1+exp((vals[e,3]-time)/vals[e,4])))
            }))
            vals <- do.call(c, lapply(1:nrow(vals), function(e){
              return(vals[e,temp$jul[e]])
            }))
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }else{
            vals <- do.call(rbind, lapply(1:nrow(vals), function(e){
              temp <- (1/(1+exp((vals[e,1]-time)/vals[e,2])))
              temp <- c((diff(temp)/diff(time)),NA)
              temp[temp < 0] <- 0
              if(scaleIRG==TRUE){
                if(all(is.na(temp))==FALSE){
                  temp <- (temp-min(temp, na.rm=TRUE))/(max(temp, na.rm=TRUE)-min(temp, na.rm=TRUE))
                }
              }
              if(NDVImetric=="sumIRG"){
                temp[is.na(temp)==TRUE] <- 0
                return(cumsum(temp))
              }else{
                return(temp)
              }
            }))
            vals <- do.call(c, lapply(1:nrow(vals), function(e){
              return(vals[e,temp$jul[e]])
            }))
            toreturn <- data.frame(unique=temp$unique, vals=vals)
            names(toreturn) <- c("unique", NDVImetric)
            return(toreturn)
          }
        }
      }
    }))
  }
  gc()
  data <- merge(data[,c("unique", "date")],vals, all.x=TRUE)
  data <- data[order(data$unique),]
  return(data[,NDVImetric])
}############### FUNCTION ENDS HERE ###########################