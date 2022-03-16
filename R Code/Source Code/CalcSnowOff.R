
# function to caluclate snowoff date, i.e., the first day there are 
# a specified number of snow free days in a row
# Jerod Merkle, Nov 2021


CalcSnowOff <- function(snow.dat=snow.dat,  # ordered vector (of length 365 or 366) representing snow depths
                        depth.thresh=0.02,  # snow depth threshold when below this value there is no snow on the ground
                        na.value = 999,   #na values set to 999, filter these out below
                        day.thresh=7        # how many consecutive dates of no snow cover to define as snow off
){
  if(length(snow.dat)>366)
    stop("Your snow.dat has more than 1 year of data in it!")
  if(any(is.na(snow.dat))==TRUE)
    stop("You have NAs in your snow.dat vector.")
  x <- ifelse((snow.dat < depth.thresh & snow.dat != na.value), 0, 1)  # turn into snow cover
  deltas <- c(-1,diff(x))   # calculate transitions or deltas over time
  switches <- which(deltas==-1)  # find the locations were snow melted (went from a 1 to a 0 in two consecutive days)
  switch.no.snow <- do.call(c, lapply(switches, function(z){
    sum(x[z:(z+day.thresh)])  # calculate how many consecutive days of snow free are after a switch
  }))
  # this is first julian day there were day.thresh consecutive snow free days afterwards
  return(switches[switch.no.snow == 0][1])  
}
