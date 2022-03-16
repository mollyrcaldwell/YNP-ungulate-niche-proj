
# function to caluclate snow on date, i.e., the first day  in fall where there are 
# a specified number of consecutive snow cover days
# Jerod Merkle, Nov 2021

CalcSnowOn <- function(snow.dat=snow.dat,  # ordered vector (of length 365 or 366) representing snow depths
                       depth.thresh=0.02,  # snow depth threshold when below this value there is no snow on the ground
                       na.value = 999,   #na values set to 999, filter these out below
                       day.thresh=7        # how many consecutive dates of snow cover to define as snow on
){
  if(length(snow.dat)>366)
    stop("Your snow.dat has more than 1 year of data in it!")
  if(any(is.na(snow.dat))==TRUE)
    stop("You have NAs in your snow.dat vector.")
  x <- ifelse((snow.dat < depth.thresh & snow.dat != na.value), 0, 1)  # turn into snow cover
  deltas <- c(0,diff(x))  # calculate transitions or deltas over time
  switches <- which(deltas==1)  # find the days that snow started to cover the ground!
  switch.full.snow <- do.call(c, lapply(switches, function(z){
    sum(x[z:(z+day.thresh)])  # calculate how many consecutive days of snow are after a switch
  }))
  # this is first julian day there were day.thresh consecutive snow cover days afterwards
  toreturn <- switches[switch.full.snow == day.thresh]
  toreturn <- toreturn[toreturn>275][1]    # first one after 1 Oct
  return(ifelse(is.na(toreturn), 365, toreturn))  
}
