#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Script to take pre-saved deterministic/SAP forecast files, compute a range of forecasts with varying anomaly weights and measure the SPS then save it to choose the best SDAP forecast.
#Sysargs are 1. Hemisphere, 2.Initialisation month to test

ccc=as.integer(args[1]) #Choice of Hemisphere
yodchoice=as.integer(args[2])  #Choice of month/initial-doy
yodlist=c(1,32,61,92,122,153,183,214,245,275,306,336)
yod=yodlist[yodchoice]

## VERY important to choose which years we are using for training
ylist=1989:1998

# Also recommended to save the Alphaset (training weights) with a descriptive name in case you wish to test/compare different years/sets, separated by the hemisphere
if(ccc == 1){
  HEM="nh"
  Alphasetname="NewTraining_2021NH_yrs89to98"}

if(ccc== 2){
  HEM="sh"
  Alphasetname="NewTraining_2021SH_yrs89to98"}

require(tictoc)
tic("whole script")

MASTERPATH="~/WORK/Data/SDAP/"
HEMPATH=paste0(MASTERPATH,"/",HEM)

binarise <-function (somearr,dlevel) {  #Function to binarise some given array  based on  this level
  ll=dim(somearr)
  if (is.null(ll)){ll=length(somearr)}  #because 1d arrays dont do 
  ibinar=array(dim =ll)
  ibinar[somearr[]>=dlevel]=1
  ibinar[somearr[]<dlevel]=0
  return (ibinar)
}
SPS_eachtimestep <-function (cell_area,Pfcst,Pobs) {  #Function to find the SPS at each timestep
  SPS = ((Pfcst-Pobs)^2)*cell_area
  return (sum(SPS,na.rm = TRUE))}
SPS_clim <-function (cell_area,sipclimHERE,siobsbinHERE){  
  #Function to compute SPS of climatology, assuming sipclimHERE (Nlat x 366) and siobsbinHERE (Nlat x 366). 
  ll=min((dim(sipclimHERE)[2]),(dim(siobsbinHERE)[2]))
  SPS=array(NA,dim=ll)
  for (j in 1:ll){
    if(all(is.na(sipclimHERE[,j]))) {SPS[j]=NA
    } else {
      SPS[j]=SPS_eachtimestep(cell_area,sipclimHERE[,j],siobsbinHERE[,j])  }
  }
  return (SPS)}



dayy=format(as.Date(yod,origin="2015-12-31"),"%d") 
monthh=format(as.Date(yod,origin="2015-12-31"),"%m")

SPSdampFcst=array(dim = c(length(ylist),20,366))
SPSdiff=array(dim = c(length(ylist),20,366))  # Diff from Climatological
SPSdiff_normalised=array(dim = c(length(ylist),20,366))

#load the file
for(ynum in 1:length(ylist)){
  ICdate=sprintf("%d%s%s",ylist[ynum],monthh,dayy)
  seekname=Sys.glob(sprintf("%s/Outputs/Forecasts/determinForecastfor%s_yod%03d",HEMPATH,ICdate,yod))  #This was the format how DetFcsts were saved
  if(length(seekname)==0) next()
  load(seekname,envir = (seekenv=new.env()))
  
  ICbinaryall=binarise(seekenv$AllInitialArr,0.15);
  
  ### An array of damping:
  alpha=1  #alpha is anomaly weight. If full damping, alpha =0 so SAP part goes to 0
  t1=seq(0,1,length.out = 20)
  for(counter in 1:20){
    alpha=t1[counter]
    Forecast_wDamping=array(dim = c(length(lat),366))
    
    for (tlead in 1:366){  
      tempcast=seekenv$Forecast[,tlead]*(alpha)  #If full damping, alpha =0 so Forecast part goes to 0
      tempcast2=seekenv$AllClimaArr[,tlead]*(1-alpha)
      Forecast_wDamping[,tlead]=tempcast+tempcast2
      remove(tempcast,tempcast2)
    }
    
    SPSForecastwdamping=SPS_clim(seekenv$grd$Nodeareas,Forecast_wDamping,ICbinaryall)
    tempdiff=(seekenv$SPS_list$SPSClimatological- SPSForecastwdamping)
    SPSdiff[ynum,counter,]=tempdiff
    SPSdiff_normalised[ynum,counter,]=(tempdiff)/seekenv$SPS_list$SPSClimatological
    remove(SPSForecastwdamping,Forecast_wDamping,tempdiff)
    invisible(gc())
  }  # counter/alpha loop.
}
SPSdiffmean=colMeans(SPSdiff,na.rm = TRUE)
SPSdiffmean_normalised=colMeans(SPSdiff_normalised,na.rm = TRUE)
grd=seekenv$grd
savename=sprintf("%s/Outputs/Alpha/%s_yod%03i",HEMPATH,Alphasetname,yod)
save(SPSdiff,SPSdiffmean,SPSdiff_normalised,ylist,yod,t1,grd, file = savename,version =2)
print(paste("Saved file: ",basename(savename),sep = ""))
