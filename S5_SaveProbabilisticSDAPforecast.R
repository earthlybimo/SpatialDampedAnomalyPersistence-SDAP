#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Script to take pre-saved deterministic/SAP forecast file and anomaly weight files and use these to create the probabilistic Spatial Damped Anmaly Forecasts (SDAP)

#Again, syst args are 1. Hemisphere, 2. Initialisation month and 3. (optional) forecast year (or all years if left blank)

ccc=as.integer(args[1]) #Choice of Hemisphere
if(ccc == 1){
  HEM="nh"
  Alphasetname="NewTraining_2021NH_yrs89to98"} #To test with diff training sets, simply change this part appropriately

if(ccc== 2){
  HEM="sh"
  Alphasetname="NewTraining_2021SH_yrs89to98"}

mc=as.integer(args[2]) #Choice of initialisation date (each doy corresponds to the start of a particular month)
yodlist=c(1,32,61,92,122,153,183,214,245,275,306,336)
yodi= yodlist[mc]

fcstylist=1999:2024  #Here, 1999 is used as the lower limit SIMPLY because we know that the training is 1989 to 1998.. if we use a diff training years set, can change this as well
if(length(args)==3) {
  temp=as.integer(args[3])
  if((temp>=1989)&(temp<=2022))  fcstylist=temp #Will be just a single year now
  remove(temp)}

require(tictoc)

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
  ll=min((dim(sipclimHERE)[2]),(dim(siobsbinHERE)[2]))
  SPS=array(NA,dim=ll)
  for (j in 1:ll){
    if(all(is.na(sipclimHERE[,j]))) {SPS[j]=NA
    } else {
      SPS[j]=SPS_eachtimestep(cell_area,sipclimHERE[,j],siobsbinHERE[,j])  }
  }
  return (SPS)}


dayy=format(as.Date(yodi,origin="2015-12-31"),"%d") 
monthh=format(as.Date(yodi,origin="2015-12-31"),"%m")

Alphafile=sprintf("%s/Outputs/Alpha/%s_yod%03i",HEMPATH,Alphasetname,yodi)
if(!file.exists(Alphafile)) stop(paste0("Did not file with anomaly weights: ",basename(Alphafile)))

load(Alphafile);lat=grd$lat
## Now let's use the pre-saved alpha to make our new forecast
SPSdiffmean_normalised=colMeans(SPSdiff_normalised,na.rm = TRUE)
maxcoll=max.col(t(SPSdiffmean_normalised))
bestalpha=t1[maxcoll]  #Results in a vector of anomaly weights, one for each leadtime 

for(ynum in 1:length(fcstylist)){
  ICdate=sprintf("%d%s%s",fcstylist[ynum],monthh,dayy)
  savename=sprintf("%s/Outputs/Forecasts/DampedForecast_%s_%s_yod%03d",HEMPATH,ICdate,Alphasetname,yod)
  if(file.exists(savename)) next()  #Already saved. 
  
  # Load SAP/deterministic forecast
  detFcfile=Sys.glob(sprintf("%s/Outputs/Forecasts/determinForecastfor%s_yod%03d",HEMPATH,ICdate,yodi))  #This was the format how Det/SAPFcsts were saved
  if(length(detFcfile)==0) next()
  load(file = detFcfile,envir = (detFcenv=new.env()))
  
  # Apply damping
  Forecast_wDamping=(detFcenv$Forecast)*0
  ICbinaryall=binarise(detFcenv$AllInitialArr,0.15);    
  for (leadtime in 1:366){  
    alpha=bestalpha[leadtime]
    tempcast=detFcenv$Forecast[,leadtime]*(alpha)  #alpha =0 when SAP weight in the Forecast = 0
    tempcast2=detFcenv$AllClimaArr[,leadtime]*(1-alpha)
    Forecast_wDamping[,leadtime]=tempcast+tempcast2
    remove(tempcast,tempcast2)
  }
  
  
  # Compute the SPS  --------
  SPSForecastwdamping=SPS_clim(grd$Nodeareas,Forecast_wDamping,ICbinaryall)
  
  ## Save:
  savename=sprintf("%s/Outputs/Forecasts/DampedForecast_%s_%s_yod%03d",HEMPATH,ICdate,Alphasetname,yod)
  save(detFcenv,SPSForecastwdamping,Forecast_wDamping, file = savename,version =2)
  #In these files, the original SAP forecasts and SPS (SAP, CLIM, PERS) are all saved inside the 'detFcenv' environment.
  print(paste("Saved file: ",basename(savename),sep = ""))
  remove(detFcenv,SPSForecastwdamping,Forecast_wDamping);invisible(gc())
}  #End of the year/detFcst
