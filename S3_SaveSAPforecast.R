#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Script to take the pre-saved inherited anomaly and future/target CLIM to create deterministic/binary Spatial Anomaly Forecasts (SAP). Takes in sys args 1. Hemisphere and 2. (optional) Which year. Then runs forecast for all initialisation files within that year
ccc=as.integer(args[1]) #Choice of Hemisphere
if(ccc == 1)  HEM="nh"
if(ccc == 2)  HEM="sh"

Ylist=1989:2021
if(length(args)==2) {
  temp=as.integer(args[2])
  if((temp>=1989)&(temp<=2022))  Ylist=temp #Will be just a single year now
  remove(temp)}


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

# yrs=2014
for(yrs in Ylist){
  filelist=sort(Sys.glob(paste(HEMPATH,"/Outputs/savedSIP/SIPinheritedetcfor",yrs,"*",sep = "")))
  if(length(filelist)==0) {next()}
  
  for(listi in 1:length(filelist))   #listi=1
  {
    tic("Deterministic forecasting, one initialisation step time: ")  
    
    load(filelist[listi])  #saved SIP and climatology of initial date
    savename=sprintf("%s/Outputs/Forecasts/determinForecastfor%s_yod%03d",HEMPATH,ICdate,yod)
    if(file.exists(savename)) next()  #If this file already exists, we don't have to redo it. 
    
    inhSIPcorrected=SIPreturn[3,]  #neighbourhood corrected inherited Anomalies
    
    ### Preparing the necessary future climatologies + initial condition to compare against
    Climaenv=new.env()
    AllClimaArr=array(dim=c(length(SIPclima),366))
    AllClimaArr[,1]=SIPclima
    AllInitialArr=array(dim=c(length(SIPclima),366))
    AllInitialArr[,1]=ICsic
    ObsExist=array(1,dim=366)
    
    for(i in 1:365) #going through the next 365 days of leadtime (such that leadtime 0 is the initial date), for forecasting and also for future SPS
    {
      seekdate=format(as.Date(ICdate,"%Y%m%d")+i,"%Y%m%d")  
      ##Here, the ICdate + automatically skips Feb29 for not leap years. That does not seem to cause any issue. 
      
      YrClimafile=Sys.glob(sprintf("%s/Outputs/Climatology/Filtered_Climatologyfor%s_*",HEMPATH,seekdate))
      #So for each leadtime, we want the climatology, to forecast, and the true condition to compare
       if(length(YrClimafile)==0){  #ie there is no Clima file with IC
        AllInitialArr[,i+1]=NaN
        ObsExist[i+1]=0
        YrClimafile=Sys.glob(sprintf("%s/Outputs/Climatology/NOICFiltered_Climatologyfor%s_*",HEMPATH,seekdate))  #Maybe there is a climafile without IC?
        if(length(YrClimafile)==0)  {AllClimaArr[,i+1]=NaN;next()}
      }
      
      load(YrClimafile,Climaenv)  
      AllClimaArr[,i+1]=Climaenv$SIPclimaGFilt
      #If we DID find the clima file with IC, this should be 1 so we should save the ICsic
      if(ObsExist[i+1]) AllInitialArr[,i+1]=Climaenv$ICsic
    }  
    
    ICbinaryall=binarise(AllInitialArr,0.15)
    ClimaMedian=binarise(AllClimaArr,0.50)
    
    ### Actual Forecasting section
    Forecast=array(dim = c(length(SIPclima),366))
    for (leadtime in 1:366){  
      fcstsicondition=inhSIPcorrected*0
      fcstsicondition[AllClimaArr[,(leadtime)]>inhSIPcorrected]=1
      Forecast[,leadtime]=fcstsicondition
    }
    
    ### Verification/ SPS section
    SPSForecast=SPS_clim(grd$Nodeareas,Forecast,ICbinaryall)
    SPSClimatological=SPS_clim(grd$Nodeareas,AllClimaArr,ICbinaryall)
    SPSClimMedian=SPS_clim(grd$Nodeareas,ClimaMedian,ICbinaryall)
    SPSPersist=SPS_persist(grd$Nodeareas,ICbinaryall)
    
    ## If we want to plot the SPS :
    # 
    # plot(SPSForecast,type = "line",xlim=c(1,120),col="blue",)
    # lines(SPSPersist,col="red")
    # lines(SPSClimMedian,col="pink")
    # lines(SPSClimatological,lty=3,col="green")
    
    SPS_list=list(SPSClimMedian=SPSClimMedian,SPSPersist=SPSPersist,SPSClimatological=SPSClimatological,SPSForecast=SPSForecast)
    save(Forecast,AllClimaArr,AllInitialArr,SIPreturn,SIPclima,SPS_list,grd,ICdate,file = savename,version = 2)
    print(paste("Saved deterministic forecast file: ",basename(savename),sep = ""))
    remove(Forecast,AllClimaArr,AllInitialArr,SIPclima)
    toc()
  }
}