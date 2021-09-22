#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
### This script uses the Climatological SIP and Initial condition for each date to find the SIP anomalies. It is designed to be run as part of a batch process, with system args giving 1. Hemisphere, 2. Month no. and 3. Year. If only 2 syst args are given, all years from 1989 to 2021 are used. Currently it is established to only initialise at the start of each month, but this can be easily changed. 

## First, to detect arguments and all!
if(length(args)==0)
  stop ("Did not get an argument in !? Try re-running with 1. Hemisphere, 2. (optional) Monthno. (1:12) and 3. (optional) year")

ccc=as.integer(args[1]) #Choice of Hemisphere

DOYlist=c(1,32,61,92,122,153,183,214,245,275,306,336)
if(length(args)>1) {
  mch=as.integer(args[2])
  temp=DOYlist #default is all months
  if((mch>=1)&(mch<=12))  temp=DOYlist[mch] #Will be just a single month now
  DOYlist=temp
  remove(temp,mch)}

Ylist=1989:2021
if(length(args)==3) {
  temp=as.integer(args[3])
  if((temp>=1989)&(temp<=2022))  Ylist=temp #Will be just a single year now
  remove(temp)}

if(ccc == 1) HEM="nh"
if(ccc== 2)  HEM="sh"

#Same as in every script:
MASTERPATH="~/WORK/Data/SDAP/"
HEMPATH=paste0(MASTERPATH,"/",HEM)
Datapath1=paste0(MASTERPATH,"/OSI450/osisaf.met.no") #OSI 450 
Datapath2=paste0(MASTERPATH,"/OSI430b/osisaf.met.no")#OSI 430b


require(ncdf4)
require(tictoc)
require(spheRlab)
library(parallel)  #Certain functions can also be parallelised to speed process

binarise <-function (somearr,dlevel) {  #Function to binarise some given array  based on  this level
  ll=dim(somearr)
  if (is.null(ll)){ll=length(somearr)}  #because 1d arrays dont do 
  ibinar=array(dim =ll)
  ibinar[somearr[]>=dlevel]=1
  ibinar[somearr[]<dlevel]=0
  return (ibinar)
}
findclosestpointindf <-function (plat,plon,df) {  #Function to find nearest point within a dataframe for any given point
  gcdistances=sl.gc.dist(c(plon,df$lon),c(plat,df$lat),sequential = FALSE)  
  return (which.min(gcdistances))
}

#Load grid
gridFilename=paste0(HEMPATH,"/Outputs/gridfile_OSISAF450",HEM,"_all_inclFlags")
load(gridFilename,envir = (gridenv=new.env()))
nodes.kept=gridenv$nodes.kept
lon=grd$lon;lat=grd$lat;grlen=length(lon)
remove(gridFilename)


# yrs=2001
for(yrs in Ylist){
  # yodi=1  #For testing
  for(yodi in DOYlist){#1:366
    Climatologyfile=Sys.glob(sprintf("%s/Outputs/Climatology/Filtered_Climatologyfor%d*_yod%03d",HEMPATH,yrs,yodi))
    if(length(Climatologyfile)==0) {warning(paste("Climafile not found for year ",yrs," day ",yodi,sep = ""));next()}
    load(Climatologyfile)
    
    savename=sprintf("%s/Outputs/savedSIP/SIPinheritedetcfor%s",HEMPATH,ICdate)
    if(file.exists(savename)) next()  #If this file already exists, we don't have to redo it.
    
    tic("One initial step")
    SIPclima=SIPclimaGFilt  #Using the gaussian filtered version. ``````````````````
    err=FALSE
    tryCatch({  #Some files give errors with the contour making so wrapping them here
      climMedcnt=sl.contours(var = SIPclima, elem=grd$elem,lon = grd$lon,lat=grd$lat,levels = c(0.5),return.edge.info = TRUE)  #Make median contour

      ICcnt0=sl.contours(var = ICsic, elem=grd$elem,lon = grd$lon,lat=grd$lat,levels = c(0.15),return.edge.info = TRUE)  #0.15 cntr from initcondtn
      
    },
    error=function(error_message){
      message("Apparently there was an error, in one of the contour making steps:::")
      message(error_message)
      message(paste("---For date: ",ICdate,sep=""))
      #next()  This next doesnt work because of the way trycatch works
      err<<-TRUE  #So we have to do this 'hack'. Also we have to use '<<-' and not '=' for this to work!
    }
    )
    if(err==TRUE) next() 
    
    cliMedncnt_df <- ldply(climMedcnt[[1]]$segments,data.frame)  #Make this into a dataframe (necessary for later)
    
    ### Initilisation steps:
    
    #First inherit climaSIP to Iniital contour
    ICcnt=anomalyatIE(ICcnt0,SIPclima)  #climatological SIP of nearest grid along each contour pnt is added
    ICcnt_df <- ldply(ICcnt[[1]]$segments,data.frame)  #Make this into a dataframe (necessary later)
    
    ## Second inherit to median contour
    climMed_inheritedSIP=array(dim=length(cliMedncnt_df$lat))
    for(i in 1:length(cliMedncnt_df$lat)){
      nearpt=findclosestpointindf(cliMedncnt_df$lat[[i]],cliMedncnt_df$lon[[i]],ICcnt_df) 
      climMed_inheritedSIP[i]=ICcnt_df$ICprobabofnearpt[[nearpt]]
    }

    ## Then inherit to all points on the grid:
    inheritedSIP=array(dim=grlen)
    for(i in 1:grlen){
      nearpt=findclosestpointindf(lat[[i]],lon[[i]],cliMedncnt_df)
      inheritedSIP[i]=climMed_inheritedSIP[nearpt]
    }
    ## Last, do some local corrections based on initial forecast (Initial state correction)
    sicICbin=binarise(ICsic,0.15)
    correctedSIP=inheritedSIP
    eps=0.05  #eps=Tolerance for correction down below.
    for(i in 1:grlen){
      #If ice, and SIPclima is NOT >= SIPlimt, then that means NOT ice fcst, so correct
      if((sicICbin[i]==1)&(SIPclima[i]<correctedSIP[i])) correctedSIP[i]=SIPclima[i]-eps
      if((sicICbin[i]==1)&(SIPclima[i]<correctedSIP[i])+eps) correctedSIP[i]=correctedSIP[i]-eps
      #If NO ice, and SIPclima is>= SIPlimt, then that means ice fcst, so correct
      if((sicICbin[i]==0)&(SIPclima[i]>=correctedSIP[i])) correctedSIP[i]=SIPclima[i]+eps
      if((sicICbin[i]==0)&(SIPclima[i]>=correctedSIP[i]-eps)) correctedSIP[i]=correctedSIP[i]+eps
    }
    
    
    #For checking, forecast for day1 is:
    # fcst=SIPinheritedI*0;fcst[SIPclima>SIPcorrected]=1;
    # Now check mismatch with which(fcst>sicICbin)  and which(fcst<sicICbin)
    
    ## Saving!
    SIPreturn=array(dim=c(3,grlen)) #In some old scripts, this format was used so this is left here
    SIPreturn[1,]=inheritedSIP
    SIPreturn[3,]=correctedSIP 
    
    Comments=paste0("First taking anomalies from IC to MP by simple nearest neighbour from MP. Then correcting points with false initial forecast, and adding a epsilon buffer to inherited corrected SIP if it does not have enough. Using epsilon: ",eps)
    
    save(SIPreturn,yod,Climatologyfile,ICsic,SIPclima,grd,ICdate,inheritedSIP,correctedSIP,Comments,file = savename,version = 2)
    print(paste("Saved file: ",basename(savename),sep = ""))
    remove(SIPreturn,SIPclima,climMedcnt,ICcnt,ICcnt0,ICcnt_df,cliMedncnt_df,ICsic)
    toc() 
  }
}