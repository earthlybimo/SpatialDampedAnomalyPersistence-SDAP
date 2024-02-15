#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
### This script uses the Climatological SIP and Initial condition for each date to find the SIP anomalies. It is designed to be run as part of a batch process, with system args giving 1. Hemisphere, 2. Month no. and 3. Year. If only 2 syst args are given, all years from 1989 to 2021 are used. Currently it is established to only initialise at the start of each month, but this can be easily changed. 

## First, to detect arguments !
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
# Datapath1=paste0(MASTERPATH,"/OSI450/osisaf.met.no") #OSI 450   Not needed anymore
# Datapath2=paste0(MASTERPATH,"/OSI430b/osisaf.met.no")#OSI 430b


require(ncdf4)
require(spheRlab)
require(tictoc)
require(plyr)
#library(parallel)  #Certain functions can also be parallelised to speed process, but optional

#Declaring some functions beforehand for ease of usage
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

anomalyatIE <-function (cnt,sip_clim) {  #Function to find anomalies (from sip_clim) at initial IceEdge contour (cnt)
  for(i in 1:length(cnt[[1]]$segments)){ #Now i is each segment of the edge
    ICprob=vector()
    for(j in 1:length(cnt[[1]]$segments[[i]]$lat)){  #For each point in this segment
      #This has to do with how contours are depicted in spheRlab. Essentially, they are somewhere between actual grid points, and the ratio of distance is also given. 
      edgepts=cnt[[1]]$segments[[i]]$edge.endpoints[j,]
      edgesic=sip_clim[edgepts]  
      wgtoffirst=cnt[[1]]$segments[[i]]$edge.relpos[j]
      sipforpt=(edgesic[[1]]*(1-wgtoffirst))+(edgesic[[2]]*(wgtoffirst)) #For properly weighing the climatologies at both points
      ICprob[j]=sipforpt
      remove(edgepts,edgesic,wgtoffirst,sipforpt)
    }
    cnt[[1]]$segments[[i]]$ICprobabofnearpt=ICprob
    remove(ICprob)
  }
  return (cnt)
}

#Load grid
gridFilename=paste0(HEMPATH,"/Outputs/gridfile_OSISAF450",HEM,"_all_inclFlags")
load(gridFilename,envir = (gridenv=new.env()));grd=gridenv$grd
lon=grd$lon;lat=grd$lat;grlen=length(lon)
remove(gridFilename)
# yrs=2001; # yodi=1  #For testing
for(yrs in Ylist){
  for(yodi in DOYlist){# Or for initialising at all days, 1:366
    Climatologyfile=Sys.glob(sprintf("%s/Outputs/Climatology/Filtered_Climatologyfor%d*_yod%03d",HEMPATH,yrs,yodi))
    if(length(Climatologyfile)==0) {warning(paste("Climafile not found for year ",yrs," day ",yodi,sep = ""));next()}
    load(Climatologyfile)
    
    savename=sprintf("%s/Outputs/savedSIP/SIP_anomaly_initialisedfor_%s",HEMPATH,ICdate)
    if(file.exists(savename)) next()  #If this file already exists, don't have to redo it.
    
    print(paste0("Saving anomalies for Day ",yodi," year ",yrs))
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
    ) #End of tryCatch  
    if(err==TRUE) next() 
  
    
    ## Here, we are using a shortcut and simply taking the CLIM_sip and comparing it later in the forecasting phase, rather than taking the anomaly against median (aka CLIM - 50%) and compare against 50% later. In effect, the result is the same. 
    
    #First inherit climaSIP to Iniital contour
    ICcnt=anomalyatIE(ICcnt0,SIPclima)  #climatological SIP of nearest grid along each contour pnt is added
    ICcnt_df <- ldply(ICcnt[[1]]$segments,data.frame)  #Make this into a dataframe (necessary later)
    
    ## Second pass to median contour
    cliMedncnt_df <- ldply(climMedcnt[[1]]$segments,data.frame)  
    #Make the CLIMMED contour into a dataframe, for easier dist measurements 
    
    climMed_inheritedSIP=array(dim=length(cliMedncnt_df$lat))
    for(i in 1:length(cliMedncnt_df$lat)){ #For every point of the CLIM-MED cntur
      nearpt=findclosestpointindf(cliMedncnt_df$lat[[i]],cliMedncnt_df$lon[[i]],ICcnt_df) 
      climMed_inheritedSIP[i]=ICcnt_df$ICprobabofnearpt[[nearpt]]}

    ## Then pass to all points on the grid: (note that we are only using ocean points)
    inheritedSIP=array(dim=grlen)
    for(i in 1:grlen){ #For each gridpoint
      nearpt=findclosestpointindf(lat[[i]],lon[[i]],cliMedncnt_df)
      inheritedSIP[i]=climMed_inheritedSIP[nearpt]
    }
    ## Last, do some local corrections based on initial forecast (Initial state correction)
    sicICbin=binarise(ICsic,0.15) #Only ice or no ice in the initial date
    correctedSIP=inheritedSIP
    eps=0.05  #eps=Tolerance for correction down below.
    for(i in 1:grlen){ #For each grid point
      #If initial ice, and SIPclima is NOT >= SIPlimt, then that means NOT ice fcst, so need to correct
      if((sicICbin[i]==1)&(SIPclima[i]<correctedSIP[i])) correctedSIP[i]=SIPclima[i]-eps
      if((sicICbin[i]==1)&(SIPclima[i]<correctedSIP[i])+eps) correctedSIP[i]=correctedSIP[i]-eps
      #If NO ice, and SIPclima is>= SIPlimt, then that means ice fcst, so need to correct
      if((sicICbin[i]==0)&(SIPclima[i]>=correctedSIP[i])) correctedSIP[i]=SIPclima[i]+eps
      if((sicICbin[i]==0)&(SIPclima[i]>=correctedSIP[i]-eps)) correctedSIP[i]=correctedSIP[i]+eps
    }
    
    
    #For checking, forecast for day1 is:
    # fcst=SIPinheritedI*0;fcst[SIPclima>SIPcorrected]=1;
    # Can check mismatch with which(fcst>sicICbin)  and which(fcst<sicICbin)
    
    ## Saving!
    SIPreturn=array(dim=c(3,grlen)) #In some old scripts, this format was used so this is left here
    SIPreturn[1,]=inheritedSIP
    SIPreturn[3,]=correctedSIP 
    
    Comments=paste0("First taking anomalies from IC to MP by simple nearest neighbour from MP. Then correcting points with false initial forecast, and adding a epsilon buffer to inherited corrected SIP if it does not have enough. Using epsilon: ",eps)
    
    save(SIPreturn,yod,Climatologyfile,ICsic,SIPclima,grd,ICdate,inheritedSIP,correctedSIP,Comments,file = savename,version = 2)
    print(paste("Saved file: ",basename(savename),sep = ""))
    remove(SIPreturn,SIPclima,climMedcnt,ICcnt,ICcnt0,ICcnt_df,cliMedncnt_df,ICsic)
    toc() 
  } #Loop over other initialisation dates and years
}