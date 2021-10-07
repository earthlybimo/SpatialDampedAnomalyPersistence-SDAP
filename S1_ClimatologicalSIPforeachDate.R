#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
### This script saves the Climatological SIP for each date using the prev 10 years of SIC. It is designed to be run as part of a batch process (slurm/cronjob), with command line arguments giving 1. Hemisphere, 2. DOY(Date) and 3. Year. If only 2 args are given, all years from 1989 to 2021 are used. 
#The script also uses the previously saved gaussian weights and grid information (saved in S0_Setup.R). The script can take some time to run, esp for multiple years. It is easy to use this with sbatch to run multiple instances for each DOY.


## First, to figure out command line arguments!
if(length(args)==0)
  stop ("Did not get an argument in !? Try re-running with 1. Hemisphere, 2. DOY (1:366) and 3. (optional) year")

ccc=as.integer(args[1]) #Choice of Hemisphere

yod=as.integer(args[2]) #DOY
Ylist=1989:2021
if(length(args)==3) {
  temp=as.integer(args[3])
  if((temp>=1989)&(temp<=2022))  Ylist=temp #Will be just a single year now
  remove(temp)}

if(ccc == 1) HEM="nh"
if(ccc== 2)  HEM="sh"

#Save paths:
MASTERPATH="~/WORK/Data/SDAP/"
HEMPATH=paste0(MASTERPATH,"/",HEM)
Datapath1=paste0(MASTERPATH,"/OSI450/osisaf.met.no") #OSI 450 
Datapath2=paste0(MASTERPATH,"/OSI430b/osisaf.met.no")#OSI 430b


require(ncdf4)
require(tictoc)
require(spheRlab)

binarise <-function (somearr,dlevel) {  #Function to binarise some given array  based on  this level
  ll=dim(somearr)
  if (is.null(ll)){ll=length(somearr)}  #because 1d arrays dont have dim, they have lengths
  ibinar=array(dim =ll)
  ibinar[somearr[]>=dlevel]=1
  ibinar[somearr[]<dlevel]=0
  return (ibinar)
}

#Gaussian weight file
gsnwtfile=sprintf("%s/Outputs/gaussianweights%s_R*",HEMPATH,HEM)  
#If there are multiple versions of the weight file, it might be worth specifying which one

#Load grid
gridFilename=paste0(HEMPATH,"/Outputs/gridfile_OSISAF450",HEM,"_all_inclFlags")
load(gridFilename,envir = (gridenv=new.env()))
nodes.kept=gridenv$nodes.kept

#### Make climatology ############
for (yearIC in Ylist){  #Which years? 
  tic("Timing for saving one climatology  ")
  dayy=format(as.Date(yod,origin="2015-12-31"),"%d") 
  monthh=format(as.Date(yod,origin="2015-12-31"),"%m")
  ICdate=sprintf("%d%s%s",yearIC,monthh,dayy)  #as.Date(ICdate,"%Y%m%d")
  
  savename=sprintf("%s/Outputs/Climatology/Filtered_Climatologyfor%s_climatology%dto%d_yod%03d",HEMPATH,ICdate,(yearIC-10),(yearIC-1),yod)
  if(!length(Sys.glob(savename))==0) next()  #Savefile already exists.
  
  #First, does the observed data for this date exist?
  ICsic=NA;ICexists=0;sic_original=NA #Default, in case the file doesn't exist for this date.
  ICfile1=Sys.glob(sprintf("%s/reprocessed/ice/conc-cont-reproc/v2p0/%d/%s/*_%s_*%d%s%s1200.nc",Datapath2,yearIC,monthh,HEM,yearIC,monthh,dayy))  #For files from newer bunch aka OSI430b
  if(length(ICfile1)==0) ICfile1=Sys.glob(sprintf("%s/reprocessed/ice/conc/v2p0/%d/%s/*_%s_*%d%s%s1200.nc",Datapath1,yearIC,monthh,HEM,yearIC,monthh,dayy))  #For files from before 2016
  
  if(!length(ICfile1)==0)  ## IF the file does exist (doing this weird way because file.exists doesn't like wildchars)
  {##open file:
    fl = nc_open(ICfile1)
    ICsic = ncvar_get(fl,"ice_conc")/100
    flag = ncvar_get(fl,"status_flag")
    ICsic = as.vector(t(ICsic))
    flag = as.vector(t(flag)) #Could be used to REcheck that we don't have a land/lake flagged point. Not sure if it is worth the extra effort.
    nc_close(fl)  
    
    ICexists=1
    sic_original=ICsic
    ICsic=ICsic[nodes.kept]
    if(any(is.na(ICsic))) warning(paste("Mismatch with the flags in this file: ",basename(ICfile1)," . This will most likely cause issues later so check!",sep = "")) # Is there any NA in the data now? After the land/lakes have been removed. So far found none!
  }
  
  if((!length(Sys.glob(sprintf("%s/Outputs/Climatology/NOICFiltered_Climatologyfor%s_climatology%dto%d_yod%03d",HEMPATH,ICdate,(yearIC-10),(yearIC-1),yod)))==0) & (ICexists==0)) next()  #Old saved NOIC file exists, and IC file still doesn't exist, so skip. Else continue!
  
  
  # Assuming all is well upto here!
  
  #Now we need to build a climatology for this date
  sicarr=array(dim=c(10,length(gridenv$grd$lat)))
  
  yrcnt=1
  for (year in (yearIC-10):(yearIC-1)){
    file1=Sys.glob(sprintf("%s/reprocessed/ice/conc/v2p0/%d/%s/*_%s_*%d%s%s1200.nc",Datapath1,year,monthh,HEM,year,monthh,dayy))
    if(length(file1)==0) file1=Sys.glob(sprintf("%s/reprocessed/ice/conc-cont-reproc/v2p0/%d/%s/*_%s_*%d%s%s1200.nc",Datapath2,year,monthh,HEM,year,monthh,dayy))
    
    if(length(file1)==0) {  #If both steps failed, the file maybe doesn't exist
      message(sprintf("File missing for : %04d-%s-%s",year,monthh,dayy))
      next()}
    tryCatch(  #Because some files give weird isssues, we are locking them in this error catching method.
      {fl = nc_open(file1)
      sic_cl = ncvar_get(fl,"ice_conc")/100
      # flag = ncvar_get(fl,"status_flag")
      sic_cl = as.vector(t(sic_cl))
      # flag = as.vector(t(flag)) 
      #Again, could be used to re-check but more extra effort.
      nc_close(fl)},
      error=function(error_message){
        message("\n Apparently there was an error, file missing or corrupted maybe. Will skip this.")
        message(sprintf("File should be for : %04d-%s-%s",year,monthh,dayy))
        message(error_message)
      })
    
    # sic_cl_original=sic_cl  #But no use keeping it really
    sic_cl=sic_cl[nodes.kept]
    if(any(is.na(sic_cl))) {
      print( "Mismatch! In the climatology, with this file:")
      print(basename(file1)) # Is there any NA in the data now? After the land/lakes have been removed. So far found none!
    }
    sicarr[yrcnt,]=binarise(sic_cl,0.15)
    yrcnt=yrcnt+1
  }
  
  SIPclima0=colMeans(sicarr,na.rm = TRUE)
  
  ## Gaussian filtering !
  
  gfile=Sys.glob(gsnwtfile)
  if(length(gfile)==0){
    warning(paste("Gaussian weight file:",gsnwtfile," is missing or couldn't be reached! Gaussian filtering is skipped. Please check.",sep = ""))
    SIPclimaGFilt=NA
    
  }else{
    load(gfile)
    print("Starting Gaussian Filtering")
    SIPclimaGFilt=sl.spatialfilter(SIPclima0, sf.gw.res)
  }
  
  ICdate=sprintf("%d%s%s",yearIC,monthh,dayy)  #as.Date(ICdate,"%Y%m%d")
  savename=sprintf("%s/Outputs/Climatology/Filtered_Climatologyfor%s_climatology%dto%d_yod%03d",HEMPATH,ICdate,(yearIC-10),(yearIC-1),yod)
  if(ICexists==0) savename=sprintf("%s/Outputs/Climatology/NOICFiltered_Climatologyfor%s_climatology%dto%d_yod%03d",HEMPATH,ICdate,(yearIC-10),(yearIC-1),yod)
  
  save(file=savename,ICexists,SIPclimaGFilt,ICsic,yod,sicarr,SIPclima0,sic_original,ICdate,version = 2)
  print(paste("Saved Climatology (and initial) for date:",ICdate," as file ",basename(savename),sep = ""))
  remove(SIPclima0,sic_original,sicarr,SIPclimaGFilt,ICsic);invisible(gc())
  toc()
} #End of the year loop
