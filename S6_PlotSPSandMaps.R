#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Script to take pre-saved SDAP forecast files and plot the SPS timeseries for the following 120 days, alongside a map of the forecasted sea-ice condition for a specific lead time if given.

#The command line system arguments are 1. Hemisphere, 2. Initialisation month and 3. Initialisation year. Optionally, a 4th integer - leadtime- can be given to plot a map of the SDAP forecast for the given combination of initial date and lead-time. 

#For a simple test, skip to line 39 and directly load DampedForecast file,  and Give HEM ("nh" or "sh"). For plotting the forecast map on a given date, give lt (leadtime) integer for plotting a map and Figpath for saving appropriately.

ccc=as.integer(args[1]) #Choice of Hemisphere
if(ccc == 1){
  HEM="nh"
  Alphasetname="NewTraining_2021NH_yrs89to98"}

if(ccc== 2){
  HEM="sh"
  Alphasetname="NewTraining_2021SH_yrs89to98"}

mc=as.integer(args[2]) #Choice of initialisation date (each doy corresponds to the start of a particular month)
yodlist=c(1,32,61,92,122,153,183,214,245,275,306,336)
yodi= yodlist[mc]

yearIC=as.integer(args[3]) #Initial year choice
if(length(args)==4) lt=as.integer(args[4]) #Choice of leadtime
# if(get0("lt"))

MASTERPATH="~/WORK/Data/SDAP/"
HEMPATH=paste0(MASTERPATH,"/",HEM)
Figpath=paste0(MASTERPATH,"/Figs") #Figures will be saved here
if(!dir.exists(Figpath)) dir.create(Figpath)

dayy=format(as.Date(yodi,origin="2015-12-31"),"%d") 
monthh=format(as.Date(yodi,origin="2015-12-31"),"%m")
ICdate=sprintf("%d%s%s",yearIC,monthh,dayy)

loadname=sprintf("%s/Outputs/Forecasts/DampedForecast_%s_%s_yod%03d",HEMPATH,ICdate,Alphasetname,yodi)
if(!file.exists(loadname)) stop("SDAP forecast file for given arguments not found. Please re-check.")

load(loadname)

if(is.null(get0("ICdate"))) ICdate=detFcenv$ICdate

# Step 1, let's plot the SPS for the first 120 days:
Fctr=((6371^2)/1000000)  #For the pre-saved data, this scale factor translates SPS values to Million square kilometers. For other data, please replace with the appropriate scale factor.

ymx=1.5 
if(HEM=="sh") ymx=3

par(mar=c(5,5,4,1)+.1)
plot(Fctr*detFcenv$SPS_list$SPSClimatological,col='grey',type='line',xlab = "Lead Time (Days)",ylab=expression(paste("SPS (10"^"6"," km"^"2",")")),ylim = c(0,ymx),xlim=c(1,120),main=paste0("SPS for forecasts initialised on ",ICdate),cex.axis=1.5,cex.lab=1.5, yaxt="n",lwd=2,yaxs="i")
lines(1:366,(Fctr*detFcenv$SPS_list$SPSForecast),col='blue',lwd=3)
lines(1:366,(Fctr*detFcenv$SPS_list$SPSPersist),col='red',lwd=2)
lines(1:366,Fctr*SPSForecastwdamping,col='black',lwd=3)
# lines(1:366,(Fctr*detFcenv$SPS_list$SPSClimMedian),col='grey',lwd=3,lty=3)
legend("bottomright",col = c('red','grey','blue','black'),legend = c('PERS','CLIM','SAP (Deterministic) Fcst','SADP (Probabilistic) Fcst'),lty=c(1),lwd=2,bty = "n")
abline(v=c(1,seq(30,360,30)),lty=3,col="grey")
if(HEM=="nh") {axis(2,at=seq(0,by = 0.25,to= 4),cex.axis=1.5);abline(h = seq(0,by=0.125,to= 4),lty=3,col="grey")}
if(HEM=="sh") {axis(2,at=seq(0,by = 0.5,to= 4),cex.axis=1.5);abline(h = seq(0,by=0.25,to= 4),lty=3,col="grey")}

## Note that the CLIM here changes a lot depending on the initialisation date. This is different than the CLIM as shown in the manuscript plot, where a constant line represents climatological SPS averaged over all dates.  
# dev.off() or save appropriately.


if(get0("lt")) {  #If lt (leadtime) exists
  
  library(RColorBrewer);library(spheRlab)
  grd=detFcenv$grd
  
  coli=rev(brewer.pal(9,name = "PuBu")) 
  testcol=sl.colbar(coli[3:9],10)  #Colorbar we are using, adjusted
  brks=seq(0.1,0.9,length.out = 9) 
  sl.plot.colbar(testcol,breaks = brks,file.name = paste0(Figpath,"SeaIceProbabilityColorbar.pdf"))
  #To save the colorbar
  
  #Actual map:
  if(HEM=="nh") pir = sl.plot.init(projection="polar",polar.latbound = 50,file.name =sprintf("%sMap_DampedFcst_%s_plusleadtime%03d_NH.pdf",Figpath,detFcenv$ICdate,lt),col.background = "grey")
  if(HEM=="sh") pir = sl.plot.init(projection="polar",polar.latbound = 50,polar.lonlatrot =c(0,-90,0),file.name =sprintf("%sMap_DampedFcst_%s_plusleadtime%03d_SH.pdf",Figpath,detFcenv$ICdate,lt),col.background = "grey")
  #Change polar.latbound to make coverage larger or smaller.
  
  pcol = sl.plot.field.elem(pir,num = (Forecast_wDamping[,(1+lt)]),lon = grd$lon, lat=grd$lat, elem=grd$elem,colbar = testcol,colbar.breaks = brks,border.lwd = 1,border = T)
  sl.plot.naturalearth(pir, what="land", resolution="coarse",fill.col = "grey")
  sl.plot.naturalearth(pir, what="coastline",resolution="medium",lines.col = "black",lwd = 2)
  res = sl.plot.lonlatgrid(plot.init.res = pir,labels = T, pole.hole=TRUE,lat.0 = 50,lat.distance = 15,lon.0 = 0,lon.distance = 30,labels.lat.every = 2,labels.lon.every = 4,lty=3,lwd = 0.5,col = "gray30")
  res = sl.plot.contours(plot.init.res = pir,contours.res = sl.contours(detFcenv$AllInitialArr[,(1+lt)],lon = grd$lon, lat=grd$lat, elem=grd$elem,levels = 0.15),col = "red",lwd = 2)  #TRUTH / observed actual contour
  SPS=SPSForecastwdamping[(lt+1)]*Fctr
  sl.plot.text(pir,lon =150 ,lat=65,labels = sprintf("SPS = %0.3g mill sq km",SPS))
  
  #If we want the CLIM-MED contour on there:
  newClimCntrs=sl.contours(detFcenv$AllClimaArr[,(1+lt)],lon = grd$lon, lat=grd$lat, elem=grd$elem,levels = c(0.1,0.5,0.9))
  res = sl.plot.contours(plot.init.res = pir,contours.res = newClimCntrs,indices = 2,col = "black",lty = 3,lwd=2)
  
  sl.plot.end(pir)
  
  print(sprintf("Map of SDAP (probabilistic forecast) saved here: %sMap_DampedFcst_%s_plusleadtime%03d.pdf",Figpath,detFcenv$ICdate,lt)) 
}

