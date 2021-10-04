### Script to setup our system and create necessary grid files to generate SDAP forecasts.



# We need these two R packages to use the method.
require(ncdf4)  #To open netcdf files
require(spheRlab) #To generate contours, plots and other grid functions.
#https://github.com/FESOM/spheRlab
require(tictoc) #Used for timing


# We will use one specific path for all our data operations. Please choose this appropriately to handle storage issues. Softlinks are also ok.
MASTERPATH="~/WORK/Data/SDAP/"
if(!dir.exists(MASTERPATH)) stop(sprintf("Data directory in MASTERPATH (%s) wasn't found. Please check",MASTERPATH))

#Within this, we will divide our outputs by hemisphere.
HEM="sh"
HEMPATH=paste0(MASTERPATH,"/",HEM)
dir.create(paste0(HEMPATH,"/Outputs/"),recursive = T)

# Now, we need the input data. For this example, we will use the SIC records from OSI SAF. Other SIC records can also be used, although some troubleshooting might be required in further steps. OSI SAF data can be downloaded from here:
# ftp://osisaf.met.no/reprocessed/ice/conc/v2p0 #OSI 450  (for years 1979 to 2015)
# ftp://osisaf.met.no/reprocessed/ice/conc-cont-reproc/v2p0/ #OSI 430b (for years 2016 on)

# For this example, we will assume that the data has already been downloaded in these locations:
Datapath1=paste0(MASTERPATH,"/OSI450/osisaf.met.no") #OSI 450 
Datapath2=paste0(MASTERPATH,"/OSI430b/osisaf.met.no")#OSI 430b
# OSI SAF 450 was continued as OSI 430b after 2016, therefore we have 2 datapaths here. They can also be merged into one directory, but here we have left them in the original structure.

if(!dir.exists(Datapath1)) stop("Data folder (OSI450/osisaf.met.no) NOT found! Check that there is a folder or a softlink.")
if(!dir.exists(Datapath2)) warning("Data folder 2 (OSI430b/osisaf.met.no) NOT found! Check that there is a folder or a softlink.")


#Now if things worked, we can start looking at the files:

#Let's look for a data/netcdf file
file1=Sys.glob(sprintf("%s/reprocessed/ice/conc/v2p0/????/??/*_%s*.nc",Datapath1,HEM))  
# Here, we depict whether it is a SH or NH file by using the HEM.
if(length(file1)==0) stop("No data file (netcdf) found. Either there is an issue with the path or something in the script is not correct. Check!")
print(paste0("Things look mostly good so far! Making the grid file now."))


##### Section to create the grid/land-lake mask, and save gridfile that can be reused repeatedly-----------------------

gridFilename=paste(HEMPATH,"/Outputs/gridfile_OSISAF450",HEM,"_all_inclFlags",sep = "")


#First, find ALL land or lake points, by combing through all the files. Data otherthan OSI SAF might have other flags to create mask, so this part should be re-checked.  
allland=NULL
alllake=NULL

N=length(file1)
if(N>1000) N=500 #If we have more than 1000 files, just check the first 500

for(i in seq(1,N,5)){  #Check every 5th file
  fl=nc_open(file1[i])  
  flag = ncvar_get(fl,"status_flag")
  nc_close(fl)
  
  flag = as.vector(t(flag))
  land = which(flag %in% c(1))
  lake = which(flag %in% c(2,6))  #2 is lake flag, 6 is lake + open water flag, so use both
  
  # If these points are not already on our list, get them.
  alllake=c(alllake,lake[!lake %in% alllake])
  allland=c(allland,land[!land %in% allland])
  
}

### Now, let's take one specific file and save the grid from that:

fl = nc_open(file1[1])
lat = ncvar_get(fl,"lat")
lon = ncvar_get(fl,"lon")
sic = ncvar_get(fl,"ice_conc")/100
flag = ncvar_get(fl,"status_flag")
nc_close(fl)

Nrow = dim(lat)[1]
Ncol = dim(lat)[2]
# sic = as.vector(t(sic))
flag = as.vector(t(flag))
###
totnodes=length(flag)

grd.full = sl.grid.curvilin2unstr(lon = lon, lat = lat, close.sides = FALSE)
grd = grd.full

# sample plot
# pir = sl.plot.init(projection = "polar", polar.latbound = 65, do.init.device = FALSE, col.background = "#F0F0F0")
# sl.plot.naturalearth(pir, what="coastline", resolution="medium")
# sl.plot.end(pir, do.close.device = FALSE)

##land/lake from just one timestep was this:
# land = which(flag %in% c(1))
# lake = which(flag %in% c(2,6))  #2 is lake flag, 6 is lake + open water flag, so use both

# sl.plot.points(pir, lon = grd$lon[allland], lat = grd$lat[allland], col="brown")
# sl.plot.points(pir, lon = grd$lon[alllake], lat = grd$lat[alllake], col="blue")

landorlake = unique(c(alllake,allland))  #This is the superset of all land or lake points from all over the year that was found earlier.

grd = sl.grid.reduce(grd, remove.points = landorlake, set.coast = TRUE)
nodes.kept = grd$reduce.kept$nodes
grd.postlandlake.prebaynoderemoval = grd

# plot the data after removal of land (including lakes)
# cb2 = sl.plot.field.elem(pir, num = sic[nodes.kept], lon = grd$lon, lat = grd$lat, elem = grd$elem,
#                          colbar = sl.colbar(cols = c("darkred","red")))


# the following is somewhat ugly and only required because sl.contours() can not handle isolated elements and
# nodes that connect more than two domains (for example two ocean and two land domains)
grd = grd.postlandlake.prebaynoderemoval
remove.pnts = "dummy"
while (length(remove.pnts) > 0) {
  remove.pnts = which(rowSums(matrix(!grd$coast[grd$neighnodes],ncol=ncol(grd$neighnodes)),na.rm=TRUE) == 0)
  # sl.plot.points(pir, lon = grd$lon[remove.pnts], lat = grd$lat[remove.pnts], col="purple")
  if (length(remove.pnts) > 0) {
    print(paste("removing",length(remove.pnts),"bay nodes"))
    grd = sl.grid.reduce(grd, remove.points = remove.pnts, set.coast = TRUE)
    nodes.kept = nodes.kept[grd$reduce.kept$nodes]
  } else {
    checkres = sl.findneighbours(grd$elem)
    if (!is.null(checkres$elems.completed)) {
      remove.pnts = unique(grd$elem[!checkres$elems.completed])
      print(paste("removing",length(remove.pnts),"multi-domain nodes"))
      grd = sl.grid.reduce(grd, remove.points = embayment.nodes, set.coast = TRUE)
      nodes.kept = nodes.kept[grd$reduce.kept$nodes]
    }
  }
}

#Now, we can save this gridfile (incl the original maybe) and re-use it for all other time-steps HOPING that the grids are about the same.
Gridfromwhichfile=file1[1]
# But note that we are using a bunch of files to look at the flags/nodes.kept 

###Calculating areas for the grid:###
Elementareas=array(dim=dim(grd$elem)[[1]])
for (i in 1:length(Elementareas)){
  nodes=grd$elem[i,]   #What nodes does this element have?
  Elementareas[i]=sl.triag.area(grd$lon[nodes],grd$lat[nodes]) #And how much is the weight here?
}

Nodeareas=array(dim=length(grd$lon))
for (i in 1:length(Nodeareas)){
  elems=grd$neighelems[i,]  #Most nodes have upto 6 elements around it, 
  areaofelems=Elementareas[elems]
  areaofelems=areaofelems[!is.na(areaofelems)]  #but some have less than 6, so some of those points could be NA
  Nodeareas[i]=sum(areaofelems)*(1/3)  #Each element area is shared by 3 nodes, and nodes get shares from all the elements its part of.
}
grd$Nodeareas=Nodeareas  #I checked that total Nodearea = Total Elementarea, so this should be correct
grd$Elementareas=Elementareas

save(file = gridFilename,Gridfromwhichfile,grd.full,alllake,allland,grd.postlandlake.prebaynoderemoval,grd,nodes.kept,version = 2)

print(paste0("Grid file saved:",basename(gridFilename)))

######### Make other folders that we need: --------------

dir.create(paste0(HEMPATH,"/Outputs/Climatology"))
dir.create(paste0(HEMPATH,"/Outputs/savedSIP"))
dir.create(paste0(HEMPATH,"/Outputs/Forecasts"))
dir.create(paste0(HEMPATH,"/Outputs/Alpha"))
dir.create(paste0(MASTERPATH,"/Figs"))

## These folders are needed if we want to compare whether NO gaussian filtering is better:
# dir.create(paste(HEMPATH,"/Outputs/NOGausssavedSIP",sep = ""))
# dir.create(paste(HEMPATH,"/Outputs/NOGaussForecast",sep = ""))
# dir.create(paste(HEMPATH,"/Outputs/NOGaussAlpha",sep = ""))


# We did a gaussian averaging of our climatology maps, for which a gaussian weight file must be generated once.
######### Make gaussian weight file: --------------
print(paste("Save the grid file! Making the gaussian weight file now.",sep = ""))

tic("Gaussian weight making time:")
Rsphere=1  #Native grid units ()
gauss.sigma=2*pi*Rsphere/360
cutoff=2 #*gauss.sigma # One can change this later!
gsnwtfile=sprintf("%s/Outputs/gaussianweights%s_R%i_cutoff%i",HEMPATH,HEM,Rsphere,cutoff)

sf.gw.res=sl.spatialfilter.getweights(lon=grd$lon,lat=grd$lat,neighmat = grd$neighnodes,areas = grd$Nodeareas,Rsphere = 1,gauss.sigma = gauss.sigma,cutoff = cutoff*gauss.sigma)#for this specific grid, the weights are the same. it results in this list: sf.gw.res
save(file = gsnwtfile,sf.gw.res,version = 2)  
print(paste0("Gaussian weight saved: ",basename(gsnwtfile)))
toc()
print("All done with Setup!")