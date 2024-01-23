args = commandArgs(trailingOnly=TRUE)
ccc=as.integer(args[1]) #Choice of Hemisphere  ccc=1 # to test

### Script to only do the gaussian weight making

require(ncdf4)  
require(spheRlab) 
require(tictoc) 

MASTERPATH="~/WORK/Data/SDAP/"
if(!dir.exists(MASTERPATH)) stop(sprintf("Data directory in MASTERPATH (%s) wasn't found. Please check",MASTERPATH))

if(ccc == 1) HEM="nh"
if(ccc== 2)  HEM="sh"

HEMPATH=paste0(MASTERPATH,HEM,"/")
dir.create(paste0(HEMPATH,"Outputs/"),recursive = T)


gridFilename=paste(HEMPATH,"Outputs/gridfile_OSISAF450",HEM,"_all_inclFlags",sep = "")
load(gridFilename)


tic("Gaussian weight making time:")
Rsphere=1  #Native grid units ()
gauss.sigma=2*pi*Rsphere/360
cutoff=1 #*gauss.sigma # One can change this later!

gsnwtfile=sprintf("%s/Outputs/gaussianweights%s_R%i_cutoff%i",HEMPATH,HEM,Rsphere,cutoff)

sf.gw.res=sl.spatialfilter.getweights(lon=grd$lon,lat=grd$lat,neighmat = grd$neighnodes,areas = grd$Nodeareas,Rsphere = 1,gauss.sigma = gauss.sigma,cutoff = cutoff*gauss.sigma)#for this specific grid, the weights are the same. it results in this list: sf.gw.res
save(file = gsnwtfile,sf.gw.res,version = 2)  
print(paste0("Gaussian weights for ",HEM," hemisphere saved: ",basename(gsnwtfile)))
toc()