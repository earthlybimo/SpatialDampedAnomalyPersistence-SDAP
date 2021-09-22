# SpatialDampedAnomalyPersistence_SDAP
Script to generate the SpatialDampedAnomalyPersistence forecast from historical and present SIC data


# Required packages:
We need these two R packages to use the method:

ncdf4  #To open netcdf files

spheRlab #To generate contours, plots and other grid functions
(can be downloaded from here: https://github.com/FESOM/spheRlab)

tictoc #Used for timing but optional

# Data Directory
We will use one specific path for all our data operations. Please choose this appropriately to handle storage issues. Softlinks are also ok.
MASTERPATH="~/WORK/Data/SDAP/"

In each script, this path has been stated at the top and can be changed with the appropriate storage path. 

Within this masterdirectory, we will divide our outputs by hemisphere.
For example, HEM="sh" or HEM ="nh"
HEMPATH=paste0(MASTERPATH,"/",HEM)
dir.create(paste0(HEMPATH,"/Outputs/"),recursive = T)

All outputs will be saved accordingly in this folder, further divided into subfolders

# Input data
We need the input data to find initial ice-edge and climatological SIP. For this example, we will use the SIC records from OSI SAF. Other SIC records can also be used, although some troubleshooting might be required in further steps. OSI SAF data can be downloaded from here:

ftp://osisaf.met.no/reprocessed/ice/conc/v2p0 #OSI 450  (for years 1979 to 2015)

ftp://osisaf.met.no/reprocessed/ice/conc-cont-reproc/v2p0/ #OSI 430b (for years 2016 on)

Note that OSI SAF 450 was continued as OSI 430b after 2016, therefore we have 2 datapaths here. They can also be merged into one directory, but here we have left them in the original structure.

In the rest of the scripts, we will assume that the data has already been downloaded in these locations:

Datapath1=paste0(MASTERPATH,"/OSI450/osisaf.met.no") #OSI 450 
Datapath2=paste0(MASTERPATH,"/OSI430b/osisaf.met.no")#OSI 430b

