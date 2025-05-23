# SpatialDampedAnomalyPersistence_SDAP
Script to generate the SpatialDampedAnomalyPersistence forecast from historical and present SIC data


# Required packages:
We need two R packages to use the method:

ncdf4  #To open netcdf files

spheRlab #To generate contours, plots and other grid functions
(can be downloaded from here: https://github.com/FESOM/spheRlab)

Optionally, it is also useful to have this package for timing the operations:
tictoc 

# Input data
We need the input data to find initial ice-edge and climatological SIP. For this example, we will use the SIC records from OSI SAF. Other SIC records can also be used, although some troubleshooting might be required in further steps. OSI SAF data can be downloaded from here:

ftp://osisaf.met.no/reprocessed/ice/conc/v2p0 #OSI 450  (for years 1979 to 2015)

ftp://osisaf.met.no/reprocessed/ice/conc-cont-reproc/v2p0/ #OSI 430b (for years 2016 on)

Note that OSI SAF 450 was continued as OSI 430b after 2016, therefore we have 2 datapaths here. They can also be merged into one directory, but here we have left them in the original structure.

In the rest of the scripts, we will assume that the data has already been downloaded in these locations:

Datapath1=paste0(MASTERPATH,"/OSI450/osisaf.met.no") #OSI 450 
Datapath2=paste0(MASTERPATH,"/OSI430b/osisaf.met.no")#OSI 430b


# Output data directory

We will use one specific path for all our data operations. Please choose this appropriately to handle storage issues. Softlinks are also ok.
MASTERPATH="~/WORK/Data/SDAP/"

In each script, this path has been stated at the top and can be changed with the appropriate storage path. 

Within this masterdirectory, we will divide our outputs by hemisphere.
For example, HEM="sh" or HEM ="nh"
HEMPATH=paste0(MASTERPATH,"/",HEM)
dir.create(paste0(HEMPATH,"/Outputs/"),recursive = T)

All outputs will be saved accordingly in this folder, further divided into subfolders.

Here, all outputs are saved in native R data format. This can be opened with the R software package and saved in other formats as necessary. Actual forecasts, saved in netcdf format, can be found here:

For the Arctic - https://zenodo.org/record/5543773
For the Antarctic - https://zenodo.org/record/5543957

These files can be opened using ncview or other netcdf enabled software, as well as R, Python, Matlab or other scientific computing language.


---------------------
Sample files from the scripts are included in the sample directory. However, they are not arranged properly into directories as expected by the script, but are to test individual steps of the process. The sample DampedForecast file (as would be generated by script S5 and used by script S6 for plotting SPS or forecast maps) is too large and could not be uploaded to github. It can be generated using script S5 or copied from here: https://drive.google.com/file/d/1zhuR26VO1Kt2esx9b2HAfIDz7LD7CfSz/view?usp=sharing/view    

In case of any issue with the script, codes or the files, please send an email to bniraula@awi.de for clarification. 

![lotus](https://github.com/user-attachments/assets/dc226cce-c620-4141-a3ef-ce8ea2eb10c4)

