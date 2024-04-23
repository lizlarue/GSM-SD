################################################################################
#Calculate NDVI and remove clouds from Landsat 8 Collection 2 - Level 2 data
################################################################################
library(rgdal)
library(raster)
#sessionInfo()
#R version 4.2.1 (2022-06-23 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows Server x64 (build 17763)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
#[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.1252    

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] raster_3.6-3 rgdal_1.5-32 sp_1.5-0    

#loaded via a namespace (and not attached):
#  [1] compiler_4.2.1   tools_4.2.1      Rcpp_1.0.9       codetools_0.2-18
#[5] grid_4.2.1       lattice_0.20-45  terra_1.6-17 

#Change last folder name of WD to reflect study site folder
setwd("./data/")
site <- "ABBY"

chm <- raster(paste0("./LULC/SD_forestmasked/", site, "/", "mask_", site, "_chm.tif"))
#chm <- projectRaster(chm, crs = (red)) #use if need to make chm match Landsat

#####################################################################################
#downloaded data is already uploaded to Azure, run for one site at a time
#####################################################################################
#import data
red = raster(list.files(path = paste0("./ndvi/raw_processed/", site, "/", site, "data/"), pattern = "\\B4.TIF$", full.names = TRUE)) #band 4
nir = raster(list.files(path = paste0("./ndvi/raw_processed/", site, "/", site, "data/"), pattern = "\\B5.TIF$", full.names = TRUE)) #band 5
QA = raster(list.files(path = paste0("./ndvi/raw_processed/", site, "/", site, "data/"), pattern = "\\PIXEL.TIF$", full.names = TRUE)) #QA_pixel 3 = clouds

#crop using chm of each site
red = crop(red, chm)
nir = crop(nir, chm)
QA = crop(QA, chm)

#for red, nir make it so that 0's = NA (0 = no data)
m <- matrix(c(0, 0, NA), ncol=3, byrow=TRUE)
red <- reclassify(red, m) 
nir <- reclassify(nir, m)

#ndvi
ndvi = (nir - red) / (nir + red)

#prepare QA_PIXEL for cloud mask from https://pages.cms.hu-berlin.de/EOL/gcg_eo/02_data_quality.html
#first must extract bit info
#bit 3 (position 4) = 1 = clouds are true
# Define function to find fill values from Landsat QA_PIXEL
fill_pixels <- function(x){
  return(intToBits(x)[4] == T) #position 4 or (bit 3) that equals true or yes cloud cover
}

#For cfmask and NDVI change the study site folder to save in appropriate space
#1 = the cloud values and 0 = non-cloud
cfmask <- calc(QA, fill_pixels, filename = paste0("./ndvi/raw_processed/", site, "/Output_Rasters/cfmask.tif"), overwrite = TRUE)
ndvi <- mask(ndvi, cfmask, filename= paste0("./ndvi/raw_processed/", site, "/Output_Rasters/ndvi.tif"), inverse=FALSE, 
             maskvalue=1, updatevalue=NA, updateNA=FALSE, overwrite=TRUE)

ndvi <- mask(ndvi, cfmask, filename= paste0("./ndvi/NDVI_rasters/", site, "_NDVI.tif"), inverse=FALSE, 
             maskvalue=1, updatevalue=NA, updateNA=FALSE, overwrite=TRUE)

