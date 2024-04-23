################################################################################
#Mask non-forest landcover from NDVI at each NEON site
################################################################################
library(rgdal)
library(raster)
setwd("./data/")

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
################################################################################
#run site by site
site <- ""
################################################################################
#list all files in each metric x site folder
output_folder <- paste0("./LULC/ndvi_forestmasked/")
ndvi <- raster(paste0("./ndvi/ndvi_rasters/", site, "_NDVI.tif")) #ndvi
maskLULC <- raster(paste0("./LULC/ESRI_reclass_resample_ndvi/", site, "ndvi_reclass30m.tif"))
ndvi <- projectRaster(ndvi, crs = crs(maskLULC))

r.i <- mask(ndvi, maskLULC, filename= paste0(output_folder, "mask_ndvi_", site, ".tif"), inverse=FALSE, 
               maskvalue=NA, updatevalue=NA, updateNA=FALSE)

#spot check, cuts to extent of land cover
plot(r.i)
################################################################################