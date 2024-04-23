################################################################################
#Resample 10 m land cover to match 30 m NDVI and clip to size of each NEON site
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
site <- "ABBY"
################################################################################
#mosaic the ESRI 2019 10 m land cover tiles together from over the USA in ArcMap
LULC <- raster(paste0("./LULC/ESRI10_clipped/", site, "/", site, ".tif"))

ndvi <- raster(paste0("./ndvi/ndvi_rasters/", site, "_NDVI.tif")) 
ndvi <- projectRaster(ndvi, crs = crs(LULC))
###############################################################################
################################################################################
#reclassify landcover to be forest (1) and non-forest (NA)
table(getValues(LULC))
#forest = 2
#non-forest = 1, 4, 5, 7, 8, 9, 10, 11
#create a reclassification matrix
reclass_df <- c(0, 1, NA, 
                2, 2, 1, 
                3, 12, NA) 
#create a matrix
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)

LULC <- reclassify(LULC, reclass_m)
writeRaster(LULC, file= paste0("./LULC/ESRI_reclass_ndvi/", site, "ndvi_reclass10m.tif"), format="GTiff", overwrite = TRUE)
################################################################################
#now resample the LULC 10m reclassified binary forest layer to 30m to match the ndvi of this site
LULC <- resample(LULC, ndvi, method = "ngb", filename = paste0("./LULC/ESRI_reclass_resample_ndvi/", site, "ndvi_reclass30m.tif"), overwrite = TRUE)

