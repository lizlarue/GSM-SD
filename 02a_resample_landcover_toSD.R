################################################################################
#Resample 10 m land cover to match 30 m SD and clip to size of each NEON site
################################################################################
library(terra)
setwd("./data/")

#sessionInfo()
#R version 4.3.1 (2023-06-16 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows Server 2019 x64 (build 17763)

#Matrix products: default


#locale:
#  [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
#[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.1252    

#time zone: America/Denver
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] terra_1.7-46

#loaded via a namespace (and not attached):
#  [1] compiler_4.3.1   tools_4.3.1      Rcpp_1.0.11      codetools_0.2-19
################################################################################
#run site by site
siteEDI <- "2019_YELL_2" #name of folder downloaded from EDI by Wang et al. (2023) that holds SD
site <- "YELL"
################################################################################
chm <- rast(paste0("./LULC/EDIraster/", siteEDI, "/", "CSD_", siteEDI, "_Cutoff05-MCH.tif"))
#mosaic the ESRI 2019 10 m land cover tiles together from over the USA in ArcMap
LULC <- rast(paste0("./LULC/ESRI10_clipped/", site, "/", site, ".tif"))
LULC <- project(LULC, crs(chm))
###############################################################################
#going to have to deal with different extents b/c LULC is square,
#whereas lidar is irregular
ext(chm)
ext(LULC)
################################################################################
#reclassify landcover to be forest (1) and non-forest (0)
#see which land use classes we have from the 2019 ESRI 10 m LULC legend
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

LULC <- classify(LULC, reclass_m)
plot(LULC)

#save
writeRaster(LULC, file= paste0("./LULC/ESRI_reclass/",site, "/", site, "_reclass10m.tif"), overwrite = TRUE)
################################################################################
#now resample the LULC 10m reclassified binary forest layer to 30m to match the chm of this site
#it makes the LULC fit the chm extent
LULC <- resample(LULC, chm, method = "near", filename = paste0("./LULC/ESRI_reclass_resample/",site, "/", site, "_reclass30m.tif"))
plot(LULC)
#a good test is to stack them and if so all good for masking step
test <- c(LULC, chm)

