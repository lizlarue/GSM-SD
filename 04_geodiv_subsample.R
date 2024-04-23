################################################################################
#Random subsampled points at NEON site raster to get geodiv and outputs as a csv
################################################################################
setwd("./data/")
################################################################################
library(terra)
library(geodiv) #sbi (dominance), sdr (edge density), sbi (heterogeneity)
library(sf)
library(dplyr)

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
#  [1] dplyr_1.1.3  sf_1.0-14    geodiv_1.1.0 terra_1.7-46

#loaded via a namespace (and not attached):
#  [1] vctrs_0.6.4        cli_3.6.1          rlang_1.1.1        DBI_1.1.3         
#[5] KernSmooth_2.23-21 generics_0.1.3     zoo_1.8-12         glue_1.6.2        
#[9] e1071_1.7-13       pracma_2.4.2       fansi_1.0.5        grid_4.3.1        
#[13] classInt_0.4-10    tibble_3.2.1       lifecycle_1.0.3    compiler_4.3.1    
#[17] codetools_0.2-19   Rcpp_1.0.11        pkgconfig_2.0.3    lattice_0.21-8    
#[21] R6_2.5.1           class_7.3-22       tidyselect_1.2.0   utf8_1.2.3        
#[25] pillar_1.9.0       parallel_4.3.1     magrittr_2.0.3     tools_4.3.1       
#[29] proxy_0.4-27       spatial_7.3-16     units_0.8-4   
###############################################################################
#change site name here for each site
site <- "YELL"
###############################################################################
#list all files in each metric x site folder
mosaic_folder <- paste0("./LULC/SD_forestmasked/", site, "/")
surface_folder <- paste0("./LULC/geodiv/", site, "/")
###############################################################################
#30 x 30 m but cannot stack to match to run geodiv functions
ndvi <- rast(paste0("./LULC/ndvi_forestmasked/mask_ndvi_", site, ".tif")) 
CVH <- rast(paste0(mosaic_folder, "mask_",site, "_CVH", ".tif"))
VCI <- rast(paste0(mosaic_folder, "mask_",site, "_VCI", ".tif"))
Q25 <-rast(paste0(mosaic_folder,"mask_", site, "_Q25", ".tif"))
chm <-rast(paste0(mosaic_folder,"mask_", site, "_MCH", ".tif"))

###############################################################################
#randomly subsample coords from site extent
#dimensions of the raster
site_extent <- ext(chm)[1:4] #xmin, xmax, ymin, ymax

buffer <- 500 #remove a 0.5 km buffer around the edge to avoid areas that can't be sampled for neighborhoods
#randomly select a subset of cells
randx <- round(runif(900, min = site_extent[1] + buffer , max = site_extent[2] - buffer), digits = 0)
randy <- round(runif(900, min = site_extent[3] + buffer, max = site_extent[4] - buffer), digits = 0)
randxy <- cbind(randx, randy)
randxy <- randxy[!duplicated(randxy),] #remove any duplicate cells
#randxy <- randxy[1:300,] 
colnames(randxy) <- c("x", "y")
randxy <- as.data.frame(randxy)
plot(randxy$x, randxy$y)
###############################################################################
#extract data from cells of random coord after clipping out 120 and 60 m square windows
all_metrics <- NULL
for(r in 1:nrow(randxy)){
  
  #turn coords into a spatial object with match crs
   cords_sf <- randxy[r,] %>% 
    st_as_sf(coords = c("x", "y"), dim = "XY") %>% 
    st_set_crs(crs(chm)) %>% 
    select()
  
  #120 m window size ------------------------------------------------------------------
  
  # Crop the raster to a 120 m grid
  extent120 <- st_buffer(cords_sf, dist = 60, endCapStyle = "SQUARE")
  chm120 <- crop(chm, extent120)
  Q25120 <- crop(Q25, extent120)
  CVH120 <- crop(CVH, extent120)
  VCI120 <- crop(VCI, extent120)
  
  if(anyNA(values(chm120))) {
    # Skip if NA values
    next
  }
     
  ndvi_r <- extract(ndvi, extent120, fun = mean, method = 'simple')
  
  # Root mean-square Roughness (sq)
  chm120sq <- sq(chm120)
  Q25120sq <- sq(Q25120)
  CVH120sq <- sq(CVH120)
  VCI120sq <- sq(VCI120)
  
  #Surface Bearing Index (sbi)
  chm120sbi <- sbi(chm120)
  Q25120sbi <- sbi(Q25120)
  CVH120sbi <- sbi(CVH120)
  VCI120sbi <- sbi(VCI120)
  
  # Surface Area Ratio (sdr)
  chm120sdr <- sdr(chm120)
  Q25120sdr <- sdr(Q25120)
  CVH120sdr <- sdr(CVH120)
  VCI120sdr <- sdr(VCI120)
  # Create texture image of chosen metric with square moving window with 60 m ---------------------
  
  # Crop the raster to a 60 m grid
  extent60 <- st_buffer(cords_sf, dist = 30, endCapStyle = "SQUARE")
  ndvi_60 <- extract(ndvi, extent60, fun = mean, method = 'simple')
  chm60 <- crop(chm, extent60)
  Q2560 <- crop(Q25, extent60)
  CVH60 <- crop(CVH, extent60)
  VCI60 <- crop(VCI, extent60)
  
  #Root mean-square Roughness (sq)
  chm60sq <- sq(chm60)
  Q2560sq <- sq(Q2560)
  CVH60sq <- sq(CVH60)
  VCI60sq <- sq(VCI60)
  
  #Surface Bearing Index (sbi)
  chm60sbi <- sbi(chm60)
  Q2560sbi <- sbi(Q2560)
  CVH60sbi <- sbi(CVH60)
  VCI60sbi <- sbi(VCI60)
  
  # Surface Area Ratio (sdr)
  chm60sdr <- sdr(chm60)
  Q2560sdr <- sdr(Q2560)
  CVH60sdr <- sdr(CVH60)
  VCI60sdr <- sdr(VCI60)
  
  all_metrics_r <- c(randxy[r,1], randxy[r,2], ndvi_r[,2], 
                     chm60sq, Q2560sq, CVH60sq, VCI60sq,
                     chm60sbi, Q2560sbi, CVH60sbi, VCI60sbi,
                     chm60sdr,Q2560sdr, CVH60sdr, VCI60sdr,
                     chm120sq, Q25120sq, CVH120sq, VCI120sq,
                    chm120sbi, Q25120sbi, CVH120sbi, VCI120sbi,
                    chm120sdr, Q25120sdr, CVH120sdr, VCI120sdr, ndvi_60[,2])
  all_metrics <- rbind(all_metrics, all_metrics_r)
}


all_metrics <- cbind.data.frame(rep(site, nrow(all_metrics)), all_metrics)
colnames(all_metrics) <- c("site", "easting","northing", "ndvi", 
                           "chm60sq", "Q2560sq", "CVH60sq", "VCI60sq",
                           "chm60sbi", "Q2560sbi", "CVH60sbi", "VCI60sbi",
                           "chm60sdr","Q2560sdr", "CVH60sdr", "VCI60sdr",
                           "chm120sq", "Q25120sq", "CVH120sq", "VCI120sq",
                           "chm120sbi", "Q25120sbi", "CVH120sbi", "VCI120sbi",
                           "chm120sdr", "Q25120sdr", "CVH120sdr", "VCI120sdr", "ndvi_60")

all_metrics <- all_metrics[!is.na(all_metrics$ndvi),]
all_metrics <- all_metrics[!is.na(all_metrics$CVH),]
all_metrics <- all_metrics[1:120, ]
#save each site separately
write.csv(all_metrics, paste0("./LULC/extracted/", site, "/", site, "_rand_points_11142023.csv"))