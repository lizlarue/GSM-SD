################################################################################
#Create a panel of example lower and higher values of GSM of structural diversity
#This is the new Fig. 1 in the revision
################################################################################
setwd("./data/") 
################################################################################
library(terra)
library(geodiv) 
library(sf)
library(dplyr)
library(neonUtilities)
library(scales)
library(lidR)
#sessionInfo()
###############################################################################
#Selected site is BART with relative low and high set of GSM-SD values for CHM from data set
#30 x 30 m masked CHM (not GSM-SD yet)
chm <-rast("./Fig1BART/mask_BART_MCH.tif") 

##############################################################################
#pull raw lidar to make figure
byTileAOP(dpID = "DP1.30003.001", site = "BART", year = "2019", #can only do 1 year at a time
          easting = c(320112, 317654, 318210, 322531, 320961, 318019), 
          northing = c(4880706, 4880775, 4878152, 4875773, 4876645, 4879949), buffer = 150, 
          savepath="./Fig1BART/", include.provisional = FALSE)

##############################################################################
#GSM of SQ ---------------------- heterogeneity
low_sq <- c(320112, 4880706) #easting, northing
low_sq <- as.data.frame(matrix(low_sq, nrow = 1, ncol = 2))
colnames(low_sq) <- c("easting", "northing")
low_sq <- low_sq %>% 
  st_as_sf(coords = c("easting", "northing"), dim = "XY") %>% 
  st_set_crs(st_crs(chm)) %>% 
  select()

high_sq <- c(317654, 4880775)
high_sq <- as.data.frame(matrix(high_sq, nrow = 1, ncol = 2))
colnames(high_sq) <- c("easting", "northing")
high_sq <- high_sq %>% 
  st_as_sf(coords = c("easting", "northing"), dim = "XY") %>% 
  st_set_crs(st_crs(chm)) %>% 
  select()

#GSM of SBI ---------------------- Dominance
low_sbi <- c(318210, 4878152)
low_sbi <- as.data.frame(matrix(low_sbi, nrow = 1, ncol = 2))
colnames(low_sbi) <- c("easting", "northing")
low_sbi <- low_sbi %>% 
  st_as_sf(coords = c("easting", "northing"), dim = "XY") %>% 
  st_set_crs(st_crs(chm)) %>% 
  select()

high_sbi <- c(322531, 4875773)
high_sbi <- as.data.frame(matrix(high_sbi, nrow = 1, ncol = 2))
colnames(high_sbi) <- c("easting", "northing")
high_sbi <- high_sbi %>% 
  st_as_sf(coords = c("easting", "northing"), dim = "XY") %>% 
  st_set_crs(st_crs(chm)) %>% 
  select()

#GSM of SDR ---------------------- Edge density
low_sdr <- c(320961, 4876645)
low_sdr <- as.data.frame(matrix(low_sdr, nrow = 1, ncol = 2))
colnames(low_sdr) <- c("easting", "northing")
low_sdr <- low_sdr %>% 
  st_as_sf(coords = c("easting", "northing"), dim = "XY") %>% 
  st_set_crs(st_crs(chm)) %>% 
  select()

high_sdr <- c(318019, 4879949)#316887, 4880796)
high_sdr <- as.data.frame(matrix(high_sdr, nrow = 1, ncol = 2))
colnames(high_sdr) <- c("easting", "northing")
high_sdr <- high_sdr %>% 
  st_as_sf(coords = c("easting", "northing"), dim = "XY") %>% 
  st_set_crs(st_crs(chm)) %>% 
  select()
###############################################################################
  #120 m window size ------------------------------------------------------------------
  
  # Crop the raster to a 120 m grid for each metric and run the GSM 
  # Root mean-square Roughness (sq)
  extent120low_sq <- st_buffer(low_sq, dist = 60, endCapStyle = "SQUARE")
  chm_lower_sq <- crop(chm, extent120low_sq) 
  round(sq(chm_lower_sq), digits = 2)
  
  extent120high_sq <- st_buffer(high_sq, dist = 60, endCapStyle = "SQUARE")
  chm_high_sq <- crop(chm, extent120high_sq) 
  round(sq(chm_high_sq), digits = 2)
  
  # -------------------------------------------------------------
  #Surface Bearing Index (sbi)
  extent120low_sbi <- st_buffer(low_sbi, dist = 60, endCapStyle = "SQUARE")
  chm_lower_sbi <- crop(chm, extent120low_sbi) 
  round(sbi(chm_lower_sbi), digits = 2)
  
  extent120high_sbi <- st_buffer(high_sbi, dist = 60, endCapStyle = "SQUARE")
  chm_high_sbi <- crop(chm, extent120high_sbi) 
  round(sbi(chm_high_sbi), digits = 2)
  
  # -------------------------------------------------------------
  # Surface Area Ratio (sdr)
  extent120low_sdr <- st_buffer(low_sdr, dist = 60, endCapStyle = "SQUARE")
  chm_lower_sdr <- crop(chm, extent120low_sdr) 
  round(sdr(chm_lower_sdr), digits = 2)
  
  extent120high_sdr <- st_buffer(high_sdr, dist = 60, endCapStyle = "SQUARE")
  chm_high_sdr <- crop(chm, extent120high_sdr) 
  round(sdr(chm_high_sdr), digits = 2)

#-----------------------------------------------------------------
forest <- gradient_n_pal(c("#0D6BFF","#0DFFFA","#FF8E0C"))(seq(0, 1, length.out = 100))
par(mfrow = c(2,3))
plot(chm_lower_sq, col = forest, range = c(1.44, 21.86), axes = FALSE)
plot(chm_lower_sbi, col = forest, range = c(1.44, 21.86), axes = FALSE)
plot(chm_lower_sdr, col = forest, range = c(1.44, 21.86), axes = FALSE)

plot(chm_high_sq, col = forest, range = c(1.44, 21.86), axes = FALSE)
plot(chm_high_sbi, col = forest, range = c(1.44, 21.86), axes = FALSE)
plot(chm_high_sdr, col = forest, range = c(1.44, 21.86), axes = FALSE)#, legend = FALSE)  

#-----------------------------------------------------------------
#extract raw lidar point cloud profiles
ctg <- readLAScatalog(folder = "./DP1.30003.001/neon-aop-products/2019/FullSite/D01/2019_BART_5/L1/DiscreteLidar/ClassifiedPointCloud/")
las_check(ctg)
par(mfrow = c(1,1))
plot(ctg)
#-------------------------------------------------
lowsq <- clip_roi(ctg, extent120low_sq)
lowsq_denoise <- classify_noise(lowsq, ivf(res=2,n=3)) 
lowsq_denoise <- filter_poi(lowsq_denoise, Classification != LASNOISE)
lowsq_denoise <- normalize_height(lowsq_denoise, tin())
plot(lowsq_denoise, color = "RGB", bg = "white")

highsq <- clip_roi(ctg, extent120high_sq)
highsq_denoise <- classify_noise(highsq, ivf(res=2,n=3)) 
highsq_denoise <- filter_poi(highsq_denoise, Classification != LASNOISE)
highsq_denoise <- normalize_height(highsq_denoise, tin())
plot(highsq_denoise, color = "RGB", bg = "white")
#-------------------------------------------------
lowsbi <- clip_roi(ctg, extent120low_sbi)
lowsbi_denoise <- classify_noise(lowsbi, ivf(res=2,n=3)) 
lowsbi_denoise <- filter_poi(lowsbi_denoise, Classification != LASNOISE)
lowsbi_denoise <- normalize_height(lowsbi_denoise, tin())
plot(lowsbi_denoise, color = "RGB", bg = "white")

highsbi <- clip_roi(ctg, extent120high_sbi)
highsbi_denoise <- classify_noise(highsbi, ivf(res=2,n=3)) 
highsbi_denoise <- filter_poi(highsbi_denoise, Classification != LASNOISE)
highsbi_denoise <- normalize_height(highsbi_denoise, tin())
plot(highsbi_denoise, color = "RGB", bg = "white")

#-------------------------------------------------
lowsdr <- clip_roi(ctg, extent120low_sdr)
lowsdr_denoise <- classify_noise(lowsdr, ivf(res=2,n=3)) 
lowsdr_denoise <- filter_poi(lowsdr_denoise, Classification != LASNOISE)
lowsdr_denoise <- normalize_height(lowsdr_denoise, tin())
plot(lowsdr_denoise, color = "RGB", bg = "white")

highsdr <- clip_roi(ctg, extent120high_sdr)
highsdr_denoise <- classify_noise(highsdr, ivf(res=2,n=3)) 
highsdr_denoise <- filter_poi(highsdr_denoise, Classification != LASNOISE)
highsdr_denoise <- normalize_height(highsdr_denoise, tin())
plot(highsdr_denoise, color = "RGB", bg = "white")
