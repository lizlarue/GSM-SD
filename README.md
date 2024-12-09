# GSM-SD
GSM-SD Manuscript
Project title: Gradient surface metrics of ecosystem structural diversity and their relationship with productivity across macrosystems

Authors: Elizabeth A. LaRue (ealarue@utep.edu), Kylie M. Rezendes, Dennis Choi, Jianmin Wang, Anna G. Downing, Songlin Fei, Brady S. Hardiman


About the data:
Structural diversity metrics were obtained from Wang et al. 2023), which were generated from NEON AOP LiDAR. Landsat 8 Collection 2 - Level 2 data were obtained from USGS Earth Explorer or from the Microsoft Planetary Computer for the download dates specified in "Landsat8DownloadDates.xlsx". The ESRI 10 m land cover is freely available from the Planetary Computer. Randomly sampled points and GSM of structural diversity metrics and NDVI used in analyses are provided in the data folder in file "NEON28GSM_11172023.csv". All raw data are freely available for download from the sources cited in the manuscript or listed above for the 28 NEON sites used. 

Wang, J., Choi, D., LaRue, E., Atkins, J., Foster, J., Matthes Hatala, J., Fahey, R., Fei, S., Hardiman, B. (2023). Structural Diversity from the NEON Discrete-Return LiDAR Point Cloud in 2013-2022. Environmental Data Initiative. https://doi.org/10.6073/pasta/e02f855d69193a46571168575b35291d

About the code steps: -----------------------------
01_NDVI.R - Generates NDVI and employs a cloud mask to Landsat 8 downloads. 

02a_resample_landcover_toSD.R; 02b_resample_landcover_toNDVI.R - Resamples ESRI 10 m landcover to match 30 m resolution of other datasets and turns non-forest cover and forest cover into a binary raster.

03a_maskSDtoforestcover.R; 03b_maskNDVItoforestcover.R - Masks for non-forest cover from structural diversity and NDVI rasters.

04_geodiv_subsample.R - Selects a random sample of coordinates within each NEON site and then runs GSMs and extracts NDVI from 60 and 120 m areas around the points. 

05_combine_data.R - Combines and cleans up data into one dataset used for the rest of the analyses "NEON28GSM_11172023.csv". 

06_spearman.R - Creates Fig. S2 output.

07_figures.R - Creates Fig. 1 and 2 (now Fig. 2 and 3), and Fig. S3-S5 output. 

08_GSM_predictors.R - Creates Table 4 output. 

09_macrosystem_site_ANCOVA.R - Creates Table 5 and Table S1 output. 

10_ExamplePanel.R - Creates a new figure of example GSM-SD rasters. Now Fig. 1

Variable metadata for "NEON28GSM_11172023.csv" for manuscript analysis: -----------------------------
site - NEON site acronym (standard use id for NEON)	
Index - Index variable from step 04
forest_type - Forest type	 
latitude - Latitude (N) of the NEON site (common coordinate for whole site). WGS1984 
longitude - Longitude (W) of the NEON site (common coordinate for whole site). WGS1984
tap - Total annual precipitation (mm)	
mat - Mean annual temperature (Celsius)	
easting	- Easting (UTM) of subsampled points within each site. Follows UTM zone used by NEON for their sites. 
northing - Northing (UTM) of subsampled points within each site. Follows UTM zone used by NEON for their sites. 	
ndvi - Normalized difference vegetation index 120 m window around subsampled points.	
chm60sq	- Root mean-square Roughness of mean outer canopy height 60 m window around subsampled points
Q2560sq	- Root mean-square Roughness of 25th quantile of canopy height 60 m window around subsampled points
CVH60sq	- Root mean-square Roughness of coefficient of variation of canopy height 60 m window around subsampled points
VCI60sq	- Root mean-square Roughness of vegetation complexity index 60 m window around subsampled points
chm60sbi - Surface Bearing Index of mean outer canopy height 60 m window around subsampled points 	
Q2560sbi - Surface Bearing Index of 25th quantile of canopy height 60 m window around subsampled points	 
CVH60sbi - Surface Bearing Index of coefficient of variation of canopy height 60 m window around subsampled points	
VCI60sbi - Surface Bearing Index of vegetation complexity index	60 m window around subsampled points
chm60sdr - Surface Area Ratio of mean outer canopy height 60 m window around subsampled points
Q2560sdr - Surface Area Ratio of 25th quantile of canopy height	60 m window around subsampled points
CVH60sdr - Surface Area Ratio of coefficient of variation of canopy height 60 m window around subsampled points	
VCI60sdr - Surface Area Ratio of vegetation complexity index 60 m window around subsampled points	
chm120sq - Root mean-square Roughness of mean outer canopy height 120 m window around subsampled points	
Q25120sq - Root mean-square Roughness of 25th quantile of canopy height	120 m window around subsampled points
CVH120sq - Root mean-square Roughness of coefficient of variation of canopy height 120 m window around subsampled points	
VCI120sq - Root mean-square Roughness of vegetation complexity index 120 m window around subsampled points	
chm120sbi - Surface Bearing Index of mean outer canopy height 120 m window around subsampled points	
Q25120sbi - Surface Bearing Index of 25th quantile of canopy height 120 m window around subsampled points	
CVH120sbi - Surface Bearing Index of coefficient of variation of canopy height 120 m window around subsampled points	
VCI120sbi - Surface Bearing Index of vegetation complexity index 120 m window around subsampled points	
chm120sdr - Surface Area Ratio of mean outer canopy height 120 m window around subsampled points	
Q25120sdr - Surface Area Ratio of 25th quantile of canopy height 120 m window around subsampled points	
CVH120sdr - Surface Area Ratio of coefficient of variation of canopy height 120 m window around subsampled points	
VCI120sdr - Surface Area Ratio of vegetation complexity index 120 m window around subsampled points	
ndvi_60	- - Normalized difference vegetation index 120 m window around subsampled points
index - Index variable for step 05. 
