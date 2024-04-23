################################################################################
#Mask non-forest land cover from SD of each NEON site
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
siteEDI <- ""
site <- ""
################################################################################
#list all files in each metric x site folder
output_folder <- paste0("./LULC/SD_forestmasked/", site, "/")

#four metrics
metrics <- c(paste0("./LULC/EDIraster/", siteEDI, "/", "CSD_", siteEDI, "_Cutoff05-MCH.tif"), 
             paste0("./LULC/EDIraster/", siteEDI, "/", "CSD_", siteEDI, "_Cutoff05-Q25.tif"),
             paste0("./LULC/EDIraster/", siteEDI, "/", "CSD_", siteEDI, "_Cutoff05-VCI.tif"),
             paste0("./LULC/EDIraster/", siteEDI, "/", "CSD_", siteEDI, "_Cutoff05-CVH.tif"))

maskLULC <- rast(paste0("./LULC/ESRI_reclass_resample/",site, "/", site, "_reclass30m.tif"))

for(i in 1:length(metrics)){
  r.i <- rast(metrics[i])
  r.i.name <- substr(metrics[i], nchar(metrics[i]) - 6, nchar(metrics[i]))
  r.i <- mask(r.i, maskLULC, filename= paste0(output_folder, "mask_", site, "_", r.i.name), inverse=FALSE, 
               maskvalues=NA, updatevalue=NA, overwrite=TRUE)
}

#spot check some of them
test <- rast(paste0(output_folder, "mask_", site, "_", r.i.name))
plot(test)
################################################################################
