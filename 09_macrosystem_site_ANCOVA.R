###############################################################################
#Table 5 - GSM-SD and NDVI
################################################################################
setwd("./data/")
library(car)
library(vegan)

#sessionInfo()
#R version 4.3.1 (2023-06-16)
#Platform: aarch64-apple-darwin20 (64-bit)
#Running under: macOS Ventura 13.6

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: America/Denver
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] vegan_2.6-4    lattice_0.21-8 permute_0.9-7  car_3.1-2      carData_3.0-5 

#loaded via a namespace (and not attached):
#  [1] MASS_7.3-60    compiler_4.3.1 Matrix_1.6-1.1 parallel_4.3.1 tools_4.3.1    mgcv_1.8-42    abind_1.4-5   
#[8] splines_4.3.1  nlme_3.1-162   grid_4.3.1     cluster_2.1.4
######################## dataset with 29 neon sites #############################################
dat <- read.csv("./NEON28GSM_11172023.csv") #forest_type, lat, long, TAP, MAT column was added through excel
dat$longitude <- dat$longitude * -1 #remove negatives from longitude
#log transform for structure and ndvi for linear relationships - not lat/long
dat <- cbind(dat[,c(1:5)], log((1+dat[,c(10:35)])))
#standardize to get effect sizes
dat <- cbind(dat[,1:3], decostand(dat[,4:31], method  = "standardize"))
#################################################################################################
# NDVI ~ sq120_chm
model <- glm(ndvi ~ chm120sq * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sbi120_chm
model <- glm(ndvi ~ chm120sbi * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sdr120_CHM
model <- glm(ndvi ~ chm120sdr * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sq120_Q25
model <- glm(ndvi ~ Q25120sq * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sbi120_Q25
model <- glm(ndvi ~ Q25120sbi * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sdr120_Q25
model <- glm(ndvi ~ Q25120sdr * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sq120_VCI
model <- glm(ndvi ~  VCI120sq * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sbi120_VCI
model <- glm(ndvi ~ VCI120sbi * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sdr120_VCI
model <- glm(ndvi ~ VCI120sdr * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sq120_CVht
model <- glm(ndvi ~ CVH120sq * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# NDVI ~ sbi120_CVht
model <- glm(ndvi ~ CVH120sbi * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 


# NDVI ~ sdr120_CVht
model <- glm(ndvi ~ CVH120sdr * forest_type + latitude + longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

#---------------------------------------------------------------------------------------------

#Table S1
###############
########################
###################################
dat2 <- dat[, c(1,6, 19:30)]
sites <- unique(dat$site)
varnames <- colnames(dat[,3:14])
OUT <- NULL

for(i in 1:length(sites)) {
  dat.i <- dat[dat$site == sites[i], ]
  name.i <- sites[i]
  
  for(n in 1:length(varnames)) {
  varnames.n <- varnames[n]
  model <- glm(dat.i$ndvi ~ dat.i[,varnames.n])
  a <- summary(model)$coefficients 
  a1 <- confint(model, method = "boot") 
  
  out.i <- c(name.i,varnames.n, 
             round(a[1,1], digits = 3), round(a1[1,1:2], digits = 3),
             round(a[2,1], digits = 3), round(a1[2,1:2], digits = 3))
  OUT <- rbind(OUT, out.i)
  }
}

colnames(OUT) <- c("site", "predictor", 
                   "intercept", "intercept_2.5", "intercept_97,5", 
                   "slope", "slope_2.5", "slope_97.5")

write.csv(OUT, "./aggregated/Individual_site_glms_12042023.csv")
