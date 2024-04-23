###############################################################################
#Table 4 - Macroscale predictors of GSM-SD
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
 # [1] MASS_7.3-60    compiler_4.3.1 Matrix_1.6-1.1 parallel_4.3.1 tools_4.3.1    mgcv_1.8-42    abind_1.4-5   
#[8] splines_4.3.1  nlme_3.1-162   grid_4.3.1     cluster_2.1.4 
######################## dataset with 29 neon sites #############################################
dat <- read.csv("./NEON28GSM_11172023.csv") #forest_type, lat, long, TAP, MAT column was added through excel
dat$longitude <- dat$longitude * -1 #remove negatives from longitude
dat$mat <- dat$mat + 3
#log transform for structure and ndvi for linear relationships
dat <- cbind(dat[,c(1:3)], log((1+dat[,c(4:35)])))
#standardize to get effect sizes
dat <- cbind(dat[,1:3], decostand(dat[,4:35], method  = "standardize"))

#################################################################################################
# sq120_chm ---------------------------------------------
model <- glm(chm120sq ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sq ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sq ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sq ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sq ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# sbi120_chm ---------------------------------------------
model <- glm(chm120sbi ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sbi ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sbi ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sbi ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sbi ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# sdr120_chm ---------------------------------------------
model <- glm(chm120sdr ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sdr ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sdr ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sdr ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(chm120sdr ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# sq120_Q25 ---------------------------------------------
model <- glm(Q25120sq ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sq ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sq ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sq ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sq ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# sbi120_Q25 ---------------------------------------------
model <- glm(Q25120sbi ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sbi ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sbi ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sbi ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sbi ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# sdri120_Q25 ---------------------------------------------
model <- glm(Q25120sdr ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sdr ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sdr ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sdr ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(Q25120sdr ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# sq120_VCI ---------------------------------------------
model <- glm(VCI120sq ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sq ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sq ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sq ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sq ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# sbi120_VCI ---------------------------------------------
model <- glm(VCI120sbi ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sbi ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sbi ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sbi ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sbi ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot")              

# sdr120_VCI ---------------------------------------------
model <- glm(VCI120sdr ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sdr ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sdr ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sdr ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(VCI120sdr ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# sq120_CVht ---------------------------------------------
model <- glm(CVH120sq ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sq ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sq ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sq ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sq ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

# sbi120_CVht ---------------------------------------------
model <- glm(CVH120sbi ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sbi ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sbi ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sbi ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sbi ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot")

# sdr120_CVht ---------------------------------------------
model <- glm(CVH120sdr ~ forest_type, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sdr ~ latitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sdr ~ longitude, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sdr ~ mat, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 

model <- glm(CVH120sdr ~ tap, data = dat)
qqnorm(residuals(model))
car::Anova(model, type = "III")
confint(model, method = "boot") 
