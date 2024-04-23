###############################################################################
#Figures
################################################################################
setwd("./data/")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(lme4)
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
#  [1] vegan_2.6-4        lattice_0.21-8     permute_0.9-7      lme4_1.1-34        Matrix_1.6-1.1    
#[6] gridExtra_2.3      RColorBrewer_1.1-3 ggplot2_3.4.2     

#loaded via a namespace (and not attached):
#  [1] gtable_0.3.3      dplyr_1.1.3       compiler_4.3.1    tidyselect_1.2.0  Rcpp_1.0.11      
#[6] parallel_4.3.1    cluster_2.1.4     splines_4.3.1     scales_1.2.1      boot_1.3-28.1    
#[11] R6_2.5.1          generics_0.1.3    MASS_7.3-60       tibble_3.2.1      nloptr_2.0.3     
#[16] car_3.1-2         munsell_0.5.0     minqa_1.2.6       pillar_1.9.0      rlang_1.1.1      
#[21] utf8_1.2.4        cli_3.6.1         mgcv_1.8-42       withr_2.5.1       magrittr_2.0.3   
#[26] grid_4.3.1        rstudioapi_0.15.0 lifecycle_1.0.3   nlme_3.1-162      vctrs_0.6.4      
#[31] glue_1.6.2        abind_1.4-5       carData_3.0-5     fansi_1.0.5       colorspace_2.1-0 
#[36] tools_4.3.1       pkgconfig_2.0.3
######################## dataset with 29 neon sites #############################################
dat <- read.csv("./NEON28GSM_11172023.csv") #forest_type, lat, long, TAP, MAT column was added through excel
#################################################################################################
#Figure 1

A <- ggplot(dat, aes(x=forest_type, y=chm120sq)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sq of CHM (120 m)") + theme_light()
B <- ggplot(dat, aes(x=forest_type, y=chm120sbi)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sbi of CHM (120 m)") + theme_light()
C <- ggplot(dat, aes(x=forest_type, y=chm120sdr)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sdr of CHM (120 m)") + theme_light()

D <- ggplot(dat, aes(x=forest_type, y=Q25120sq)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sq of Q25 (120 m)") + theme_light()
E <- ggplot(dat, aes(x=forest_type, y=Q25120sbi)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sbi of Q25 (120 m)") + theme_light()
Fi <- ggplot(dat, aes(x=forest_type, y=Q25120sdr)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sdr of Q25 (120 m)") + theme_light()

G <- ggplot(dat, aes(x=forest_type, y=VCI120sq)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sq of VCI (120 m)") + theme_light()
H <- ggplot(dat, aes(x=forest_type, y=VCI120sbi)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sbi of VCI (120 m)") + theme_light()
I <- ggplot(dat, aes(x=forest_type, y=VCI120sdr)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sdr of VCI (120 m)") + theme_light()

J <- ggplot(dat, aes(x=forest_type, y=CVH120sq)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sq of CV(ht) (120 m)") + theme_light()
K <- ggplot(dat, aes(x=forest_type, y=CVH120sbi)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sbi of CV(ht) (120 m)") + theme_light()
L <- ggplot(dat, aes(x=forest_type, y=CVH120sdr)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="purple") + geom_boxplot(width=0.1) +
  labs(x="Forest type", y = "Sdr of CV(ht) (120 m)") + theme_light()

grid.arrange(A, B, C, D, E, Fi, G, H, I, J, K, L, nrow=4, ncol=3)

#############################################################################################
#Figure S3
A <- ggplot(dat, aes(x=site, y=chm120sq)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") +
  labs(x="Site", y = "Sq CHM") + coord_flip()
B <- ggplot(dat, aes(x=site, y=Q25120sq)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") + 
  labs(x="Site", y = "Sq Q25") + coord_flip()
C <- ggplot(dat, aes(x=site, y=VCI120sq)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") + 
  labs(x="Site", y = "Sq VCI") + coord_flip()
D <- ggplot(dat, aes(x=site, y=CVH120sq)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") + 
  labs(x="Site", y = "Sq CV(ht)") + coord_flip()

grid.arrange(A, B, C, D, nrow=2, ncol=2)

#############################################################################################
#Figure S4
A <- ggplot(dat, aes(x=site, y=chm120sbi)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") +
  labs(x="Site", y = "Sbi CHM") + coord_flip()
B <- ggplot(dat, aes(x=site, y=Q25120sbi)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") + 
  labs(x="Site", y = "Sbi Q25") + coord_flip()
C <- ggplot(dat, aes(x=site, y=VCI120sbi)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") + 
  labs(x="Site", y = "Sbi VCI") + coord_flip()
D <- ggplot(dat, aes(x=site, y=CVH120sbi)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") + 
  labs(x="Site", y = "Sbi CV(ht)") + coord_flip()

grid.arrange(A, B, C, D, nrow=2, ncol=2)

#############################################################################################
#Figure S5
A <- ggplot(dat, aes(x=site, y=chm120sdr)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") +
  labs(x="Site", y = "Sdr CHM") + coord_flip()
B <- ggplot(dat, aes(x=site, y=Q25120sdr)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") + 
  labs(x="Site", y = "Sdr Q25") + coord_flip()
C <- ggplot(dat, aes(x=site, y=VCI120sdr)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") + 
  labs(x="Site", y = "Sdr VCI") + coord_flip()
D <- ggplot(dat, aes(x=site, y=CVH120sdr)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="orange") + 
  labs(x="Site", y = "Sdr CV(ht)") + coord_flip()

grid.arrange(A, B, C, D, nrow=2, ncol=2)

#############################################################################################
#Figure 2
dat <- cbind(dat[,1:9], log((1+dat[,10:35])))
#standardize to get effect sizes
dat <- cbind(dat[,1:9], decostand(dat[,10:35], method  = "standardize"))

modelb <- glm(ndvi ~ chm120sq * forest_type, data = dat) 
B <- ggplot(dat, aes(x = chm120sq, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sq CHM", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modelc <- glm(ndvi ~  chm120sbi * forest_type, data = dat) 
C <- ggplot(dat, aes(x = chm120sbi, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sbi CHM", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modeld <- glm(ndvi ~  chm120sdr * forest_type, data = dat) 
D <- ggplot(dat, aes(x = chm120sdr, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sdr CHM", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modelf <- glm(ndvi ~  Q25120sq * forest_type, data = dat) 
Fi <- ggplot(dat, aes(x =  Q25120sq, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sq Q25", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modelg <- glm(ndvi ~  Q25120sbi * forest_type, data = dat) 
G <- ggplot(dat, aes(x = Q25120sbi, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sbi Q25", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modelh <- glm(ndvi ~  Q25120sdr * forest_type, data = dat) 
H <- ggplot(dat, aes(x = Q25120sdr, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sdr Q25", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modelj <- glm(ndvi ~ VCI120sq * forest_type, data = dat) 
J <- ggplot(dat, aes(x = VCI120sq, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sq VCI", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modelk <- glm(ndvi ~  VCI120sbi * forest_type, data = dat) 
K <- ggplot(dat, aes(x = VCI120sbi, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sbi VCI", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modell <- glm(ndvi ~  VCI120sdr * forest_type, data = dat) 
L <- ggplot(dat, aes(x = VCI120sdr, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sdr VCI", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modeln <- glm(ndvi ~  CVH120sq * forest_type, data = dat) 
N <- ggplot(dat, aes(x = CVH120sq, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sq CV(ht)", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

model0 <- glm(ndvi ~ CVH120sbi * forest_type, data = dat) 
O <- ggplot(dat, aes(x = CVH120sbi, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sbi CV(ht)", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

modelp <- glm(ndvi ~  CVH120sdr * forest_type, data = dat) 
P <- ggplot(dat, aes(x = CVH120sdr, y= ndvi, group=forest_type)) +
  geom_point(aes(colour = factor(forest_type)), size=1, shape=19) + theme(legend.position="none") +
  labs(x="Sdr CV(ht)", y="NDVI") + geom_point(alpha = .05) + geom_line(aes(y = predict(modell), colour = factor(forest_type)), size = 1) + 
  scale_colour_manual(values= c("#0D6BFF", "#0DFFFA", "#FF8E0C"))

grid.arrange(B, C, D, Fi, G, H, J, K, L, N, O, P, nrow=4, ncol=3)

