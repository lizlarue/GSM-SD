###############################################################################
#spearman correlations
################################################################################
setwd("./data/")
library("PerformanceAnalytics")
library(ggplot2)
library(ggcorrplot)
library(Hmisc)
library(stringr)
library(tidyr)
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
#  [1] vegan_2.6-4                lattice_0.21-8             permute_0.9-7             
#[4] tidyr_1.3.0                stringr_1.5.0              Hmisc_5.1-0               
#[7] ggcorrplot_0.1.4           ggplot2_3.4.2              PerformanceAnalytics_2.0.4
#[10] xts_0.13.1                 zoo_1.8-12                

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.4        generics_0.1.3    stringi_1.7.12    digest_0.6.33     magrittr_2.0.3   
#[6] evaluate_0.21     grid_4.3.1        fastmap_1.1.1     Matrix_1.6-1.1    nnet_7.3-19      
#[11] backports_1.4.1   Formula_1.2-5     gridExtra_2.3     mgcv_1.8-42       purrr_1.0.1      
#[16] fansi_1.0.5       scales_1.2.1      cli_3.6.1         rlang_1.1.1       splines_4.3.1    
#[21] munsell_0.5.0     base64enc_0.1-3   withr_2.5.1       parallel_4.3.1    tools_4.3.1      
#[26] checkmate_2.3.0   htmlTable_2.4.1   dplyr_1.1.3       colorspace_2.1-0  vctrs_0.6.4      
#[31] R6_2.5.1          rpart_4.1.19      lifecycle_1.0.3   htmlwidgets_1.6.2 MASS_7.3-60      
#[36] foreign_0.8-84    cluster_2.1.4     pkgconfig_2.0.3   pillar_1.9.0      gtable_0.3.3     
#[41] glue_1.6.2        data.table_1.14.8 xfun_0.39         tibble_3.2.1      tidyselect_1.2.0 
#[46] rstudioapi_0.15.0 knitr_1.43        nlme_3.1-162      htmltools_0.5.5   rmarkdown_2.23   
#[51] compiler_4.3.1    quadprog_1.5-8 
#################################################################################################
dat <- read.csv("./NEON28GSM_11172023.csv")
#log transform for linear relationships
dat <- cbind(dat[,1:9], log((1+dat[,10:35])))
#standardize to get effect sizes
dat <- cbind(dat[,1:9], decostand(dat[,10:35], method  = "standardize"))

dat2 <- dat[,15:34]
data_long <- dat2 %>%                         
  pivot_longer(colnames(dat2)) %>% 
  as.data.frame()

ggplot(data_long, aes(x = value)) +    
  geom_histogram() + 
  facet_wrap(~ name, scales = "free")

################################################################################
#everthing for SI
dat3 <- as.matrix(dat[,11:34])
corr <- rcorr(dat3, type = "spearman")
p <- round(corr$r, digits = 2)
ggcorrplot(p,
           type = "lower",
           outline.color = "white",
           colors = c("mediumpurple1", "white", "orange"))