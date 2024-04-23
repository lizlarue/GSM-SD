###############################################################################
#Combine site files into one dataset
################################################################################
setwd("./data/")
library(stringr)
library(tidyr)
library(sf)
#sessionInfo()

######################## combine 28 site files ##############################################
#first check to see if there are overlapping areas
sites <- read.csv("./sites/Neon28sites.csv") #epsg code for UTM WGS84 zone
randpoints <- list.files("./extracted_11172023", full.names = TRUE)
#################################
OUT <- NULL
for(i in 1:length(randpoints)){
  
  site.i <- read.csv(randpoints[i])  
  site.i <- site.i[,-1]
  site.i$index <- row.names(site.i)
  site.name <- str_sub(randpoints[i], start=22, end = 25)

  points <- st_as_sf(site.i, coords = c("easting", "northing"), crs = sites[sites$site == site.name,3])
  
  # Create a buffer around the points
  points <- st_buffer(points, endCapStyle="SQUARE", dist = 60)
  #plot(points[1])
  overlaps <- st_overlaps(points)
  
  list_length <- dim(overlaps)[1]
  X <- NULL
  for(f in 1:list_length){
    x <- overlaps[[f]]
    if (length(x) == 0) {next}
    X <- c(X, x)
  }
  site.i <- site.i[!(site.i$index %in% X),]
  points <- points[!(points$index %in% X),]
  plot(points[1])
  
  OUT <- rbind(OUT, site.i)
}

OUT <- OUT[!OUT$chm < 3, ]
table(OUT$site)

#take first 50
OUT2 <- NULL
for(n in 1:nrow(sites)){
  site.name <- sites[n, "site"]

  OUT.n <- OUT[OUT$site == site.name, ]
  OUT.n <- OUT.n[1:50,]
  OUT2 <- rbind(OUT2, OUT.n)
}

table(OUT2$site)

write.csv(OUT2, "./aggregated/NEON28GSM_11172023.csv")