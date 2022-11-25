## commented so it doesn't get included in the script itself!!
## library(knitr)
## knitr::purl(input = "ploidy_project_loop.Rmd",
##             output = paste("ploidy_loop_",
##                            format(Sys.time(),
##                                   "%Y-%m-%e_%H-%M"), ".R",
##                            sep = ""),
##             documentation = 0)

##install.packages(pkgs = "~/hacking/dismo/", repos = NULL)
#library(devtools)
#load_all(path = "~/hacking/dismo/")
#install_github("danlwarren/ENMTools")
#install.packages(pkgs = "~/hacking/dismo", repos = NULL)

library(ecospat)
library(raster)
library(maptools)
library(dismo)
library(rgeos)
library(ENMTools)
library(ade4)
library(geodata)
library(geosphere)

devtools::session_info()

comb <- read.csv("data/Hybrid_parent_comb.csv") #csv with
  #hybrid and its parents 
hybrids <- comb[,1]
hybrids <- gsub(" ", "", hybrids) #remove extra spaces
parent1s <- comb[,2]
parent1s <- gsub(" ", "", parent1s)
parent2s <- comb[,3]
parent2s <- gsub(" ", "", parent2s)
parent2s[parent2s == ""] <- NA

#set date string
date.str <- Sys.Date()

if(!dir.exists("./data/world")) {
  dir.create("./data/world")
}

if(!dir.exists("./Plots")) {
  dir.create("./Plots")
}


## get worldclim data 

wc <- getData(name = "worldclim", var = "bio", res = 10,
              path = "data")  

## 111 km per degree (approximate)
res_km <- res(wc) * 111 

w2 <- world(path = "./data/world", resolution = 2)
w2sp <- as(w2, "Spatial")

geoGridR <- 200
geoGrid <- expand.grid(longitude =
                         seq(-180, 180, length.out = geoGridR),
                       latitude =
                         seq(-90, 90, length.out = geoGridR))



#function to get the coordinates from the csv files
getcoords <- function(csv){
  coords <- csv[csv$decimalLongitude !=0 & csv$decimalLatitude !=0, ] #remove zero coordinates
  coords <- subset(coords, !is.na(decimalLongitude)) #remove NAs
  coords <- coords[, c("decimalLongitude", "decimalLatitude")] #just the coordinates
}

#function to extract the climatic conditions where the species occur
envextract <- function(sp.coords){
  env <- extract(wc, sp.coords)
  env <- env[complete.cases(env), ]
}

#function to calculate min, max and temperature breadth for where the species occurs
tempbreadth <- function(env){
  bio5 <- env[, "bio5"]
  maxtemp <- max(bio5, na.rm=T)
  
  bio6 <- env[, "bio6"]
  mintemp <- min(bio6, na.rm=T)
  
  tempbreadth <- maxtemp - mintemp
  
  results <- c(maxtemp, mintemp, tempbreadth)
  return(results)
}

precipbreadth <- function(env){
  Bio13 <- env[, "bio13"]
  maxprec <- max(Bio13, na.rm=T)
  
  Bio14 <- env[, "bio14"]
  minprec <- min(Bio14, na.rm=T)
  
  precbreadth <- env[, "bio15"]
  precbreadth <- max(precbreadth, na.rm = T)
  
  results <- c(maxprec, minprec, precbreadth)
}

env.range.area <- function(pca.scores){
  ## Create a raster for the PCA space ##
  ## create an empty raster:

  ## extent is based on the scores for the species, but the scale is set by
  ## the global PCA.
  ## Starts as a 10 x 10 grid:
  mask.rast <- raster(extent(range(pca.scores[,1]), range(pca.scores[,2])))

  ## snap "out" to make sure extend always increases the extent!
  mask.rast <- extend(mask.rast, extent(mask.rast) + 1, snap = "out")
  res(mask.rast) <- 0.2 ## resolution is a fixed value for all species
  
  #set the background cells in the raster to 0
  mask.rast[!is.na(mask.rast)] <- 0 
  
  #set the cells that contain points to 1 
  env_range <- rasterize(pca.scores, mask.rast, field = 1)
  env_range <- merge(env_range, mask.rast)
  
  #determine resolution of the raster
  res_rast <- res(env_range)
  
  #determine number of cells that have value of 1
  ncells <- freq(env_range, value= 1, useNA= "no")
  
  # calculate area. Since scores are from the global PCA, the total
  # environmental range for every species should be on the same scale (I
  # think) 
  area <- res_rast[1] * res_rast[2] * ncells
  return(area)
}


## Global PCA
glob <- getValues(wc)
glob <- glob[complete.cases(glob), ]

pca.env <- dudi.pca(glob, scannf=F, nf = 2)
  
## PCA scores for the whole study area
scores.globclim <- pca.env$li

Summary.File <- NULL

recNums <- numeric()
collections <- character()

for (i in 1:nrow(comb)){
##for (i in 108:109){
  #i = 108 ## for testing
  message("**** I = ", i, " of ", nrow(comb), " ****")
  message("     ", date())
  ## get the names of the species
  hybrid_name <- hybrids[i]
  parent1_name <- parent1s[i]
  if (!is.na(parent2s[i])) {
    parent2_name <- parent2s[i]
  } else {
    parent2_name <- NA
  }

  message("     Processing: ", hybrid_name)
  
  ## read in data for the species and extract coordinates 
  hybrid <- read.csv(paste("data/records/", hybrids[i], ".csv", sep="")) 
  parent1 <- read.csv(paste("data/records/", parent1s[i], ".csv", sep=""))

  hybridC <- getcoords(hybrid)

  ## Exclude strange high-arctic record that is clearly wrong:
  if(hybrid_name == "Leucaena_leucocephala"){
    hybridC <- hybridC[hybridC[, 2] < 89, ]
  }

  hybridC <- SpatialPoints(hybridC, proj4string=CRS(proj4string(w2sp)))
  hybridC <- hybridC[w2sp] #remove points not on land
  hybridC <- data.frame(coordinates(hybridC)) #extract coordinates and convert to dataframe

  parent1C <- getcoords(parent1)
  parent1C <- SpatialPoints(parent1C, proj4string=CRS(proj4string(w2sp)))
  parent1C <- parent1C[w2sp] #remove points not on land
  parent1C <- data.frame(coordinates(parent1C))

  ## Record observation numbers:

  recNumsNew <- c(nrow(hybridC), nrow(parent1C))
  names(recNumsNew) <- c(hybrid_name, parent1_name)
  recNums <- c(recNums, recNumsNew)

  collections <- c(collections, hybrid$institutionCode,
                   parent1$institutionCode) 

  #Determine max and min latitude 
  maxlat_hybrid <- max(abs(hybridC$decimalLatitude))
  maxlat_parent1 <- max(abs(parent1C$decimalLatitude))
  
  minlat_hybrid <- min(abs(hybridC$decimalLatitude))
  minlat_parent1 <- min(abs(parent1C$decimalLatitude))

  hybrid_lat_range <- abs(diff(range(hybridC$decimalLatitude)))
  parent1_lat_range <- abs(diff(range(parent1$decimalLatitude)))
  
  #Determine the centroid
  hybrid_centroid <- centroid(hybridC)
  centroid_latitude <- hybrid_centroid[2]
  centroid_longitude <- hybrid_centroid[1]
  parent1_centroid <- centroid(parent1C)
  P1_centroid_latitude <- parent1_centroid[2]
  P1_centroid_longitude <- parent1_centroid[1]
  
  ## Extract climate conditions where the species occurs 

  message("      extracting climate data")
  hybrid_env <- envextract(hybridC)
  parent1_env <- envextract(parent1C)

  ## determining max,min temp and precipitation breadths for each species and parent 1 
  message("      temp")
  hybrid_temptest <- tempbreadth(hybrid_env)
  hybrid_tempmax <- hybrid_temptest[1]
  hybrid_tempmin <-hybrid_temptest[2]
  hybrid_tempbreadth <- hybrid_temptest[3]
  
  parent1_temptest <- tempbreadth(parent1_env)
  parent1_tempmax <- parent1_temptest[1]
  parent1_tempmin <- parent1_temptest[2]
  parent1_tempbreadth <- parent1_temptest[3]

  message("      precip")
  hybrid_prectest <- precipbreadth(hybrid_env)
  hybrid_precmax <- hybrid_prectest[1]
  hybrid_precmin <-hybrid_prectest[2]
  hybrid_precbreadth <- hybrid_prectest[3]
  
  parent1_prectest <- precipbreadth(parent1_env)
  parent1_precmax <- parent1_prectest[1]
  parent1_precmin <- parent1_prectest[2]
  parent1_precbreadth <- parent1_prectest[3]

  message("      geographic ranges")

  ## Geographical ranges
  ## Don't need to create the same raster twice, so combine this with the
  ## next section:
  ## area_hybridOld <- geo.range.area(hybridC)
  ## area_parent1Old <- geo.range.area(parent1C)
  
  message("      geographic overlap")
  bufferWidth <- 10000 # in meters, 10000 = 10km
  ## Geographical range overlap
  hybrid_enm <- enmtools.species(presence.points = hybridC)
  hybrid_enm$range <-
    background.buffer(points = hybrid_enm$presence.points,
                      buffer.width = bufferWidth,
                      buffer.type = "circles",
                      mask = wc,
                      return.type = "raster") 

  ncells <- freq(hybrid_enm$range, value = 1, useNA= "no")
  area_hybrid <- res_km[1] * res_km[2] * ncells 
  ncells <- NULL ## erase before we reuse, just in case?

  parent1_enm <- enmtools.species(presence.points = parent1C)
  parent1_enm$range <-
    background.buffer(points = parent1_enm$presence.points,
                      buffer.width = bufferWidth,
                      mask = wc,
                      return.type = "raster",
                      buffer.type = "circles") 
  ncells <- freq(parent1_enm$range, value = 1, useNA= "no")
  area_parent1 <- res_km[1] * res_km[2] * ncells 

  comp1_overlap_geo <- geog.range.overlap(hybrid_enm, parent1_enm)
  
  message("     Ecospat")

  message("     ... PCA")
  ## PCA scores for parent1
  scores.parent1 <- suprow(pca.env, parent1_env)$li
  ## PCA scores for hybrid
  scores.hybrid <- suprow(pca.env, hybrid_env)$li

  message("     ... env.range.area")
  ## Environmental ranges
  env_range_area <- env.range.area(scores.hybrid)
  env_range_area_P1 <- env.range.area(scores.parent1)
  
  message("     ... grid.clim.dyn")
  ## Environmental overlap
  grid.clim.parent1 <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                             glob1 = scores.globclim,
                                             sp = scores.parent1, R=100,
                                             th.sp=0)
  
  grid.clim.hybrid <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                            glob1 = scores.globclim,
                                            sp = scores.hybrid, R=100,
                                            th.sp=0)


  message("     ... ecospat.niche.overlap")  
  D.overlap <- ecospat.niche.overlap(grid.clim.parent1, grid.clim.hybrid,
                                     cor = TRUE)
  Schoeners_metric <- D.overlap$D
  Hellinger_metric <- D.overlap$I
  
  jpeg(filename = paste("Plots/niche_", hybrid_name, "+",
                        parent1_name,".jpeg", sep= "")) 
  
  message("     ... ecospat.plot.niche.dyn")  
  ecospat.plot.niche.dyn(grid.clim.parent1, grid.clim.hybrid, quant = 0.25,
                         interest = 2, title = "Niche Overlap",
                         name.axis1 = "PC1", name.axis2 = "PC2")
  dev.off() #prevents the display of the plots 

  message("     ... ecospat.niche.dyn.index")    
  ## determining niche overlap b/w parent 1 and hybrid
  index_parent1comp <-
    ecospat.niche.dyn.index(grid.clim.parent1,
                            grid.clim.hybrid) 
  comp1_overlap <-
    index_parent1comp$dynamic.index.w["stability"] 
  comp1_hybrid_only <-
    index_parent1comp$dynamic.index.w["expansion"] 
  parent1_only <-
    index_parent1comp$dynamic.index.w["unfilling"] 

  message("     Geographic Ecospat")  
  ## Ecospat geographic grids. Need to increase the grid
  ## size to bring resolution below 1 degree

  grid.geo.hybrid <- ecospat.grid.clim.dyn(glob = geoGrid,
                                           glob1 = geoGrid,
                                           sp = hybridC,
                                           R = geoGridR,
                                           th.sp=0, geomask = w2sp)

  grid.geo.parent1 <- ecospat.grid.clim.dyn(glob = geoGrid,
                                            glob1 = geoGrid,
                                            sp = parent1C,
                                            R = geoGridR,
                                            th.sp = 0, geomask = w2sp) 
  
  D.geo.overlap <- ecospat.niche.overlap(grid.geo.parent1, grid.geo.hybrid,
                                         cor = TRUE)
  Schoeners_metricGeo <- D.geo.overlap$D
  Hellinger_metricGeo <- D.geo.overlap$I
  
  jpeg(filename = paste("Plots/geog_", hybrid_name, "+_",
                        parent1_name,".jpeg", sep= "")) 
  
  ecospat.plot.niche.dyn(grid.geo.parent1, grid.geo.hybrid, quant = 0,
                         interest = 2, title = "Niche Overlap",
                         name.axis1 = "PC1", name.axis2 = "PC2")
  dev.off() #prevents the display of the plots 
  
  ## determining niche overlap b/w parent 1 and hybrid
  index_parent1geo <-
    ecospat.niche.dyn.index(grid.geo.parent1, grid.geo.hybrid) 
  comp1_overlapGeo <- index_parent1geo$dynamic.index.w["stability"]
  comp1_hybrid_onlyGeo <- index_parent1geo$dynamic.index.w["expansion"]
  parent1_onlyGeo <- index_parent1geo$dynamic.index.w["unfilling"]


  ## If the hybrid has a parent 2 then redo the above calculations with parent 2:
  if (!is.na(parent2s[i])) {
    message("     Parent 2")  
    parent2_name <- parent2s[i]
    parent2 <- read.csv(paste("data/records/", parent2s[i],
                              ".csv", sep="")) 
    parent2C <- getcoords(parent2)
    parent2C <- SpatialPoints(parent2C, proj4string = CRS(proj4string(w2sp)))
    parent2C <- parent2C[w2sp] #remove points not on land
    parent2C <- data.frame(coordinates(parent2C))

    recNumsNew <- nrow(parent2C)
    names(recNumsNew) <- parent2_name
    recNums <- c(recNums, recNumsNew)

    collections <- c(collections, parent2$institutionCode)
    
    #Determine max and min latitude 
    maxlat_parent2 <- max(abs(parent2C$decimalLatitude))
    minlat_parent2 <- min(abs(parent2C$decimalLatitude))
    parent2_lat_range <- abs(diff(range(parent2$decimalLatitude)))
    
    #determine centroid for parent 2 and position for species with two parents
    parent2_centroid <- centroid(parent2C)
    P2_centroid_latitude <- parent2_centroid[2]
    P2_centroid_longitude <- parent2_centroid[1]

    
    message("     ... extracting env")
    ## Extract environmental occurrence conditions
    parent2_env <- envextract(parent2C)

    ## max, min and breadths for temperature and precipitation
    parent2_temptest <- tempbreadth(parent2_env)
    parent2_tempmax <- parent2_temptest[1]
    parent2_tempmin <- parent2_temptest[2]
    parent2_tempbreadth <- parent2_temptest[3]
    
    parent2_prectest <- precipbreadth(parent2_env)
    parent2_precmax <- parent2_prectest[1]
    parent2_precmin <- parent2_prectest[2]
    parent2_precbreadth <- parent2_prectest[3]

    message("     ... range size & overlap")    
    ## geographical range size and overlap

    parent2_lat <- abs(diff(range(parent2$decimalLatitude)))
    ## area_parent2 <- geo.range.area(parent2C)
    ncells <- NULL
    parent2_enm <- enmtools.species(presence.points = parent2C)
    parent2_enm$range <-
      background.buffer(points = parent2_enm$presence.points,
                        buffer.width = bufferWidth,
                        mask = wc,
                        buffer.type = "circles",
                        return.type = "raster") 

    ncells <- freq(parent2_enm$range, value = 1, useNA= "no")
    area_parent2 <- res_km[1] *res_km[2] * ncells 

    comp2_overlap_geo <- geog.range.overlap(hybrid_enm, parent2_enm)
    
    ## PCA scores for parent2
    scores.parent2 <- suprow(pca.env, parent2_env)$li
    
    ## Environmental range
    env_range_area_P2 <- env.range.area(scores.parent2)

    message("     ... Ecospat")        
    ## Environmental range overlap 
    grid.clim.parent2 <-
      ecospat.grid.clim.dyn(glob = scores.globclim,
                            glob1 = scores.globclim,
                            sp = scores.parent2, R = 100,
                            th.sp=0) 
    
    D.overlap.2 <- ecospat.niche.overlap(grid.clim.parent2,
                                         grid.clim.hybrid, 
                                         cor = TRUE)
    Schoeners_metric_2 <- D.overlap.2$D
    Hellinger_metric_2 <- D.overlap.2$I
    
    jpeg(filename = paste("Plots/niche_", hybrid_name, "+_",
                          parent2_name,".jpeg", sep= "")) 
    
    ecospat.plot.niche.dyn(grid.clim.parent2,
                           grid.clim.hybrid,
                           quant = 0.25, interest = 2,
                           title = "Niche Overlap",
                           name.axis1 = "PC1", name.axis2 = "PC2")
    dev.off()
    
    ## determining niche overlap b/w parent 2 and hybrid
    index_parent2comp <-
      ecospat.niche.dyn.index(grid.clim.parent2,
                              grid.clim.hybrid) 
    comp2_overlap <- index_parent2comp$dynamic.index.w["stability"]
    comp2_hybrid_only <- index_parent2comp$dynamic.index.w["expansion"]
    parent2_only <- index_parent2comp$dynamic.index.w["unfilling"]

    ###between parents
    index_parentscomp <- ecospat.niche.dyn.index(grid.clim.parent1, grid.clim.parent2)
    comp3_overlap <- index_parentscomp$dynamic.index.w["stability"]
    comp3_parent2_only <- index_parentscomp$dynamic.index.w["expansion"]
    comp3_parent1_only <- index_parentscomp$dynamic.index.w["unfilling"]

    message("     ... Geographic ecospat")        
    ## Geographic grids
    grid.geo.parent2 <-
      ecospat.grid.clim.dyn(glob = geoGrid, glob1 = geoGrid,
                            sp = parent2C, R = geoGridR,
                            th.sp=0, geomask = w2sp)
    
    D.geo.overlap.2 <- ecospat.niche.overlap(grid.geo.parent2,
                                             grid.geo.hybrid,
                                             cor = TRUE)
    Schoeners_metricGeo_2 <- D.geo.overlap.2$D
    Hellinger_metricGeo_2 <- D.geo.overlap.2$I
    
    jpeg(filename = paste("Plots/geog_", hybrid_name, "+_",
                          parent2_name,".jpeg", sep= "")) 
    
    ecospat.plot.niche.dyn(grid.geo.parent2, grid.geo.hybrid, quant = 0,
                           interest = 2, title = "Niche Overlap",
                           name.axis1 = "PC1", name.axis2 = "PC2")
    dev.off() #prevents the display of the plots 
    
    ## determining niche overlap b/w parent 1 and hybrid
    index_parent2geo <-
      ecospat.niche.dyn.index(grid.geo.parent2,
                              grid.geo.hybrid)  
    comp2_overlapGeo <- index_parent2geo$dynamic.index.w["stability"]
    comp2_hybrid_onlyGeo <- index_parent2geo$dynamic.index.w["expansion"]
    parent2_onlyGeo <- index_parent2geo$dynamic.index.w["unfilling"]

    ###between parents 
    index_parentsgeo <- ecospat.niche.dyn.index(grid.geo.parent1, grid.geo.parent2) 
    comp3_overlapGeo <- index_parentsgeo$dynamic.index.w["stability"]
    comp3_parent2_onlyGeo <- index_parentsgeo$dynamic.index.w["expansion"]
    comp3_parent1_onlyGeo <- index_parentsgeo$dynamic.index.w["unfilling"]

  }
  
  ## If the hybrid does not have a parent 2 then define all of the variables as "NA"
  if (is.na(parent2s[i])){
    parent2_lat <- NA
    maxlat_parent2 <- NA
    minlat_parent2 <- NA
    parent2_tempmax <- NA
    parent2_tempmin <- NA
    parent2_tempbreadth <- NA
    parent2_precmax <- NA
    parent2_precmin <- NA
    parent2_precbreadth <- NA
    area_parent2 <- NA
    comp2_overlap_geo <- NA
    env_range_area_P2 <- NA
    Schoeners_metric_2 <- NA
    Hellinger_metric_2 <- NA
    comp2_overlap <- NA
    comp2_hybrid_only <- NA
    parent2_only <- NA
    Schoeners_metricGeo_2 <- NA
    Hellinger_metricGeo_2 <- NA
    comp2_overlapGeo <- NA
    parent1_onlyGeo <- NA
    comp2_hybrid_onlyGeo <- NA
    parent2_onlyGeo <- NA
    comp3_overlap <- NA
    comp3_parent2_only <- NA
    comp3_parent1_only <- NA
    comp3_overlapGeo <- NA
    comp3_parent2_onlyGeo <- NA
    comp3_parent1_onlyGeo <- NA
    P2_centroid_longitude <- NA
    P2_centroid_latitude <- NA
  }
  
  data.line <-
    data.frame(allo = hybrid_name, p1 = parent1_name, p2 = parent2_name, 

               ## Latitudinal range
               lat_range_hy = hybrid_lat_range,
               lat_range_p1 = parent1_lat_range,
               lat_range_p2 = parent2_lat_range,

               max_lat_hy = maxlat_hybrid,
               max_lat_p1 = maxlat_parent1,
               max_lat_p2 = maxlat_parent2,
               min_lat_hy = minlat_hybrid,
               min_lat_p1 = minlat_parent1,
               min_lat_p2 = minlat_parent2,

               temp_max_hy = hybrid_tempmax,
               temp_max_p1 = parent1_tempmax,
               temp_max_p2 = parent2_tempmax, 
               temp_min_hy = hybrid_tempmin,
               temp_min_p1 = parent1_tempmin,
               temp_min_p2 = parent2_tempmin,
               temp_breadth_hy = hybrid_tempbreadth,
               temp_breadth_p1 = parent1_tempbreadth, 
               temp_breadth_p2 = parent2_tempbreadth,

               precip_max_hy = hybrid_precmax,
               precip_max_p1 = parent1_precmax,
               precip_max_p2 = parent2_precmax,
               precip_min_hy = hybrid_precmin,
               precip_min_p1 = parent1_precmin,
               precip_min_p2 = parent2_precmin,
               precip_breadth_hy = hybrid_precbreadth,
               precip_breadth_p1 = parent1_precbreadth, 
               precip_breadth_p2 = parent2_precbreadth,

               range_km2_hy = area_hybrid,
               range_km2_p1 = area_parent1,
               range_km2_p2 = area_parent2,
               
               geo_stability_p1h = comp1_overlap_geo,
               geo_stability_p2h = comp2_overlap_geo, 

               env_range_hy = env_range_area,
               env_range_p1 = env_range_area_P1,
               env_range_p2 = env_range_area_P2,

               schoeners_overlap_metric_p1h = Schoeners_metric,
               schoeners_overlap_metric_p2h = Schoeners_metric_2,

               hellinger_metric_p1h = Hellinger_metric,
               hellinger_metric_p2h = Hellinger_metric_2,

               env_stability_p1h = comp1_overlap, 
               env_stability_p2h = comp2_overlap,

               env_expansion_p1h = comp1_hybrid_only,
               env_unfilling_p1h = parent1_only, 

               env_expansion_p2h = comp2_hybrid_only,
               env_unfilling_p2h = parent2_only,

               schoeners_geo_p1h = Schoeners_metricGeo,
               schoeners_geo_p2h = Schoeners_metricGeo_2,

               hellinger_geo_p1h = Hellinger_metricGeo,
               hellinger_geo_p2h = Hellinger_metricGeo_2,

               geo_stability_p1h = comp1_overlapGeo,
               geo_stablity_p2h = comp2_overlapGeo,

               geo_expansion_p1h = comp1_hybrid_onlyGeo,
               geo_unfilling_p1h = parent1_onlyGeo,

               geo_expansion_p2h = comp2_hybrid_onlyGeo,
               geo_unfilling_p2h = parent2_onlyGeo,

               env_stability_p1p2 = comp3_overlap,
               env_expansion_p1p2 = comp3_parent2_only,
               env_unfilling_p1p2 = comp3_parent1_only,

               geo_stability_p1p2 = comp3_overlapGeo,
               geo_expansion_p1p2 = comp3_parent2_onlyGeo,
               geo_unfilling_p1p2 = comp3_parent1_onlyGeo, 

               centroid_longitude_hy = centroid_longitude,
               centroid_latitude_hy = centroid_latitude,
               centroid_longitude_p1 = P1_centroid_longitude, 
               centroid_latitude_p1 = P1_centroid_latitude, 
               centroid_longitude_p2 = P2_centroid_longitude,
               centroid_latitude_p2 = P2_centroid_latitude)
  
  Summary.File <- rbind(Summary.File, data.line)
}

row.names(Summary.File) <- Summary.File$Species

write.table(Summary.File,
            file = (paste("Ploidy_project_analysis_", date.str, ".csv",
                          sep=""))) 

write.table(recNums, "record_numbers.csv")
write.table(collections, "collections.csv")

