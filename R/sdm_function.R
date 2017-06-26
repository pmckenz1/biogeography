#' Gets map of species distribution from maxent model
#' @param directory.path path to an empty directory to store data and maxent outputs
#' @param coords a vector of four boundary coordinates: min lon, max lon, min lat, max lat
#' @param species.name the scientific name of the species to be mapped
#' @param climate.map.path map to a global climate map. The path should be to the folder containing 19 raster layers. If NA, the map will be downloaded to directory.path.
#' @param threshold Default is "balanced," but other options include "10pct" and "minimum"

library(raster)
library(rgbif)

sdm_function <- function(directory.path,coords = NULL,species.name,climate.map.path = NULL,threshold = "balanced") {
print("Getting climate data...")
if (is.null(climate.map.path)) {
climate_data <- getData("worldclim", var="bio", res=2.5,path = directory.path)
} else {
  rasterlist <- list()
  for (i in 1:19) {
    rasterlist[[i]] <- raster(paste0(climate.map.path,"/bio",i,".bil"))
  }
  climate_data <- stack(rasterlist)
}
names(climate_data) <- c("Annual Mean Temp", "Mean Diurnal Range", "Isothermality", "Temp Seasonality", 
                           "Max Temp Warmest Month", "Min Temp Coldest Month", "Temp Annual Range", "Mean Temp Wettest Quarter",
                           "Mean Temp Driest Quarter", "Mean Temp Warmest Quarter", "Mean Temp Coldest Quarter", "Annual Precip", 
                           "Precip Wettest Month", "Precip Driest Month", "Precip Seasonality", "Precip Wettest Quarter", 
                           "Precip Driest Quarter", "Precip Warmest Quarter", "Precip Coldest Quarter")
print("......Finished")

if (is.null(coords)) {
  coords <- c(-180,180,-90,90)
}

ext<-extent(coords) 
clim_map_US <- crop(climate_data,ext)
print("Getting occurrence data...")
bluebirds <- occ_data(scientificName = species.name, decimalLatitude = paste0(coords[3],", ",coords[4]), decimalLongitude = paste0(coords[1],", ",coords[2]), hasCoordinate = TRUE, limit = 1000)
bluebirdlat <- bluebirds$data$decimalLatitude
bluebirdlon <- bluebirds$data$decimalLongitude
write.csv(cbind(species = bluebirds$data$name,lon = bluebirdlon,lat = bluebirdlat),paste0(directory.path,"/occurrence_data.csv"),row.names = FALSE)
print("......Finished")

layers<-c(1:19)
print("Writing raster inputs...")
layerpaths <- paste0(directory.path,"/bio",layers)
writeRaster(stack(clim_map_US[[layers]]), layerpaths, 
            bylayer=TRUE, format='ascii', overwrite=T)
print("......Finished")

maxent.jar.path <- "maxent.jar"
environmental.layers.path <- directory.path
sample.data.path <- paste0(directory.path,"/occurrence_data.csv")
output.path <- paste0(directory.path,"/outputs")
dir.create(file.path(directory.path, "/outputs"))
print("Starting maxent...")
system(paste0('java -mx512m -jar ',maxent.jar.path,' nowarnings noprefixes -E "" -E ',gsub(" ","_",species.name),' outputformat=logistic outputdirectory=',output.path,' samplesfile=',sample.data.path,' environmentallayers=',environmental.layers.path,' autorun'))
print("......Finished")

print("Reading maxent outputs...")
maxent_output <-raster(paste0(output.path,"/",gsub(" ","_",species.name),".asc"))
results <- read.csv(paste0(output.path,"/maxentResults.csv"))
if (threshold == "balanced") {
mapthreshold <- results$Balance.training.omission..predicted.area.and.threshold.value.Logistic.threshold
} else if (threshold == "10pct") {
  mapthreshold <- results$X10.percentile.training.presence.Logistic.threshold
} else if (threshold == "minimum") {
  mapthreshold <- results$Minimum.training.presence.Logistic.threshold
}
m <- c(0,mapthreshold,0,mapthreshold,1,1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(maxent_output, rclmat)
plot(crop(rc,ext))
print("......Finished")
list(raw_output = maxent_output,threshold_map = rc, maxent_results = results)
}

