library(ape)
genbank <- read.csv("/Volumes/My\ Passport/biogeography_phylogeny/DONEgene_table15.csv",stringsAsFactors = FALSE)

genbankspeciesnames <- unique(genbank$species)

species_genbank <- integer(0)
for (i in 1:length(genbankspeciesnames)) {
  the_species <- genbank[(genbank$species ==  genbankspeciesnames[i]),]
  the_species <- the_species[sample(nrow(the_species),1),]
  species_genbank <- rbind(species_genbank,the_species)
}
species_genbank <- species_genbank[-139:-140,]

sequences <- read.GenBank(species_genbank$id)

# make this a character list
char.seqs <- as.character(sequences)

#adjust the names of the list
for (i in 1:length(names(char.seqs))) {
  id <- names(char.seqs)[i]
  species.name <- species_genbank$species[(species_genbank$id == id)]
  newid <- paste0(species.name,"_",id)
  names(char.seqs)[i] <- newid
}

library(rgbif)
names_gbif <- paste0("Pedicularis ", species_genbank$species)

# get the occurrence data for each species
occurrences <- list()
for (i in 1:nrow(species_genbank)) {
  occurrences[[i]] <- occ_data(scientificName = names_gbif[i], limit = 5000, hasCoordinate = TRUE)
  print(i)
}

# get the number of records found for each species
counts <- integer(0)
for (i in 1:length(occurrences)) {
  counts <- c(counts,occurrences[[i]]$meta$count)
}

# overwrite our sequence list to only include those with at least 10 occurrences
char.seqs <- char.seqs[counts >= 10]

#Now write the file
#write.dna(char.seqs,"/Volumes/My\ Passport/biogeography_phylogeny/sequences.fasta",format = "fasta")

# Now align it and make a phylogeny
#############
############


library(raster)
sdm_function <- function(directory.path,coords = NULL,species.name,climate.map.path = NULL, observation_file_name = "occurrence_data", output_folder = "outputs", threshold = "balanced", asc_present = FALSE, plotoutputs = TRUE) {
  if (!(asc_present)) {
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
    }
    
    if (is.null(coords)) {
      coords <- c(-180,180,-90,90)
  }
  

  print("Getting occurrence data...")
  bluebirds <- occ_data(scientificName = species.name, decimalLatitude = paste0(coords[3],", ",coords[4]), decimalLongitude = paste0(coords[1],", ",coords[2]), hasCoordinate = TRUE, limit = 1000)
  bluebirdlat <- bluebirds$data$decimalLatitude
  bluebirdlon <- bluebirds$data$decimalLongitude
  write.csv(cbind(species = bluebirds$data$name,lon = bluebirdlon,lat = bluebirdlat),paste0(directory.path,"/",observation_file_name,".csv"),row.names = FALSE)
  print("......Finished")
  
  ext<-extent(coords) 
  if (!(asc_present)) {
    clim_map_US <- crop(climate_data,ext)
    layers<-c(1:19)
    print("Writing raster inputs...")
    layerpaths <- paste0(directory.path,"/bio",layers)
    writeRaster(stack(clim_map_US[[layers]]), layerpaths, 
                bylayer=TRUE, format='ascii', overwrite=T)
    print("......Finished")
  }
  

  terminalpath <- gsub("\\ ","\\\\ ",directory.path)
  maxent.jar.path <- "maxent.jar"
  environmental.layers.path <- terminalpath
  sample.data.path <- paste0(terminalpath,"/",observation_file_name,".csv")
  output.path <- paste0(terminalpath,"/",output_folder)
  dir.create(file.path(directory.path, "/",output_folder))
  print("Starting maxent...")
  system(paste0('java -mx512m -jar ',maxent.jar.path,' nowarnings noprefixes -E "" -E ',gsub(" ","_",species.name),' outputformat=logistic outputdirectory=',output.path,' samplesfile=',sample.data.path,' environmentallayers=',environmental.layers.path,' autorun'))
  print("......Finished")
  
  output.path <- paste0(directory.path,"/",output_folder)
  if (plotoutputs == TRUE) {
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
}

underscore_indices <- stringi::stri_locate(names(char.seqs),regex = "_")[,1]
species<- substr(names(char.seqs),1,underscore_indices-1)
species[9] <- "sudetica"
scientificname <- paste0("Pedicularis ",species)

# now run maxent for all 98 species

#for (i in 1:length(species)) {
#  sdm_function(directory.path = "/Volumes/My\ Passport/biogeography_phylogeny/species_dist", 
#             species.name = scientificname[i], 
#             climate.map.path = "/Users/pmckenz1/Desktop/projects/example_species_dist/climate_map/wc2-5",
#             observation_file_name = paste0("obs_",species[i]), 
#             output_folder = paste0("outputs_",species[i]), 
#             threshold = "10pct", 
#             asc_present = TRUE,
#             plotoutputs = FALSE
#             )
#  print(i)
#}

# Now collect outputs

rasterstack <- stack()
for (i in 1:98) {
  indiv.species <- species[i]
  current.directory <- paste0("/Volumes/My\ Passport/biogeography_phylogeny/species_dist/outputs_",indiv.species)
  rasterlayer <- raster(paste0(current.directory,"/Pedicularis_",indiv.species,".asc"))
  results <- read.csv(paste0(current.directory,"/maxentResults.csv"))
  mapthreshold <- results$X10.percentile.training.presence.Logistic.threshold
  m <- c(0,mapthreshold,0,mapthreshold,1,1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  rc <- reclassify(rasterlayer, rclmat)
  rasterstack <- stack(rasterstack,rc)
  print(i)
}
names(rasterstack) <- paste0("P.",species[1:98])

plot(rasterstack)

