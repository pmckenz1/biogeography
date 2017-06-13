# This is summarized in collecting_inputs_example.pdf

library(rgbif)

lons<- "80,112"
lats <- "18,40"

species <- c("salviiflora","muscoides","trichoglossa","ingens","axillaris","armata","rhinanthoides","rex",
"integrifolia","verticillata","cheilanthifolia","lyrata","lutescens","kansuensis")

for (i in species) {
  dir.create(paste0("example_species_dist/species_dist_maps/",i))
}

dir <- paste0("example_species_dist/species_dist_maps/",species)
speciesnames <- paste0("Pedicularis ",species)
for (i in 1:length(species)) {
  sdm_function(
    directory.path = dir[i],
    coords = c(80,112,18,40),
    species.name = speciesnames[i],
    climate.map.path = "example_species_dist/climate_map/wc2-5"
  )
}

rasterfilenames <- paste0(gsub(" ","_",speciesnames),".asc")
specieslists <- list()
for (i in 1:length(species)) {
  raster_path <- paste0("example_species_dist/species_dist_maps/",species[i],"/outputs/",rasterfilenames[i])
  results_file <- paste0("example_species_dist/species_dist_maps/",species[i],"/outputs/maxentResults.csv")
  specieslists[[i]] <- conv_sdm_rasters_threshold(raster_path = raster_path,
                           results_file = results_file,
                           "balanced")
}
environmental_matrix <- sum(stack(specieslists))
write.csv(as.matrix(environmental_matrix),"example_species_dist/phyloproject_inputs/environmental_matrix.csv",row.names = FALSE)

rasterfilenames <- paste0(gsub(" ","_",speciesnames),".asc")
for (i in 1:length(species)) {
  raster_path <- paste0("example_species_dist/species_dist_maps/",species[i],"/outputs/",rasterfilenames[i])
  results_file <- paste0("example_species_dist/species_dist_maps/",species[i],"/outputs/maxentResults.csv")
  specieslayer <- conv_sdm_rasters_threshold(raster_path = raster_path,
                                                  results_file = results_file,
                                                  "10pct")
  write.csv(as.matrix(specieslayer),paste0("example_species_dist/phyloproject_inputs/",species[i],".csv"),row.names = FALSE)
}


