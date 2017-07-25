taxa <- c('Acrocomia_aculeata', 'Acrocomia_crispa',
'Aiphanes_horrida', 'Aiphanes_minima',
'Allagoptera_caudescens', 'Astrocaryum_acaule',
'Astrocaryum_aculeatum', 'Hexopetion_alatum',
'Astrocaryum_campestre', 'Astrocaryum_carnosum',
'Astrocaryum_chambira', 'Astrocaryum_chonta',
'Astrocaryum_faranae', 'Astrocaryum_farinosum',
'Astrocaryum_ferrugineum', 'Astrocaryum_gratum',
'Astrocaryum_gynacanthum', 'Astrocaryum_huaimi',
'Astrocaryum_huicungo', 'Astrocaryum_jauari',
'Astrocaryum_javarense', 'Astrocaryum_macrocalyx',
'Astrocaryum_malybo', 'Hexopetion_mexicanum',
'Astrocaryum_minus', 'Astrocaryum_murumuru',
'Astrocaryum_paramaca', 'Astrocaryum_perangustatum',
'Astrocaryum_rodriguesii', 'Astrocaryum_sciophilum',
'Astrocaryum_scopatum', 'Astrocaryum_sociale',
'Astrocaryum_standleyanum', 'Astrocaryum_ulei',
'Astrocaryum_urostachys', 'Astrocaryum_vulgare',
'Attalea_phalerata', 'Bactris_bifida',
'Bactris_gasipaes', 'Barcella_odora',
'Beccariophoenix_madagascariensis', 'Butia_eriospatha',
'Cocos_nucifera','Desmoncus_orthacanthos',
'Desmoncus_polyacanthos', 'Elaeis_oleifera',
'Jubaea_chilensis', 'Jubaeopsis_caffra',
'Lytocaryum_weddellianum', 'Parajubaea_cocoides',
'Reinhardtia_simplex', 'Voanioala_gerardii')

for (i in 1:length(taxa)) {
  dir.create(paste0("/Users/pmckenz1/Desktop/projects/matrixranger/astrocaryum_proj/",taxa)[i])
}




for (i in taxa[8:length(taxa)]) {
  sdm_function(directory.path = paste0("/Users/pmckenz1/Desktop/projects/matrixranger/astrocaryum_proj/",i),
               coords = c(-130,-25,-45,35),
               species.name = gsub("_"," ",i),
               climate.map.path = "/Users/pmckenz1/Desktop/projects/example_species_dist/climate_map/wc2-5",
               threshold = "10pct",
               raster.location = "/Users/pmckenz1/Desktop/projects/matrixranger/astrocaryum_proj/0environ_rasters",
               do.plot = FALSE)
  print(i)
}


library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-130, -25), ylim = c(-45, 35), asp = 1)

points(acro_occur$lon,acro_occur$lat,pch = 18, col = "red")

test <- conv_sdm_rasters_threshold(paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/",taxa[1],"/outputs/",taxa[1],".asc"),
                           paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/",taxa[1],"/outputs/maxentResults.csv"),
                           threshold = "10pct")
astro_rasters <- raster(paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/",taxa[1],"/outputs/",taxa[1],".asc"))
for( i in 2:length(taxa)) {
  rastername<- paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/",taxa[i],"/outputs/",taxa[i],".asc")
  if(file.exists(rastername)) {
    result<- raster(rastername)
    astro_rasters <- stack(astro_rasters,result)
    print(i)
  }
}
enviro_matrix <- sum(astro_rasters)
enviro_matrix <- as.matrix(enviro_matrix)

enviro_matrix[!is.na(enviro_matrix)] <- (enviro_matrix[!is.na(enviro_matrix)] - min(enviro_matrix[!is.na(enviro_matrix)]))/max(enviro_matrix[!is.na(enviro_matrix)])
image(enviro_matrix)
write.csv(file = "/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0enviro_matrix.csv",enviro_matrix,row.names = FALSE)

for(i in 1:length(taxa)) {
  rastername <- paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/",taxa[i],"/outputs/",taxa[i],".asc")
  if(file.exists(rastername)) {
  startingrange <- conv_sdm_rasters_threshold(rastername,
                                   paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/",taxa[i],"/outputs/maxentResults.csv"),
                                   threshold = "10pct")
  startingrange <- as.matrix(startingrange)
  write.csv(file = paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0startingranges/",taxa[i],".csv"),startingrange,row.names = FALSE)
  }
  print(i)
}

library(ape)
dated.tree <- read.nexus("/Users/pmckenz1/Desktop/projects/matrixranger/data/Astrocaryum_dated_tree_in_BEAST.nex")


species_names <- list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0startingranges")
species_names <- gsub(".csv","",species_names)

pruned.dated.tree <- drop.tip(dated.tree,dated.tree$tip.label[!(dated.tree$tip.label %in% species_names)])

temp <- pruned.dated.tree$edge.length[1]
pruned.dated.tree$edge.length[1]<- temp


 temp <- pruned.dated.tree$edge.length[77]
pruned.dated.tree$edge.length[77] <- temp
plot(pruned.dated.tree)
pruned.dated.tree$tip.label
