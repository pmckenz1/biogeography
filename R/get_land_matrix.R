get_land_matrix <- function(nrows,ncols, lats, lons) {
  # function requires htmltab package
  lons <- seq(lons[1],lons[2],length.out = ncol)
  lats <- seq(lats[1],lats[2],length.out = nrow)
  lat_lon_mat <- matrix(nrow = nrow,ncol = ncol)
  for (i in 1:ncol) {
    for (q in 1:nrow) {
      newlink <- paste0("http://geonetwork3.fao.org/aglw/climatex.php?xcoord=",lons[i],"&ycoord=",lats[q],"&dddms=dd")
      if (length(suppressMessages(htmltab(newlink,which=1))) > 4) {
        lat_lon_mat[q,i] <- 1  # if the cell coordinate is on land
      } else {lat_lon_mat[q,i] <- 0}  # if the cell coordinate is in water
      print(paste0("col = ",i,", row = ",q))
    }
  }
}