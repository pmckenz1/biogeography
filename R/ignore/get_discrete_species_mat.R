#' Get a matrix of cells with observations of a species
#' 
#' @param datalatitude a vector of latitudes of observations
#' @param datalongitude a vector of longitudes of observations
#' @param nrows the number of rows desired in the output matrix
#' @param ncols the number of columns desired in the output matrix
#' @param borderlats the boundary latitudes of the matrix
#' @param borderlons the boundary longitudes of the matrix
#' 
#' @details This function produces a matrix with discrete cells of value 1 representing points present in the data.
#' @details Vectors of latitude and longitude for observations can be acquired from gbif, and matrix resolution and boundary values should be synced with get_land_matrix().
#' 
get_discrete_species_mat <- function(datalatitude,datalongitude,nrows,ncols,borderlats,borderlons) {
  species_mat <- matrix(data = rep(0,(nrow*ncol)),nrow = nrow,ncol = ncol)
  lons <- seq(borderlons[1],borderlons[2],length.out = ncol)
  lats <- seq(borderlats[1],borderlats[2],length.out = nrow)
  for (i in 1:length(datalatitude)) {
    closeness_lon <- abs(lons - datalongitude[i])
    closeness_lat <- abs(lats - datalatitude[i])
    species_mat[(closeness_lat == min(closeness_lat)),(closeness_lon == min(closeness_lon))] <- 1
  }
  species_mat
}
