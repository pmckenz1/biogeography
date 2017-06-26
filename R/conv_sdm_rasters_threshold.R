#' Convert maxent output rasters to binary by threshold values
#' @param raster_path path to the .asc file
#' @param results_file path to the results file
#' @param threshold either "balanced", "10pct", or "minimum"

conv_sdm_rasters_threshold <- function(raster_path,results_file,threshold) {
  maxent_output <-raster(raster_path)
  results <- read.csv(results_file)
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
  rc
}
