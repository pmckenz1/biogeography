#' Implement box movement method
#' @param Ps Occurrence matrix. Values between (inclusive) 0 and 1.
#' @param Es The environmental matrix. Same dimensions as Ps, values between (inclusive) 0 and 1.
#' @param branch_length from phylogeny
#' @param box_samples the number of boxes to sample from within the range of the Ps/Es matrices

box_movement_method <- function(Ps, Es, branch_length, box_samples = 10000) {
  # Save a copy of the occurrence matrix to use as reinforcement later on
  Ps_copy <- Ps
  for (i in 1:branch_length) {
    # make landmass of zero values. This will keep track of each time a cell
    # has been selected, for denominator when averaging at the end
    counter <- (matrix(nrow = nrow(Ps),
                       ncol = ncol(Ps),
                       rep(0,nrow(Ps)*ncol(Ps))) +
                  Ps) * 0
    # This is landmass of zeroes as well, to be filled with averages from the random boxes.
    empty_filler <- (matrix(nrow = nrow(Ps),
                            ncol = ncol(Ps),
                            rep(0,nrow(Ps)*ncol(Ps))) +
                       Ps) * 0

    # Set limits on our sampling scheme. Limits defined by matrix size.
    sample_limit <- round(min(c(nrow(Ps),ncol(Ps)))/10)
    possible_sizes <- seq(1,sample_limit)[c(T,F)]

    for (q in 1:(box_samples)) {
      # sample size of our box
      boxsize <- sample(possible_sizes,1)
      # sample center for box
      center <- c(sample(nrow(Ps),1),sample(ncol(Ps),1))
      # box bounds for indexing
      bounds <- c((center + (boxsize-1)/2),(center - (boxsize-1)/2))
      rowbounds <- bounds[c(3,1)]
      colbounds <- bounds[c(4,2)]

      # "if our box falls completely within the matrix"
      if ((sum(center > (boxsize+1)/2) == 2) && (sum(center < c((nrow(Ps) - (boxsize+1)/2),
                                                                (ncol(Ps) - (boxsize+1)/2))) == 2)) {
        # fill the box with starting occurrence data
        box <- Ps[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]]

        # add one to the counter for each cell within the box
        counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] <- counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] + 1

        # calculate the mean value across all non-NA cells in the box
        boxmean <- mean(box,na.rm = T)

        # "if the box isn't in the middle of the ocean"
        if (!is.na(boxmean)) {
          # replace values within the box with a new adjusted mean (this will very by cell as different cells will have different counter values)
          box <- (mean(box[!is.na(box)]) + box*(counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] - 1))/counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]]
          # adjust these values to favor primarily starting occurrence data, but also environmental data
          box <- box + (box*mean(Es[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]],na.rm = T))^2 + (box*mean(Ps_copy[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]],na.rm = T))
          # finally add the values from the box to our empty (but filling up!) matrix
          empty_filler[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] <- box
        }
      }
    }
    # I included this initially to keep a good distribution of values, but I'm realizing
    # that the sqrt part may need to be relocated inside the loop (only happen to box values)
    empty_filler <- sqrt(empty_filler/max(empty_filler,na.rm = TRUE))

    # You can print/plot stuff to keep track of values through the run if you want
    #plot(as.raster(empty_filler))
    #vals <- empty_filler[!is.na(empty_filler)]
    #sum(vals)
    #vals <- vals[vals > 0]
    #hist(vals)

    # Finally, replace the starting occurrence data and run it again.
    Ps <- empty_filler
    print(paste0("Done with ",i," in ", round(branch_length)))
  }
  empty_filler
}
