get_prob_matrix <- function(Ps, Es, alpha, beta) {
  if (!(is.matrix(Ps) && is.matrix(Es))) {
    stop("Your data should be input as matrices.")
  }
  #add border of 0 values
  Ps <- cbind(rep(0,nrow(Ps)),Ps,rep(0,nrow(Ps)))
  Ps <- rbind(rep(0,ncol(Ps)),Ps,rep(0,ncol(Ps)))

  Es <- cbind(rep(0,nrow(Es)),Es,rep(0,nrow(Es)))
  Es <- rbind(rep(0,ncol(Es)),Es,rep(0,ncol(Es)))

  # multiply environmental matrix
  Ps <- Ps*Es

  #progress bar!
  pb <- txtProgressBar(min = 0, max = nrow(Ps)*ncol(Ps))

  # save a copy that will become the new one (so cell-by-cell changes don't affect calcs)
  Ps_copy <- Ps

  #movement
  for (q in 2:(nrow(Ps)-1)) {
    for (w in 2:(ncol(Ps)-1)) {
      currentval <- Ps[q,w]
      if (!(is.na(currentval))) {
        # calculate ratio of occupied cells in 9*9 square surrounding each cell
        numpossible<- length(na.omit(unlist(Ps[((q-1):(q+1)),((w-1):(w+1))])))
        numpresent <- sum(na.omit(unlist(Ps[((q-1):(q+1)),((w-1):(w+1))])))
        Nbar_ratio <- numpresent/numpossible

        Ps_copy[q,w] <- ( ( (1-currentval) * (Nbar_ratio)) +
                       (currentval * (Nbar_ratio)) )

      }
      #update progress
      if (q*w %% 1000 == 0) {
        setTxtProgressBar(pb,q*w)
      }
    }
  }

  # finish progress bar
  close(pb)
  #remove border
  Ps_copy <- Ps_copy[2:(nrow(Ps_copy)-1),2:(ncol(Ps_copy)-1)]

  # normalize
  #Ps_copy <- Ps_copy/max(Ps,na.rm = TRUE)

  # output
  Ps_copy
}
