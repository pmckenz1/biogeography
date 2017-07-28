for (w in c(100,500,1000,5000,10000,50000)) {

# sampling
first.branch.starting <- as.matrix(read.csv(startingfiles[1]))
Ps_copy <- first.branch.starting

# make landmass of zero values
counter <- (matrix(nrow = nrow(first.branch.starting),
                  ncol = ncol(first.branch.starting),
                  rep(0,nrow(first.branch.starting)*ncol(first.branch.starting))) +
              first.branch.starting) * 0

sample_limit <- round(min(c(nrow(first.branch.starting),ncol(first.branch.starting)))/10)
possible_sizes <- seq(1,sample_limit)[c(T,F)]

for (q in 1:w) {

boxsize <- sample(possible_sizes,1)

center <- c(sample(nrow(first.branch.starting),1),sample(ncol(first.branch.starting),1))

bounds <- c((center + (boxsize-1)/2),(center - (boxsize-1)/2))

rowbounds <- bounds[c(3,1)]
colbounds <- bounds[c(4,2)]

if ((sum(center > (boxsize+1)/2) == 2) && (sum(center < c((nrow(first.branch.starting) - (boxsize+1)/2),
                                         (ncol(first.branch.starting) - (boxsize+1)/2))) == 2)) {


  box <- first.branch.starting[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]]
  box
  counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] <- counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] + 1

  Ps_copy[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] <- (mean(box[!is.na(box)]) +  box*(counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] - 1))/(counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]])


}
}

hist(Ps_copy[(!is.na(Ps_copy))][Ps_copy[(!is.na(Ps_copy))] > 0])

}
