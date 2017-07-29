
# sampling
first.branch.starting <- as.matrix(read.csv(startingfiles[1]))

box_movement_method <- function(Ps, Es, branch_length,sample_multiplier) {
  Ps_copy <- Ps
  for (i in 1:branch_length) {
# make landmass of zero values
counter <- (matrix(nrow = nrow(Ps),
                  ncol = ncol(Ps),
                  rep(0,nrow(Ps)*ncol(Ps))) +
              Ps) * 0
empty_filler <- (matrix(nrow = nrow(Ps),
                        ncol = ncol(Ps),
                        rep(0,nrow(Ps)*ncol(Ps))) +
                   Ps) * 0

sample_limit <- round(min(c(nrow(Ps),ncol(Ps)))/10)
possible_sizes <- seq(1,sample_limit)[c(T,F)]

for (q in 1:(sample_multiplier)) {

boxsize <- sample(possible_sizes,1)

center <- c(sample(nrow(Ps),1),sample(ncol(Ps),1))

bounds <- c((center + (boxsize-1)/2),(center - (boxsize-1)/2))

rowbounds <- bounds[c(3,1)]
colbounds <- bounds[c(4,2)]

if ((sum(center > (boxsize+1)/2) == 2) && (sum(center < c((nrow(Ps) - (boxsize+1)/2),
                                         (ncol(Ps) - (boxsize+1)/2))) == 2)) {

  box <- Ps[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]]
  box
  counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] <- counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] + 1
  boxmean <- mean(box,na.rm = T)

  if (!is.na(boxmean)) {
    box <- (mean(box[!is.na(box)]) + box*(counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] - 1))/counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]]
    box <- box + (box*mean(Es[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]],na.rm = T))^2 + (box*mean(Ps_copy[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]],na.rm = T))
    empty_filler[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] <- box
  }
#  Ps_copy[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] <- (mean(box[!is.na(box)]) +  box*(counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] - 1))/(counter[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]])
#  Ps_copy[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] <- Ps_copy[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] + 1*Ps_copy[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]] * enviro_matrix[rowbounds[1]:rowbounds[2],colbounds[1]:colbounds[2]]

}
}
empty_filler <- sqrt(empty_filler/max(empty_filler,na.rm = TRUE))
plot(as.raster(empty_filler))

vals <- empty_filler[!is.na(empty_filler)]
sum(vals)
vals <- vals[vals > 0]
#hist(vals)
Ps <- empty_filler
print(paste0("Done with ",i," in ", round(branch_length)))
}
empty_filler
}


###############################

ntips <- length(pruned.dated.tree$tip.label)

edgelengths <- pruned.dated.tree$edge.length
#standardize edge lengths between 1 and 100
edgelengths <- (edgelengths * (99 / (max(edgelengths) - min(edgelengths))))
edgelengths <- edgelengths - min(edgelengths) + 1

nedges <- length(pruned.dated.tree$edge.length)

edges.to.tips <- pruned.dated.tree$edge[,2] <= ntips

edge.from.to.tips <- pruned.dated.tree$edge[edges.to.tips,]
edgelength.tips<- edgelengths[edges.to.tips]

tips.labels <- (pruned.dated.tree$tip.label[edge.from.to.tips[,2]])

#just for memory, associate each tip value with the species it represents
tip.numbers <- cbind(tips.labels,edge.from.to.tips[,2])

startingfiles <- paste0(tips.labels,".csv")
startingfiles <- paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0startingranges/",startingfiles)

enviro_matrix <- as.matrix(read.csv("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0enviro_matrix.csv"))

for (p in 37:38) {
  first.branch.starting <- as.matrix(read.csv(startingfiles[p]))
  first.branch.starting <- box_movement_method(Ps = first.branch.starting, Es = enviro_matrix, branch_length = edgelength.tips[p], sample_multiplier = 10000)
  write.csv(first.branch.starting,
            paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/",paste0(edge.from.to.tips[p,], collapse = "_"),".csv"),
            row.names = FALSE)
  print(p)
}

test1 <- as.matrix(read.csv("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/79_23.csv"))
test2 <- as.matrix(read.csv("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/79_24.csv"))

plot(as.raster(test1))
plot(as.raster(test2))
plot(as.raster(sqrt(test1*test2)/max(sqrt(test1*test2),na.rm = T)))
