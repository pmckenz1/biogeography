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



#time testing
ptm <- proc.time()
enviro_matrix <- read.csv("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0enviro_matrix.csv")
enviro_matrix <- as.matrix(enviro_matrix)
first.branch.starting <- read.csv(startingfiles[1])
first.branch.starting <- as.matrix(first.branch.starting)
next_step <- get_prob_matrix(Ps = first.branch.starting, Es = enviro_matrix, alpha = .5, beta = .8)
proc.time() - ptm

#first run, no progress bar: 1097.350
#second run, with progress bar: 1165.247
#third run, with progress bar and as matrices: 111.005

# compute tips to first interior node
enviro_matrix <- as.matrix(read.csv("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0enviro_matrix.csv"))
for (p in 2:3) {
  first.branch.starting <- as.matrix(read.csv(startingfiles[p]))
  for (i in 1:round(edgelength.tips[p])) {
    first.branch.starting <- get_prob_matrix(Ps = first.branch.starting, Es = enviro_matrix, alpha = .5, beta = .8)
    print(paste0("done with ",i," in ", round(edgelength.tips[p])))
  }
  write.csv(first.branch.starting,
            paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/",paste0(edge.from.to.tips[p,], collapse = "_"),".csv"),
            row.names = FALSE)
}

test1 <- as.matrix(read.csv("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/45_34.csv"))
test2 <- as.matrix(read.csv("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/47_40.csv"))

test <- test1*test2/max(test1*test2,na.rm = T)
plot(as.raster(test1))
