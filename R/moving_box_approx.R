
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
#plot(as.raster(empty_filler))

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

library(ape)
dated.tree <- read.nexus("/Users/pmckenz1/Desktop/projects/matrixranger/data/Astrocaryum_dated_tree_in_BEAST.nex")


species_names <- list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0startingranges")
species_names <- gsub(".csv","",species_names)

pruned.dated.tree <- drop.tip(dated.tree,dated.tree$tip.label[!(dated.tree$tip.label %in% species_names)])
pruned.dated.tree$edge.length <- (pruned.dated.tree$edge.length * (99 / (max(pruned.dated.tree$edge.length) - min(pruned.dated.tree$edge.length))))
pruned.dated.tree$edge.length <- pruned.dated.tree$edge.length - min(pruned.dated.tree$edge.length) + 1

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

#running movement up branches
for (p in 1:length(tips.labels)) {
  first.branch.starting <- as.matrix(read.csv(startingfiles[p]))
  first.branch.starting <- box_movement_method(Ps = first.branch.starting, Es = enviro_matrix, branch_length = edgelength.tips[p], sample_multiplier = 10000)
  write.csv(first.branch.starting,
            paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/",paste0(edge.from.to.tips[p,], collapse = "_"),".csv"),
            row.names = FALSE)
  print(p)
}

# saving new nodes!
for (p in 1:length(unique(edge.from.to.tips[,1]))) {
  if (sum(edge.from.to.tips[,1] == unique(edge.from.to.tips[,1])[p]) == 2) {
    internal.node.num <- unique(edge.from.to.tips[,1])[p]
    tip_nums <- edge.from.to.tips[,2][(edge.from.to.tips[,1] == internal.node.num)]
    filenames <- paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/",internal.node.num,"_",tip_nums,".csv")
    first <- as.matrix(read.csv(filenames[1]))
    second <- as.matrix(read.csv(filenames[2]))
    node <- sqrt(first*second)
    write.csv(node,paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes/",internal.node.num,".csv"),row.names = FALSE)
    print(p)
  }
}

for (p in list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes")) {
  plot(as.raster(as.matrix(read.csv(paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes/",p)))))
}

# ONE TIME ONLY get the starting ranges in the nodes folder as well, saved with numbers
for (p in list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0startingranges")) {
  node <- as.matrix(read.csv(paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0startingranges/",p)))
  species <- gsub(".csv","",p)
  nodenumber <- tip.numbers[(tip.numbers[,1] == species),2]
  write.csv(node,
            paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes/",nodenumber,".csv"),
            row.names = FALSE)
}


while(length(list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes")) < (length(pruned.dated.tree$edge.length) + 1)) {

# figure out which new nodes are new and need movement to be run
eligible_branches <- pruned.dated.tree$edge[(pruned.dated.tree$edge[,2] %in% as.numeric(gsub(".csv","",list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes")))),]

torun_branches <- eligible_branches[!(paste0(eligible_branches[,1],"_",eligible_branches[,2]) %in% gsub(".csv","",list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps"))),]

#save the new info all together
if (!is.null(nrow(torun_branches))) {
  torun_branch_info <- cbind(torun_branches,pruned.dated.tree$edge.length[(pruned.dated.tree$edge[,2] %in% torun_branches[,2])])
} else {
  torun_branch_info <- c(torun_branches,pruned.dated.tree$edge.length[(pruned.dated.tree$edge[,2] %in% torun_branches[2])])
  dim(torun_branch_info) <- c(1,3)
  }
#run movement
for (p in 1:nrow(torun_branch_info)) {
  nodenumber <- torun_branch_info[p,2]
  first.branch.starting <- as.matrix(read.csv(paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes/",nodenumber,".csv")))
  first.branch.starting <- box_movement_method(Ps = first.branch.starting, Es = enviro_matrix, branch_length = torun_branch_info[p,3], sample_multiplier = 10000)
  write.csv(first.branch.starting,
            paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/",paste0(torun_branch_info[p,1]),"_",paste0(torun_branch_info[p,2]),".csv"),
            row.names = FALSE)
  print(p)
}

done_branches <- matrix(unlist(strsplit(gsub(".csv","",list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/")),"_")),
       ncol = 2,
       byrow = TRUE)

for (p in 1:length(unique(done_branches[,1]))) {
  if (sum(done_branches[,1] == unique(done_branches[,1])[p]) == 2) {
    prospect <- unique(done_branches[,1])[p]
    if (!(prospect %in% gsub(".csv","",list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes/")))) {
      files <- done_branches[(done_branches[,1] == prospect),]
      filenames <- paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0nodemaps/",prospect,"_",files[,2],".csv")
      first <- as.matrix(read.csv(filenames[1]))
      second <- as.matrix(read.csv(filenames[1]))
      node <- sqrt(first*second)
      write.csv(node,paste0("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes/",prospect,".csv"),row.names = FALSE)
      print(p)
    }
  }
}

}

gsub(".csv","",list.files("/Users/pmckenz1/Desktop/projects/astrocaryum_proj_data/0internal_nodes/"))
