# correcting to make this give probabilities

get_Nbar <- function(Ps,i,q) {
  if (i > 1 && i < nrow(Ps) && q > 1 && q < ncol(Ps)) {
    Nbar <- sum(Ps[(i-1):(i+1),(q-1):(q+1)]) - Ps[i,q]
  }
  if (i == 1 && q > 1 && q < ncol(Ps)) {
    Nbar <- sum(Ps[(i):(i+1),(q-1):(q+1)]) - Ps[i,q]
  }
  if (i == nrow(Ps) && q > 1 && q < ncol(Ps)) {
    Nbar <- sum(Ps[(i-1):(i),(q-1):(q+1)]) - Ps[i,q]
  }
  if (i > 1 && i < nrow(Ps) && q == 1) {
    Nbar <- sum(Ps[(i-1):(i+1),(q):(q+1)]) - Ps[i,q]
  }
  if (i > 1 && i < nrow(Ps) && q == ncol(Ps)) {
    Nbar <- sum(Ps[(i-1):(i+1),(q-1):(q)]) - Ps[i,q]
  }
  if (i == nrow(Ps) && q == ncol(Ps)) {
    Nbar <- sum(Ps[(i-1):(i),(q-1):(q)]) - Ps[i,q]
  }
  if (i == 1 && q == ncol(Ps)) {
    Nbar <- sum(Ps[(i):(i+1),(q-1):(q)]) - Ps[i,q]
  }
  if (i == nrow(Ps) && q == 1) {
    Nbar <- sum(Ps[(i-1):(i),(q):(q+1)]) - Ps[i,q]
  }
  if (i == 1 && q == 1) {
    Nbar <- sum(Ps[(i):(i+1),(q):(q+1)]) - Ps[i,q]
  }
  Nbar
}

get_prob_matrix <- function(Ps,Es, alpha, beta) {
  probs_older <- Ps
  for (i in 1:nrow(Ps)) { # rows
    for (q in 1:ncol(Ps)) { # columns
      Pcell <- Ps[i,q]
      Ecell <- Es[i,q]
      Nbar <- get_Nbar(Ps,i,q)
      Polder <- Ecell * ( ( (1-Pcell) * (Nbar/8) * alpha) + (Pcell * (Nbar/8) * beta) )
      probs_older[i,q] <- Polder
    }
  }
#  probs_older <- probs_older/max(probs_older)
  probs_older
}






Es <- matrix(nrow = 100, ncol = 100)
Es[1:(100*100)] <- 1
Es[c(30,50),] <- .9
Es[,c(43,60)] <- .9

Ps1 <- matrix(nrow = 100,ncol = 100)
Ps1[1:(100*100)] <- 0
Ps1[45:55,45:55] <- 1

Ps2 <- matrix(nrow = 100,ncol = 100)
Ps2[1:(100*100)] <- 0
Ps2[25:35,25:35] <- 1

Ps3 <- matrix(nrow = 100,ncol = 100)
Ps3[1:(100*100)] <- 0
Ps3[40:60,30:45] <- 1

Ps4 <- matrix(nrow = 100,ncol = 100)
Ps4[1:(100*100)] <- 0
Ps4[60:90,60:80] <- 1

Ps5 <- matrix(nrow = 100,ncol = 100)
Ps5[1:(100*100)] <- 0
Ps5[10:15,20:25] <- 1


image(Ps1)
image(Ps2)
image(Ps3)
image(Ps1+Ps2+Ps3)

for (w in 1:200) {
  
  probs <- get_prob_matrix(Ps4,Es,alpha = .5,beta = 0.8)
#  step_older_Ps <- get_new_ps(Ps, probs)

  Ps4 <- probs

  probs <- get_prob_matrix(Ps5,Es,alpha = .5,beta = 0.8)
  #  step_older_Ps <- get_new_ps(Ps, probs)
  
  Ps5 <- probs
  print(w)
}
Ps4 <- Ps4/max(Ps4)
Ps5 <- Ps5/max(Ps5)
anc45 <- Ps4*Ps5
for (w in 1:400) {
  
  probs <- get_prob_matrix(Ps2,Es,alpha = .5,beta = 0.8)
  #  step_older_Ps <- get_new_ps(Ps, probs)
  
  Ps2 <- probs
  
  probs <- get_prob_matrix(Ps3,Es,alpha = .5,beta = 0.8)
  #  step_older_Ps <- get_new_ps(Ps, probs)
  
  Ps3 <- probs
  print(w)
}
Ps2 <- Ps2/max(Ps2)
Ps3 <- Ps3/max(Ps3)
anc23 <- Ps2*Ps3

for (w in 1:400) {
  
  probs <- get_prob_matrix(anc12,Es,alpha = .5,beta = 0.8)
  #  step_older_Ps <- get_new_ps(Ps, probs)
  
  anc12 <- probs
  print(w)
}

for (w in 1:700) {
  
  probs <- get_prob_matrix(Ps3,Es,alpha = .5,beta = 0.8)
  #  step_older_Ps <- get_new_ps(Ps, probs)
  
  Ps3 <- probs
  print(w)
}
anc12norm <- anc12/max(anc12)
Ps3norm <- Ps3/max(Ps3)
image(anc12)
image(Ps3norm)
image(Ps3norm*anc12norm)
anc123 <- Ps3norm*anc12norm

Ps2 <- Ps2/max(Ps2)

Ps <- Ps/max(Ps)
image(Ps)
