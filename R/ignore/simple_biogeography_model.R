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
      Polder <- Ecell * ( ( (1-Pcell) * Nbar * alpha) + (Pcell * Nbar * beta) )
      probs_older[i,q] <- Polder
    }
  }
  probs_older <- probs_older/max(probs_older)
  probs_older
}

get_new_ps <- function(Ps, prob_matrix) {
  new_Ps <- Ps
  for (i in 1:10) {
    for (q in 1:10) {
      probsample <- runif(1)
      if (probsample > prob_matrix[i,q]) {
        new_Ps[i,q] <- 0
      }
      else {
        new_Ps[i,q] <- 1
      }
    }
  }
  new_Ps
}



Ps <- matrix(nrow = 10,ncol = 10)
Es <- matrix(nrow = 10,ncol = 10)
Es[1:100] <- 1
#Es[1:100] <- 0
#Ps[1:100] <- 0

Es[25:70] <- 1
Ps[25:70] <- 1

Ps[,1:5] <- 0
Ps[c(1:2,7:10),] <- 0
Ps[c(2:7),c(5,8)] <- 0
Ps[c(2,7),6:7] <- 0
Ps[,9:10] <- 0
Es[,c(5,8)] <- .5







for (w in 1:200) {

probs <- get_prob_matrix(Ps,Es,alpha = .001,beta = .9999)
step_older_Ps <- get_new_ps(Ps, probs)
Ps <- step_older_Ps
print(w)
}












