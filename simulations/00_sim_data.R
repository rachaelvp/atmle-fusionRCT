sim_data <- function(p_subgroup, n, AS_interaction = FALSE) {
  
  # covariates for all participants, including subpopulation indicator (S)
  W1 <- round(rnorm(n, mean = 50, sd = 12))
  W2 <- round(rnorm(n, mean = 28, sd = 5), 2)
  S <- rbinom(n, 1, p_subgroup)
  
  # randomized treatment assignment (A)
  A <- rbinom(n, 1, 0.5)
  
  # outcome regression w nonlinearities
  base_effect <- 120 - 0.1*W1 + 0.05*W2^2 - 5*S
  if(AS_interaction){
    treatment_effect <- -10 * A - 5*A*S
  } else {
    treatment_effect <- -10 * A
  }
  
  random_effect <- rnorm(n)
  Y <- base_effect + treatment_effect + random_effect
  
  d <- data.frame(
    S = S, W1 = W1, W2 = W2, A = A, Y = Y
  )
  return(d)
}

get_psi0 <- function(p_subgroup, N = 1e6, AS_interaction = FALSE){
  d <- sim_data(p_subgroup, N, AS_interaction)
  psi0_pooled <- mean(d[d$A == 1,"Y"])-mean(d[d$A == 0,"Y"])
  
  d_S <- d[d$S == 1,]
  psi0_subgroup <- mean(d_S[d_S$A == 1,"Y"])-mean(d_S[d_S$A == 0,"Y"])
  
  psi0 <- c("psi0_subgroup"=psi0_subgroup, "psi0_pooled"=psi0_pooled)
  
  print(paste0("subgroup-specific effect: ", psi0_subgroup))
  print(paste0("pooled effect: ", psi0_pooled))
  
  return(psi0)
}

