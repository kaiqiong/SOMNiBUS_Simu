# Pointwise coverage in the simulation study.

load("Scenario 1/Exp_With_Error/Exp-REML-knots5p10.9Samp40Simu1000EDF_Z1_Z2_Ind.RData")


em.coverage <- matrix(NA, nrow = length(beta.0), ncol = dim(sum.est)[1])

colnames(em.coverage) <- paste("beta", 0:(dim(sum.est)[1]-1), sep =".")

for ( i in 1:dim(sum.est)[1]){ # for beta.0, beta.1, beta.2 and beta.3
  beta.tr <- eval(as.name(paste("beta", i-1, sep=".")))
  beta.est <- eval(as.name(paste("beta", i-1, "est", sep=".")))
  se.est <- eval(as.name(paste("beta", i-1, "se", sep=".")))
  
  for ( m in 1:M){
    beta.est[,m] <- beta.est[match(pos,pos.index[,m]),m]
    se.est[,m] <- se.est[match(pos,pos.index[,m]),m]
  }
  
  for ( j in 1:length(beta.0)){ # for each tj
    ll <- beta.est[j,]-1.96*se.est[j,]; hh <-  beta.est[j,]+ 1.96*se.est[j,]
    em.coverage[j, i] <- sum((ll < beta.tr[j]) & (hh > beta.tr[j]))/M
  }
}

save(em.coverage, file = paste("Scenario 1/Coverage_results/CI-coverage_sample",samp.size, "Simu_", M, ".RData", sep = "" ))

#-----------------------------------------------------------

load("Scenario 1/Exp_With_Error/Exp-REML-knots5p10.9Samp100Simu1000EDF_Z1_Z2_Ind.RData")


em.coverage <- matrix(NA, nrow = length(beta.0), ncol = dim(sum.est)[1])

colnames(em.coverage) <- paste("beta", 0:(dim(sum.est)[1]-1), sep =".")

for ( i in 1:dim(sum.est)[1]){ # for beta.0, beta.1, beta.2 and beta.3
  beta.tr <- eval(as.name(paste("beta", i-1, sep=".")))
  beta.est <- eval(as.name(paste("beta", i-1, "est", sep=".")))
  se.est <- eval(as.name(paste("beta", i-1, "se", sep=".")))
  
  for ( m in 1:M){
    beta.est[,m] <- beta.est[match(pos,pos.index[,m]),m]
    se.est[,m] <- se.est[match(pos,pos.index[,m]),m]
  }
  
  for ( j in 1:length(beta.0)){ # for each tj
    ll <- beta.est[j,]-1.96*se.est[j,]; hh <-  beta.est[j,]+ 1.96*se.est[j,]
    em.coverage[j, i] <- sum((ll < beta.tr[j]) & (hh > beta.tr[j]))/M
  }
}

save(em.coverage, file = paste("Scenario 1/Coverage_results/CI-coverage_sample",samp.size, "Simu_", M, ".RData", sep = "" ))

#-----------------------------------------------------------

load("Scenario 1/Exp_With_Error/Exp-REML-knots5p10.9Samp150Simu1000EDF_Z1_Z2_Ind.RData")


em.coverage <- matrix(NA, nrow = length(beta.0), ncol = dim(sum.est)[1])

colnames(em.coverage) <- paste("beta", 0:(dim(sum.est)[1]-1), sep =".")

for ( i in 1:dim(sum.est)[1]){ # for beta.0, beta.1, beta.2 and beta.3
  beta.tr <- eval(as.name(paste("beta", i-1, sep=".")))
  beta.est <- eval(as.name(paste("beta", i-1, "est", sep=".")))
  se.est <- eval(as.name(paste("beta", i-1, "se", sep=".")))
  
  for ( m in 1:M){
    beta.est[,m] <- beta.est[match(pos,pos.index[,m]),m]
    se.est[,m] <- se.est[match(pos,pos.index[,m]),m]
  }
  
  for ( j in 1:length(beta.0)){ # for each tj
    ll <- beta.est[j,]-1.96*se.est[j,]; hh <-  beta.est[j,]+ 1.96*se.est[j,]
    em.coverage[j, i] <- sum((ll < beta.tr[j]) & (hh > beta.tr[j]))/M
  }
}

save(em.coverage, file = paste("Scenario 1/Coverage_results/CI-coverage_sample",samp.size, "Simu_", M, ".RData", sep = "" ))

#-----------------------------------------------------------

load("Scenario 1/Exp_With_Error/Exp-REML-knots5p10.9Samp400Simu1000EDF_Z1_Z2_Ind.RData")


em.coverage <- matrix(NA, nrow = length(beta.0), ncol = dim(sum.est)[1])

colnames(em.coverage) <- paste("beta", 0:(dim(sum.est)[1]-1), sep =".")

for ( i in 1:dim(sum.est)[1]){ # for beta.0, beta.1, beta.2 and beta.3
  beta.tr <- eval(as.name(paste("beta", i-1, sep=".")))
  beta.est <- eval(as.name(paste("beta", i-1, "est", sep=".")))
  se.est <- eval(as.name(paste("beta", i-1, "se", sep=".")))
  
  for ( m in 1:M){
    beta.est[,m] <- beta.est[match(pos,pos.index[,m]),m]
    se.est[,m] <- se.est[match(pos,pos.index[,m]),m]
  }
  
  for ( j in 1:length(beta.0)){ # for each tj
    ll <- beta.est[j,]-1.96*se.est[j,]; hh <-  beta.est[j,]+ 1.96*se.est[j,]
    em.coverage[j, i] <- sum((ll < beta.tr[j]) & (hh > beta.tr[j]))/M
  }
}

save(em.coverage, file = paste("Scenario 1/Coverage_results/CI-coverage_sample",samp.size, "Simu_", M, ".RData", sep = "" ))
