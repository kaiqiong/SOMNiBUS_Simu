
library(mgcv)

source("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Functions_better_approx_Quasi_error_Rand_effect.R")
load("/project/6002088/zhaokq/Simulation/FunctionsToLoad/BANK1betas.RData")
# Load the population parameters
load("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Settings_Final_MD_one_direction.RData")


source("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Weights.R")
source("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Smooth.MSC.R")


beta.0 <- beta.0/4 + 2

my.samp <- 100

#----------------------------------------------#
n.knots = 5; #my.p1 = 0.9
#----------------------------------------------#

M = 1000 # MC size


beta.3 <- rep(0, length(beta.1))




k=20
weights=compute.weights(pos,length(pos),k)

sum.est.SMSC <- matrix(NA, nrow = M, ncol =1)
sum.est.BS <- sum.est.SMSC

B = 1000

time.0 <- Sys.time()
for(mm in 1: M){
  set.seed(3432421+mm)
  
  #------------------- Use bootstrap to build the Z matrix-------------------------------------------#
  Z <- data.frame(matrix(NA, nrow= my.samp, ncol = 2)); colnames(Z) <- c("disease", "cell_type")
  
  #Z$disease <- sample(my.pheno$disease, size = my.samp, replace = T)
  #Z$cell_type<- sample(my.pheno$cell_type, size = my.samp, replace = T)
  #table(my.pheno$disease)
  
  Z$disease <- sample(c(0,1), size = my.samp, replace = T, 
                      prob = c(sum(my.pheno$disease =="control"),sum(my.pheno$disease =="RA") )/nrow(my.pheno) )
  Z$cell_type <- sample(c(0,1), size = my.samp, replace = T, 
                        prob = c(sum(my.pheno$cell_type =="MONO"),sum(my.pheno$cell_type =="TCELL") )/nrow(my.pheno) )
  Z$NullZ <-  sample(c(0,1), size = nrow(Z), replace = T)
  #Z$disease <- ifelse(Z$disease=="RA",1,0)
  #Z$cell_type <- ifelse(Z$cell_type=="TCELL",1,0)
  
  Z <-as.matrix(Z);rownames(Z)<- NULL
  
  samp.Z <- Z
  # Use bootstrap to build the read-depth matrix 
  my.X <- matrix(sample(as.vector(dat.use.total), size = nrow(Z)*length(pos), replace = T) ,
                 nrow = nrow(Z), ncol = length(pos))
  #--- Simulate the data ---#
  
  sim.dat<-BSMethGammSim(n= nrow(Z), posit = pos, theta.0 =beta.0, beta= cbind(beta.1,beta.2, beta.3), 
                         X = my.X, Z =Z,p0 = 0, p1 = 1)
  
  
  X <-my.X; Y <- sim.dat$Y 
  samp.size <- nrow(Y); my.p <- ncol(Y)
  samp.name <- paste0("ID", 1:samp.size)
  
  #------------------------------------#
  #--- use Karim's code for BSmooth ---#
  #------------------------------------#
  
  pi.hat.smsc = matrix(0, ncol = my.p, nrow = samp.size)
  pi.hat.bs = pi.hat.smsc
  for ( ii in 1:samp.size){
    dat_now <- data.frame(sim.dat$S[ii,],X[ii,], sim.dat$Y[ii,])
    colnames(dat_now) <- c("S", "CT", "Ccount")
    
    pi.hat.smsc[ii,] <- Smooth.MSC(dat_now, weights,error = 10^(-6) )$pi
    pi.hat.bs[ii,] <- BSmooth_r(dat_now, weights)
  }
  
  
  stat.1=apply(pi.hat.smsc[Z[,3]==1,],2,mean)-apply(pi.hat.smsc[Z[,3]==0,],2,mean)
  stat.2=apply(pi.hat.bs[Z[,3]==1,],2,mean, na.rm = T)-apply(pi.hat.bs[Z[,3]==0,],2,mean, na.rm = T)
  
  stat.1.aggr=sum(stat.1^2)
  stat.2.aggr=sum(stat.2^2)
  
  stat.1.aggr.boot=0
  stat.2.aggr.boot=0
  
  for(b in 1:B)
  {
    indices=sample(samp.size)
    Z_now <- Z[indices,3]
    
    stat.1.boot=apply(pi.hat.smsc[Z_now==1,],2,mean)-apply(pi.hat.smsc[Z_now==0,],2,mean)
    stat.2.boot=apply(pi.hat.bs[Z_now==1,],2,mean, na.rm = T)-apply(pi.hat.bs[Z_now==0,],2,mean, na.rm = T)
    
    stat.1.aggr.boot[b]=sum(stat.1.boot^2)
    stat.2.aggr.boot[b]=sum(stat.2.boot^2)
  }
  
  sum.est.SMSC[mm,]=mean(stat.1.aggr.boot>stat.1.aggr)
  sum.est.BS[mm,]=mean(stat.2.aggr.boot>stat.2.aggr)
  
  
  
 print(Sys.time() - time.0)
  
  
  if(mm-trunc(mm/200)*200 ==0 ){
    print(mm)
    print(Sys.time()-time.0)
    save.image(paste0("NoError_Three_Z", "Samp", my.samp, "Simu", M, ".RData"))
    
  }
}
print(Sys.time()-time.0)
save.image(paste0("NoError_Three_Z", "Samp", my.samp, "Simu", M, ".RData"))

