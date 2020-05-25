


source("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Functions-V12.R")
load("/project/6002088/zhaokq/Simulation/FunctionsToLoad/BANK1betas.RData")
# Load the population parameters
load("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Settings_Final_MD_one_direction.RData")



source("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Weights.R")
source("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Smooth.MSC.R")

#source("~/scratch/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Calculate SE/SEcalculation/Functions-V12.R")

#load("~/scratch/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Region near BANK1/BANK1betas.RData")
# Load the population parameters
#load("~/scratch/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/SOMNiBUS_Simu/Scenario 2/Settings_Final_MD_one_direction.RData")


sss <- commandArgs(TRUE)

sss <- as.numeric(sss)

beta.0 <- BETA.0[,sss]; beta.1 <- BETA.1[,sss]

my.samp <- 100

#----------------------------------------------#
n.knots = 5; my.p1 = 0.9
#----------------------------------------------#

M = 100 # MC size
pos.index.BS <- pos.index.DM <-beta.1.est.BS <- beta.1.est.DM <- matrix(NA, nrow=length(pos), ncol = M)

GamObj <- vector("list", M)

sum.est.BS <- matrix(NA, nrow = M, ncol = 2)
colnames(sum.est.BS) <- c("Stats", "p-value")

sum.est.DM <- matrix(NA, nrow = M, ncol=3)
colnames(sum.est.DM) <- c("Mul_Regions", "MinP", "Fisher")
# array(NA, c( 2,M ,2), dimnames = list(c("s(Posit)","s(Posit):Z[, 1]"), NULL, c("chi.sq", "pvalue")))

#sum.est.em <- array(NA, c( 2,M ,3), dimnames = list(c("s(Posit)","s(Posit):Z[, 1]"), NULL, c("chi.sq", "pvalue", "logPvalue")))

pval_dmrseq <- vector("list", M)

beta.1.se.BS <- beta.1.se.DM <-  matrix(NA, nrow=length(pos), ncol = M)
time.0 <- Sys.time()

k=20
weights=compute.weights(pos,length(pos),k)

sum.est.SMSC <- matrix(NA, nrow = M, ncol =1)
sum.est.BS <- sum.est.SMSC

B = 1000

see0 <- Sys.time()
for(mm in 1: M){
  set.seed(3432421+mm)
  
  #------------------- -------------------------------------------#
  Z <- data.frame(matrix(NA, nrow= my.samp, ncol = 1)); colnames(Z) <- c( "cell_type")
  
  Z$cell_type <- sample(c(0,1), size = my.samp, replace = T ) # simulate Z from binomial distribution
  Z <-as.matrix(Z);rownames(Z)<- NULL
  samp.Z <- Z
  
  # Use bootstrap to build the read-depth matrix 
  my.X <- matrix(sample(as.vector(dat.use.total), size = nrow(Z)*length(pos), replace = T) ,
                 nrow = nrow(Z), ncol = length(pos))
  #--- Simulate the data ---#
  
  sim.dat<-BSMethGammSim(n= nrow(Z), posit = pos, theta.0 =as.matrix(beta.0, ncol=1), beta= cbind(beta.1), 
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
  
  
  stat.1=apply(pi.hat.smsc[Z==1,],2,mean)-apply(pi.hat.smsc[Z==0,],2,mean)
  stat.2=apply(pi.hat.bs[Z==1,],2,mean, na.rm = T)-apply(pi.hat.bs[Z==0,],2,mean, na.rm = T)
  
  stat.1.aggr=sum(stat.1^2)
  stat.2.aggr=sum(stat.2^2)
  
  stat.1.aggr.boot=0
  stat.2.aggr.boot=0
  
  for(b in 1:B)
  {
    indices=sample(samp.size)
    Z_now <- Z[indices]
    
    stat.1.boot=apply(pi.hat.smsc[Z_now==1,],2,mean)-apply(pi.hat.smsc[Z_now==0,],2,mean)
    stat.2.boot=apply(pi.hat.bs[Z_now==1,],2,mean, na.rm = T)-apply(pi.hat.bs[Z_now==0,],2,mean, na.rm = T)
    
    stat.1.aggr.boot[b]=sum(stat.1.boot^2)
    stat.2.aggr.boot[b]=sum(stat.2.boot^2)
  }
  
  sum.est.SMSC[mm,]=mean(stat.1.aggr.boot>stat.1.aggr)
  sum.est.BS[mm,]=mean(stat.2.aggr.boot>stat.2.aggr)
  
  
  
  
  print(Sys.time() - see0)
  
  
  
  
}



save.image(paste0 ("S", sss, "MD_", round(MD[sss]*100,2),  "Samp", my.samp, "Simu", M, "EDF.RData"))







#plot(pos, beta.0)
#par(mfrow=c(1,2), mar=c(4,4,1,1)) # 10 width X 3 height
#plot(pos[order(pos)], beta.0[order(pos)], type="l", xlab="Position", 
#     ylab=expression(beta[0]), xaxt ="n")
#axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
#axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)

#summary(out$FinalGamObj)$s.table


