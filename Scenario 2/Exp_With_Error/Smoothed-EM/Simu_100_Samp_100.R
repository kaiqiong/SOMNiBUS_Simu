

source("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Functions-V12.R")
load("/project/6002088/zhaokq/Simulation/FunctionsToLoad/BANK1betas.RData")
# Load the population parameters
load("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Settings_Final_MD_one_direction.RData")

sss<-SLURM_ARRAY_TASK_ID <- commandArgs(TRUE)

sss <- as.numeric(sss)

beta.0 <- BETA.0[,sss]; beta.1 <- BETA.1[,sss]

my.samp <- 100

#----------------------------------------------#
n.knots = 5; my.p1 = 0.9
#----------------------------------------------#

M = 100 # MC size
pos.index <-beta.0.est <- beta.1.est <- matrix(NA, nrow=length(pos), ncol = M)

GamObj <- vector("list", M)

sum.est <- array(NA, c( 2,M ,2), dimnames = list(c("s(Posit)","s(Posit):Z[, 1]"), NULL, c("chi.sq", "pvalue")))

sum.est.em <- array(NA, c( 2,M ,3), dimnames = list(c("s(Posit)","s(Posit):Z[, 1]"), NULL, c("chi.sq", "pvalue", "logPvalue")))


beta.0.se <- beta.1.se <-  matrix(NA, nrow=length(pos), ncol = M)
time.0 <- Sys.time()

save.image(file = "see.RData")
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
                         X = my.X, Z =Z,p0 = 0.003, p1 = my.p1)
  
  
  #--- Organize the data before EM-smooth ---#
  X <-my.X; Y <- sim.dat$Y 
  samp.size <- nrow(Y); my.p <- ncol(Y)
  my.span.dat <- data.frame(Y=as.vector(t(Y)), X=as.vector(t(X)), Posit = rep(pos, samp.size),
                            ID = rep(1:samp.size, each=my.p),
                            sapply(1:ncol(Z), function(i){rep(Z[,i], each = my.p)}))
  colnames(my.span.dat)[-(1:4)] <- colnames(Z)
  my.span.dat <- my.span.dat[my.span.dat$X>0,]
  
  # my.span.dat<- data.frame(my.span.dat, null = sample(c(0,1), size = nrow(my.span.dat), replace = T))
  if(ncol(Z)==1){
    Z <- matrix(my.span.dat[,-(1:4)], nrow = nrow(my.span.dat))
  }else{
    Z <- my.span.dat[,-(1:4)]
  }
  
  # -- Fit the model use the proposed method ---#
  out <- BSMethEM(data=my.span.dat, n.k = rep(n.knots,ncol(Z)+1 ),epsilon = 10^(-6)/300,p0 = 0.003,
                  p1 =my.p1,maxStep = 500, method="REML", detail=F)
  # -- Extract the functional parameter estimates ---#
  
  beta.0.est[,mm] <- out$Beta.out[,"Intercept"]
  beta.1.est[,mm] <- out$Beta.out[,2]
  
  pos.index[,mm] <- out$SE.pos
  
  beta.0.se [,mm] <-  out$SE.out[,1]
  beta.1.se [,mm] <- out$SE.out[,2]
  
  
  GamObj[[mm]]<-out$FinalGamObj
  
  sum.est[,mm,"chi.sq"] <- summary(out$FinalGamObj)$s.table[,"Chi.sq"]
  sum.est[,mm,"pvalue"] <- summary(out$FinalGamObj)$s.table[,"p-value"]
  
  
  sum.est.em [ ,mm, "chi.sq"] <-  out$chi.sq
  sum.est.em [ ,mm, "pvalue"] <- out$pvalue
  sum.est.em [ ,mm, "logPvalue"] <- out$pvalue.log
  
}
print(Sys.time()-time.0)


save.image(paste0 ("S", sss, "MD_", round(MD[sss]*100,2),  "Exp-REML-knots", n.knots, "p1", my.p1, "Samp", my.samp, "Simu", M, "EDF.RData"))







#plot(pos, beta.0)
#par(mfrow=c(1,2), mar=c(4,4,1,1)) # 10 width X 3 height
#plot(pos[order(pos)], beta.0[order(pos)], type="l", xlab="Position", 
#     ylab=expression(beta[0]), xaxt ="n")
#axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
#axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)

#summary(out$FinalGamObj)$s.table


