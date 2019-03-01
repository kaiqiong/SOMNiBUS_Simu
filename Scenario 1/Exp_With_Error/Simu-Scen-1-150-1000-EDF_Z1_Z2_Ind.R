# Load the functions
source("To load/Functions-V12.R")

# Load the results from the real data analysis
load("To load/BANK1betas.RData")

beta.0 <- beta.0/4 + 2

my.samp <- 150

#----------------------------------------------#
n.knots = 5; my.p1 = 0.9
#----------------------------------------------#

M = 1000 # MC size
pos.index <- beta.3.est <- beta.0.est <- beta.1.est <- beta.2.est <- matrix(NA, nrow=length(pos), ncol = M)

GamObj <- vector("list", M)

sum.est <- array(NA, c( 4,M ,2), dimnames = list(c("s(Posit)","s(Posit):Z[, 1]", 
                                                   "s(Posit):Z[, 2]",
                                                   "s(Posit):Z[, 3]"), NULL, c("chi.sq", "pvalue")))

sum.est.em <- sum.est

# Another way to get beta.0.est by multiple the alpha with the design matrix


beta.3.se <- beta.0.se <- beta.1.se <- beta.2.se <-  matrix(NA, nrow=length(pos), ncol = M)
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
  
  #Z$disease <- ifelse(Z$disease=="RA",1,0)
  #Z$cell_type <- ifelse(Z$cell_type=="TCELL",1,0)
  
  Z <-as.matrix(Z);rownames(Z)<- NULL
  
  # Use bootstrap to build the read-depth matrix 
  my.X <- matrix(sample(as.vector(dat.use.total), size = nrow(Z)*length(pos), replace = T) ,
                 nrow = nrow(Z), ncol = length(pos))
  #--- Simulate the data ---#
 
  sim.dat<-BSMethGammSim(n= nrow(Z), posit = pos, theta.0 =beta.0, beta= cbind(beta.1,beta.2), 
                         X = my.X, Z =Z,p0 = 0.003, p1 = 0.9)
  
  
  #--- Organize the data before EM-smooth ---#
  X <-my.X; Y <- sim.dat$Y 
  samp.size <- nrow(Y); my.p <- ncol(Y)
  my.span.dat <- data.frame(Y=as.vector(t(Y)), X=as.vector(t(X)), Posit = rep(pos, samp.size),
                            ID = rep(1:samp.size, each=my.p),
                            sapply(1:ncol(Z), function(i){rep(Z[,i], each = my.p)}))
  colnames(my.span.dat)[-(1:4)] <- colnames(Z)
  my.span.dat <- my.span.dat[my.span.dat$X>0,]
  
  my.span.dat<- data.frame(my.span.dat, null = sample(c(0,1), size = nrow(my.span.dat), replace = T))
  
  Z <- my.span.dat[,-c(1:4)]
# -- Fit the model use the proposed method ---#
  out <- BSMethEM(data=my.span.dat, n.k = rep(n.knots,ncol(Z)+1 ),epsilon = 10^(-6)/300,p0 = 0.003,
                  p1 = my.p1,maxStep = 500, method="REML")

  beta.0.est[,mm] <- out$Beta.out[,"Intercept"]
  beta.1.est[,mm] <- out$Beta.out[,2]
  beta.2.est[,mm] <- out$Beta.out[,3]
  beta.3.est[,mm] <- out$Beta.out[,4]
  
  pos.index[,mm] <- out$SE.pos
  
  
  beta.0.se [,mm] <-  out$SE.out[,1]
  beta.1.se [,mm] <-  out$SE.out[,2]
  beta.2.se [,mm] <-  out$SE.out[,3]
  beta.3.se [,mm] <-  out$SE.out[,4]
  
  GamObj[[mm]]<-out$FinalGamObj
  
  sum.est[,mm,"chi.sq"] <- summary(out$FinalGamObj)$s.table[,"Chi.sq"]
  sum.est[,mm,"pvalue"] <- summary(out$FinalGamObj)$s.table[,"p-value"]
  
  
  sum.est.em [ ,mm, "chi.sq"] <-  out$chi.sq
  
  sum.est.em [ ,mm, "pvalue"] <- out$pvalue

}

print(Sys.time()-time.0)

save.image(paste0("Exp-REML-knots", n.knots, "p1", my.p1, "Samp", my.samp, "Simu", M, "EDF_Z1_Z2_Ind.RData"))







#plot(pos, beta.0)
#par(mfrow=c(1,2), mar=c(4,4,1,1)) # 10 width X 3 height
#plot(pos[order(pos)], beta.0[order(pos)], type="l", xlab="Position", 
#     ylab=expression(beta[0]), xaxt ="n")
#axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
#axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)

#summary(out$FinalGamObj)$s.table


