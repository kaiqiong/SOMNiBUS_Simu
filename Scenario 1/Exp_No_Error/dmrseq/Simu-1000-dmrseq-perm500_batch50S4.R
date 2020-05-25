# Data are generated without error 
pacman::p_load(globaltest)
pacman::p_load(bsseq)
pacman::p_load(dmrseq)
pacman::p_load(BiSeq)
pacman::p_load(metap)

source("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Functions-V12.R")

# Load the population parameters
load("/project/6002088/zhaokq/Simulation/FunctionsToLoad/BANK1betas.RData")

beta.0 <- beta.0/4 + 2


sss <- commandArgs(TRUE)

sss <- as.numeric(sss)

samp.v <- c(40, 100, 150, 400)

my.samp <- 400

#----------------------------------------------#
n.knots = 5; my.p1 = 0.9
#----------------------------------------------#

M = 1000/100 # MC size
pos.index <-beta.0.est <- beta.1.est <- matrix(NA, nrow=length(pos), ncol = M)

GamObj <- vector("list", M)

sum.est <- array(NA, c( 4,M ,3), dimnames = list(c("s(Posit)","s(Posit):Z[, 1]","s(Posit):Z[, 2]","s(Posit):Z[, 3]"), NULL, c("chi.sq", "pvalue","pvalue_log")))

sum.est.em <- sum.est

beta.3 <- rep(0, length(beta.1))



beta.0.se <- beta.1.se <-  matrix(NA, nrow=length(pos), ncol = M)

# store results from the other two methods
pos.index.BS <-beta.1.est.BS <- beta.1.est.DM <- matrix(NA, nrow=length(pos), ncol = M)


sum.est.BS <- matrix(NA, nrow = M, ncol = 2)
colnames(sum.est.BS) <- c("Stats", "p-value")

sum.est.DM <- matrix(NA, nrow = M, ncol=3)
colnames(sum.est.DM) <- c("Mul_Regions", "MinP", "Fisher")
# array(NA, c( 2,M ,2), dimnames = list(c("s(Posit)","s(Posit):Z[, 1]"), NULL, c("chi.sq", "pvalue")))

#sum.est.em <- array(NA, c( 2,M ,3), dimnames = list(c("s(Posit)","s(Posit):Z[, 1]"), NULL, c("chi.sq", "pvalue", "logPvalue")))


beta.1.se.BS <- beta.1.se.DM <-  matrix(NA, nrow=length(pos), ncol = M)

time.0 <- Sys.time()
pval_dmrseq <- vector("list", M)





total_index <- 1:1000
sep_index <- split(total_index, ceiling(total_index/10))


for(mm in sep_index[[sss]]){

  imm <- match(mm, sep_index[[sss]])
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
  
  
  
  #--- Organize the data before EM-smooth ---#
  X <-my.X; Y <- sim.dat$Y 
  samp.size <- nrow(Y); my.p <- ncol(Y)
  my.span.dat <- data.frame(Y=as.vector(t(Y)), X=as.vector(t(X)), Posit = rep(pos, samp.size),
                            ID = rep(1:samp.size, each=my.p),
                            sapply(1:ncol(Z), function(i){rep(Z[,i], each = my.p)}))
  colnames(my.span.dat)[-(1:4)] <- colnames(Z)
  my.span.dat <- my.span.dat[my.span.dat$X>0,]
  
  
  Z <- my.span.dat[,-c(1:4)]
  
  # my.span.dat<- data.frame(my.span.dat, null = sample(c(0,1), size = nrow(my.span.dat), replace = T))
  if(ncol(Z)==1){
    Z <- matrix(my.span.dat[,-(1:4)], nrow = nrow(my.span.dat))
  }else{
    Z <- my.span.dat[,-(1:4)]
  }
  #------------------#
  # -- SmoothEM ---#
  #-----------------#
  #out <- BSMethEM(data=my.span.dat, n.k = rep(n.knots,ncol(Z)+1 ),epsilon = 10^(-6)/300,p0 = 0.003,
  #                p1 = my.p1,maxStep = 500, method="REML")
  #gam.int = out$FinalGamObj
  # -- Extract the functional parameter estimates ---#
  # -- Extract the functional parameter estimates ---#
  
  # --- Extract the (approximate) p-value and Chi.sq estimate -- these inference ignored the uncertainty in the smoothing parameter
 # sum.est.em [ ,mm, "chi.sq"] <-  out$chi.sq
 # sum.est.em [ ,mm, "pvalue"] <- out$pvalue
 # sum.est.em[, mm, "pvalue_log"] <- out$pvalue.log
  
  #
  samp.name <- paste0("ID", 1:samp.size)
  
  #---------------------#
  #--- Method dmrseq ---#
  #---------------------#
  
  see.0 <- Sys.time()
  
  BStmp <- BSseq(chr = rep("chr4", length(pos)), pos = pos, M = t(Y), Cov = t(X), sampleNames =samp.name )
  pData(BStmp)$disease <- as.vector(samp.Z[,1])
  pData(BStmp)$cell_type <- as.vector(samp.Z[,2])
  pData(BStmp)$NullZ <- as.vector(samp.Z[,3])
  regions <- dmrseq(bs = BStmp, testCovariate = "NullZ", cutoff = 0.00001, adjustCovariate = c("disease", "cell_type"),maxPerms=500)
  
  temp.p <- regions@elementMetadata@listData$pval

  imm <- match(mm, sep_index[[sss]])
    pval_dmrseq[[imm]] <- temp.p
  if (length(temp.p)>1){
    sum.est.DM[imm,1] <- 1
    sum.est.DM[imm,2] <-  min(temp.p)
    sum.est.DM[imm,3] <-  sumlog(temp.p)$p
  }
  if (length(temp.p)==1){
    sum.est.DM[imm,1] <- 0
    sum.est.DM[imm,2] <-  temp.p
    sum.est.DM[imm,3] <-  temp.p
  }
  
  print(Sys.time()-see.0)
  
  
  if(mm-trunc(mm/20)*20 ==0 ){
    print(mm)
    print(Sys.time()-time.0)
    save.image(paste0("Batch", sss, "WithError_Three_Z", "Samp", my.samp, "Simu", M, ".RData"))
  }
  
  
}
print(Sys.time()-time.0)
save.image(paste0("Batch", sss,"WithError_Three_Z", "Samp", my.samp, "Simu", M, ".RData"))


