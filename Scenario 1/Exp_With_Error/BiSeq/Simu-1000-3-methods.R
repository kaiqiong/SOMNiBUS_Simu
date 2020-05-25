# Data are generated without error 
pacman::p_load(globaltest)
#pacman::p_load(bsseq)
#pacman::p_load(dmrseq)
pacman::p_load(BiSeq)
#pacman::p_load(metap)

source("../FunctionsToLoad/Functions-V12.R")

# Load the population parameters
load("../FunctionsToLoad/BANK1betas.RData")

beta.0 <- beta.0/4 + 2


sss <- commandArgs(TRUE)

sss <- as.numeric(sss)

samp.v <- c(40, 100, 150, 400)

my.samp <- samp.v[sss]

#----------------------------------------------#
n.knots = 5; my.p1 = 0.9
#----------------------------------------------#

M = 1000 # MC size
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



beta.1.se.BS <- beta.1.se.DM <-  matrix(NA, nrow=length(pos), ncol = M)

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
                         X = my.X, Z =Z,p0 = 0.003, p1 = my.p1)
  
  
  
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
  
  #-------------#
  # -- BiSeq ---#
  #-------------#
  samp.name <- paste0("ID", 1:samp.size)
  
  
  
  # create a BSraw data type
  metadata <- list(Sequencer = "Sequencer", Year = "2017")
  rowRanges <- GRanges(seqnames = "chr4",
                       ranges = IRanges(start = pos, end = pos ))
  colData <- DataFrame(group = samp.Z,
                       row.names = samp.name )
  totalReads <- data.frame(t(X))
  methReads <- data.frame(t(Y))
  colnames(totalReads) <- colnames(methReads) <- samp.name
  rownames(totalReads) <- rownames(methReads) <- NULL
  for (j in 1:ncol(totalReads)){
    totalReads[, j] <- as.integer (totalReads[, j])
    methReads[,j] <- as.integer (methReads[,j])
  }
  
  dda <- BSraw(metadata = metadata,
               rowRanges = rowRanges,
               colData = colData,
               totalReads = as.matrix(totalReads),
               methReads = as.matrix(methReads))
  
  rowRanges(dda)$cluster.id <- "Chr4:INT15331"
  
  # Individual CpG estimates
  ind.cov <- totalReads(dda) > 0
  quant <- quantile(totalReads(dda)[ind.cov], 0.9)
  rrbs.clust.lim <- limitCov(dda, maxCov = quant)
  predictedMeth <- predictMeth(object=rrbs.clust.lim)
  betaResults <- betaRegression(formula = ~group.NullZ+group.cell_type+group.disease, link = "probit",
                                object = predictedMeth, type="BR")
  
  #cbind( 2*pnorm(abs(betaResults$estimate/betaResults$std.error), lower.tail = F), betaResults$p.val)
  
  beta.1.est.BS[,mm]<-betaResults$estimate
  beta.1.se.BS[,mm] <- betaResults$std.error
  pos.index.BS[,mm] <- betaResults$pos
  
  # Region-based statistic
  # -- averaged standardized Z score
  
  # the BiSeq region-based statistics --- combining the pvalue(transformed into a Z sore) 
  # by calculating the mean of Z score
  
  # to estimate the SD of the Z scores we have to estimate the correlation and hence the variogram of methylation between two CpG sites within a cluster.
  
  
  predictedMethNull <- predictedMeth
  
  colData(predictedMethNull)$group.NullZ <- colData(predictedMeth)$group.NullZ[sample(1:my.samp, size = my.samp, replace = T)]
  betaResultsNull <- betaRegression(formula = ~ group.NullZ + group.cell_type+group.disease, link = "probit",
                                object = predictedMethNull, type="BR")
  
  vario <- makeVariogram(betaResultsNull)
  vario.sm <- smoothVariogram(vario, sill = 0.9)
  plot(vario$variogram$v)
  lines(vario.sm$variogram[,c("h", "v.sm") ], col = "red", lwd = 1.5)
   
  ## auxiliary object to get the pValsList for the test results of interest
  
  vario.aux <- makeVariogram(betaResults, make.variogram = FALSE)
  vario.sm$pValsList <- vario.aux$pValsList
  locCor <- estLocCor(vario.sm)
  clusters.rej <- testClusters(locCor, FDR.cluster = 1)
  sum.est.BS[mm,2] <-  clusters.rej$clusters.reject$p.value
  #rrbs <- rawToRel(dda)
  #ttt <- globalTest(group.NullZ + group.disease + group.cell_type ~ group.disease + group.cell_type ,rrbs)
  
  #ttt@result[,"p-value"]
  #sum.est.BS[mm, 1]<- ttt@result[,"Statistic"]
  
  
  print(Sys.time()-time.0)
  
  
  if(mm-trunc(mm/200)*200 ==0 ){
    print(mm)
    print(Sys.time()-time.0)
    save.image(paste0("WithError_Three_Z", "Samp", my.samp, "Simu", M, ".RData"))
  }
  
  
}
print(Sys.time()-time.0)
save.image(paste0("WithError_Three_Z", "Samp", my.samp, "Simu", M, ".RData"))


