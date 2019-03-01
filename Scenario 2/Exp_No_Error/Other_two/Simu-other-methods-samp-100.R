pacman::p_load(globaltest)
pacman::p_load(bsseq)
pacman::p_load(dmrseq)
pacman::p_load(BiSeq)
pacman::p_load(metap)

source("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Functions-V12.R")
load("/project/6002088/zhaokq/Simulation/FunctionsToLoad/BANK1betas.RData")
# Load the population parameters
load("/project/6002088/zhaokq/Simulation/FunctionsToLoad/Settings_Final_MD_one_direction.RData")

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


beta.1.se.BS <- beta.1.se.DM <-  matrix(NA, nrow=length(pos), ncol = M)
time.0 <- Sys.time()
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
  
  #---------------------#
  #--- Method Biseq ---#
  #---------------------#
  
  see0<-Sys.time()
  
  # create a BSraw data type
  metadata <- list(Sequencer = "Sequencer", Year = "2017")
  rowRanges <- GRanges(seqnames = "chr4",
                       ranges = IRanges(start = pos, end = pos ))
  colData <- DataFrame(group = samp.Z,
                       row.names = samp.name )
  totalReads <- data.frame(t(X))
  methReads <- data.frame(t(Y))
  colnames(totalReads) <- colnames(methReads) <- NULL
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
  betaResults <- betaRegression(formula = ~group.cell_type, link = "probit",
                                object = predictedMeth, type="BR")
  
  #cbind( 2*pnorm(abs(betaResults$estimate/betaResults$std.error), lower.tail = F), betaResults$p.val)
  
  beta.1.est.BS[,mm]<-betaResults$estimate
  beta.1.se.BS[,mm] <- betaResults$std.error
  pos.index.BS[,mm] <- betaResults$pos
  
  # Region-based statistic
  rrbs <- rawToRel(dda)
  ttt <- globalTest(group.cell_type~1,rrbs)
  
  sum.est.BS[mm,2] <-  ttt@result[,"p-value"]
  sum.est.BS[mm, 1]<- ttt@result[,"Statistic"]
  print(Sys.time() - see0)
  #---------------------#
  #--- Method dmrseq ---#
  #---------------------#
  see0<-Sys.time()
  BStmp <- BSseq(chr = rep("chr4", length(pos)), pos = pos, M = t(Y), Cov = t(X), sampleNames =samp.name )
  pData(BStmp)$CellType <- as.vector(samp.Z)
  regions <- dmrseq(bs = BStmp, testCovariate = "CellType", cutoff = 0.00001)
  
  temp.p <- regions@elementMetadata@listData$pval
  if (length(temp.p)>1){
    sum.est.DM[mm,1] <- 1
    sum.est.DM[mm,2] <-  min(temp.p)
    sum.est.DM[mm,3] <-  sumlog(temp.p)$p
  }
  if (length(temp.p)==1){
    sum.est.DM[mm,1] <- 0
    sum.est.DM[mm,2] <-  temp.p
    sum.est.DM[mm,3] <-  temp.p
  }
    print(Sys.time() - see0)
  
}
print(Sys.time()-time.0)


save.image(paste0 ("S", sss, "MD_", round(MD[sss]*100,2),  "Samp", my.samp, "Simu", M, "EDF.RData"))







#plot(pos, beta.0)
#par(mfrow=c(1,2), mar=c(4,4,1,1)) # 10 width X 3 height
#plot(pos[order(pos)], beta.0[order(pos)], type="l", xlab="Position", 
#     ylab=expression(beta[0]), xaxt ="n")
#axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
#axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)

#summary(out$FinalGamObj)$s.table


