
# Load the functions
source("To load/Functions-V12.R")

# Load the results from the real data analysis
load("To load/BANK1betas.RData")


# The setting used in the one-region simulation experiment
beta.0 <- beta.0/4 + 2


# Settings
my.samp <- 400  
mm <- 1
set.seed(3432421+mm)
# these does not influence the settings for pi.0 and pi.1 at all

Z <- data.frame(matrix(NA, nrow= my.samp, ncol = 1)); colnames(Z) <- c( "cell_type")

Z$cell_type <- sample(c(0,1), size = my.samp, replace = T ) # simulate Z from binomial distribution
Z <-as.matrix(Z);rownames(Z)<- NULL
samp.Z <- Z

# Use bootstrap to build the read-depth matrix 
my.X <- matrix(sample(as.vector(dat.use.total), size = nrow(Z)*length(pos), replace = T) ,
               nrow = nrow(Z), ncol = length(pos))
# Data are simulated without error

sim.dat<-BSMethGammSim(n= nrow(Z), posit = pos, theta.0 =beta.0, beta= cbind(beta.2), 
                       X = my.X, Z = samp.Z,p0 = 0.003, p1 = 0.9)
# There are no randomness in the generation of pi.
pi.0 <- sim.dat$pi[which(Z==0)[1],]
pi.1 <- sim.dat$pi[which(Z==1)[1],]

diff <- c(seq(0.52, 0.54, 0.005), seq(0.545, 0.565, 0.0025) )
#diff <- seq(0.52, 0.605, 0.005) # the difference between pi.0 and pi.0.new
PI.0 <- BETA.0 <- BETA.1 <- matrix(NA, nrow = length(pos), ncol = length(diff))
MD <- rep(NA, length(diff))

bb <- c(102711743, 102712100-110)
for ( i in diff){
  
  pi.0.new <- pi.0 + i
  
  if(any(pi.0.new>1)){
    
    #if(any(pi.0.new[pos<bb[1]] > 1)){
      pi.0.new[pos<bb[1]] <- pi.1[pos<bb[1]]
    #}
    
    if(any(pi.0.new[pos>bb[2]] > 1)){
      pi.0.new[pos>bb[2]] <- pi.1[pos>bb[2]]
    }
    
    if(any(pi.0.new>=1)){
      
     temp.id <-  intersect(which(pos >= bb[1]) , which(pos <= bb[2]))
     temp.vec <- pi.0.new[temp.id]
     id <- which(temp.vec >  max(pi.1[(pos>bb[2])]))
     temp.vec[id] <- max(pi.1[(pos>bb[2])])
     pi.0.new[temp.id] <- temp.vec
     
    }
    #pi.0.new[id] <- max(pi.0.new[-id])
    beta.0.new <- log(pi.0.new/(1-pi.0.new))
    
    lo <- predict(loess(beta.0.new~ pos))
    pi.0.new <- exp(lo)/(1+exp(lo))
  }else{
    beta.0.new <- log(pi.0.new/(1-pi.0.new))
  }

  PI.0[,match(i, diff)] <- pi.0.new
  BETA.0[,match(i, diff)] <- beta.0.new
  BETA.1[,match(i, diff)] <- log((pi.1)/(1-pi.1)) - beta.0.new
  
  MD[match(i, diff)] <- max((pi.0.new - pi.1)) # 0.08
}

save(BETA.0, BETA.1, MD, PI.0,
     file ="Scenario 2/Settings_Final_MD_one_direction.RData")



# Generate the plot



pdf(file = "Power-Setting-Finner.pdf", width = 6, height = 4)
par(mfrow = c(1, 1),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
)
bb <- c(102711743, 102712100)
plot(pos[order(pos)], pi.0[order(pos)], type ="l", ylim = c(0.9,1), 
     xlab="Genomic Position", ylab =expression(pi),xaxt ="n")
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)
abline(v = bb[1], lty = 2)
abline(v = bb[2], lty = 2)
lines(pos[order(pos)], pi.1[order(pos)], type ="l", col = "red", lwd = 2, lty=4)
for( i in 1:ncol(PI.0)){
  lines(pos[order(pos)],PI.0[ order(pos), i])
}

legend("bottomright", legend = c("Z=0", "Z=1"), col = c(1,2), lty= c(1,4), bty="n", cex = 1)
dev.off()


