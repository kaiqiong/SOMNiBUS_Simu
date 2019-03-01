
# Load the results from smoothed-EM

ff <- list.files( path = "Scenario 2/Exp_With_Error/Smoothed-EM/",pattern="RData")


my.power <- NULL
my.md.gam <- NULL
my.ss.gam <- NULL
for ( i in 1:length(ff)){
  load(paste0("Scenario 2/Exp_With_Error/Smoothed-EM/",ff[i]))
  my.md.gam <- c(my.md.gam, MD[sss])
  my.power <- c(my.power, sum(sum.est[2,,2]<0.05)/100)
  my.ss.gam <- c(my.ss.gam, sss)
}




# Load the results from BiSeq and dmrseq

ff <- list.files( path = "Scenario 2/Exp_With_Error/Other_two/",pattern="RData")


diff <- c(seq(0.52, 0.54, 0.005), seq(0.545, 0.565, 0.0025) )

my.power.BS <-my.power.DM.m <- my.power.DM.f <- NULL
my.ss <- my.md <- NULL

for ( i in 1:length(ff)){
  load(paste0("Scenario 2/Exp_With_Error/Other_two/",ff[i]))
  my.ss <- c(my.ss, sss)
  my.md <- c(my.md, MD[sss])
  my.power.BS <- c(my.power.BS, sum(sum.est.BS[,2]<0.05)/100)
  my.power.DM.m <- c(my.power.DM.m,sum(sum.est.DM[,2]<0.05)/100 )
  my.power.DM.f <- c(my.power.DM.f,sum(sum.est.DM[,3]<0.05)/100 )
}




pdf(file = "Power-compare-With-Error.pdf", width = 6, height = 5)
par(mfrow = c(1, 1),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
)

diff <- c(seq(0.52, 0.54, 0.005), seq(0.545, 0.565, 0.0025) )
my.ss.gam.u = 15-my.ss.gam; my.ss.u <- 15-my.ss
plot((0.5803883-diff)[(my.ss.gam.u)[order(my.ss.gam.u)]], (my.power)[rev(order(my.ss.gam.u))], ylim = c(0, 1), pch = 19, type = "b",
     ylab = "Power", xlab = "Maximum Deviance", main = "With Error")
lines((0.5803883-diff)[my.ss.u[order(my.ss.u)]], my.power.BS[rev(order(my.ss.u))], col = "blue", pch = 19, type = "b")
lines((0.5803883-diff)[my.ss.u[order(my.ss.u)]], my.power.DM.m[rev(order(my.ss.u))], col = "red", pch = 19, type = "b")
legend("bottomright", c("smoothed-EM", "BiSeq", "dmrseq"), pch = 19, lty = 1, col = c("black", "blue", "red"), bty = "n", cex = 1)

dev.off()



# might save the points in the plot
#res.gam <- cbind( (0.5803883-diff)[(my.ss.gam.u)[order(my.ss.gam.u)]],
#                  (my.power)[rev(order(my.ss.gam.u))])
#
#res.BS <- cbind((0.5803883-diff)[my.ss.u[order(my.ss.u)]],  my.power.BS[rev(order(my.ss.u))])
#res.DM <- cbind((0.5803883-diff)[my.ss.u[order(my.ss.u)]],  my.power.DM.m[rev(order(my.ss.u))])

#colnames(res.gam) <- colnames(res.BS) <- colnames(res.DM) <- c("MD", "Power")
#setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/kaiqiong.zhao/scratch/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Manu-simulation/Plots-in-manuscript/Results")

#save(res.BS, res.DM, res.gam, file="Power-With-Error-Samp-100-Simu-100.RData")


