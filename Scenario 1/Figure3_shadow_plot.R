
# Load the output from the simulation
load("Scenario 1/Exp_With_Error/Exp-REML-knots5p10.9Samp40Simu100EDF_Z1_Z2_Ind.RData")



# Plot type is line

pdf(file="Figure3_Simu_1_shadow-plot_sample_40Simu_100-1.pdf", width = 7, height = 5)
par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
) 
plot(pos[order(pos)], beta.0[order(pos)], col="red", xaxt ="n",  type="l", ylim=c(-1, 1.7),
     xlab="Genomic Position", ylab=" ", main = expression(paste(beta[0], "(t)")))
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)
for( i in 1:M){
  lines(  pos.index[order(pos.index[,i]),i], beta.0.est[order(pos.index[,i]),i ],pch = as.character(i),col="gray") 
}
lines(pos[order(pos)], beta.0[order(pos)], col="red")

plot(pos[order(pos)], beta.1[order(pos)],col=2, xaxt ="n", ylim = c(-1,1), type="l",
     xlab="Genomic Position", ylab=" ", main = expression(paste(beta[1], "(t)")))
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)
for( i in 1:M){
  lines(pos.index[order(pos.index[,i]),i], beta.1.est[order(pos.index[,i]),i ], pch = as.character(i), col="gray") 
}
lines(pos[order(pos)], beta.1[order(pos)], col="red")
abline(h=0, lty=4)



plot(pos[order(pos)], beta.2[order(pos)], col=2, xaxt ="n",  type="l",
     xlab="Genomic Position", ylab=" ",  main = expression(paste(beta[2], "(t)")) , ylim = c(0,5))  #ylim = c(min(beta.2.est), max(beta.2.est))
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)

axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)

for( i in 1:M){
  lines(pos.index[order(pos.index[,i]),i], beta.2.est[order(pos.index[,i]),i], pch = as.character(i), col="gray") 
}
lines(pos[order(pos)], beta.2[order(pos)], col="red")
abline(h=0, lty=4)


beta.3 = rep(0, length(pos))
plot(pos[order(pos)], beta.3[order(pos)], col=2, xaxt ="n", ylim = c(-1,1), type="l",
     xlab="Genomic Position", ylab=" ", main = expression(paste(beta[3], "(t)")))
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)

for( i in 1:M){
  lines(pos.index[order(pos.index[,i]),i], beta.3.est[order(pos.index[,i]),i], pch = as.character(i), col="gray") 
}
lines(pos[order(pos)], beta.3[order(pos)], col="red")
abline(h=0, lty=4)
dev.off()