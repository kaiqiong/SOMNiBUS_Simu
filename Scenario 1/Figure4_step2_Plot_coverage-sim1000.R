
load("Scenario 1/Coverage_results/CI-coverage_sample40Simu_1000.RData")
res1 <- em.coverage

load("Scenario 1/Coverage_results/CI-coverage_sample100Simu_1000.RData")
res2 <-em.coverage

load("Scenario 1/Coverage_results/CI-coverage_sample150Simu_1000.RData")
res3 <-em.coverage

load("Scenario 1/Coverage_results/CI-coverage_sample400Simu_1000.RData")
res4 <-em.coverage



pdf(file = "coverage.pdf", width = 7, height = 5)
par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
)
plot(pos[order(pos)], res1[order(pos),1], ylim = c(0.7, 1), pch = 20, col = 1, xlab="Genomic Position (t)", ylab="CI Coverage",
     main = expression(paste(beta[0], "(t)")), type="l", lty = 1, lwd = 1.5, xaxt ="n")
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)
lines(pos[order(pos)], res2[order(pos),1],  col = 2, lwd = 2)
lines(pos[order(pos)], res3[order(pos),1],  col = 3, lwd = 2)
lines(pos[order(pos)], res4[order(pos),1],  col = 4, lwd = 2)
abline(h = 0.95, lty = 4, lwd = 2)
legend("bottomright", legend = c("N = 40", "N=100", "N=150", "N=400"), col = 1:4, lty=1, bty="n", cex = 0.8)

plot(pos[order(pos)], res1[order(pos),2], ylim = c(0.7, 1), pch = 20, col = 1, xlab="Genomic Position (t)", ylab="CI Coverage",
     main = expression(paste(beta[1], "(t)")), type="l", lty = 1, lwd = 1.5, xaxt ="n")
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)
lines(pos[order(pos)], res2[order(pos),2],  col = 2, lwd = 2)
lines(pos[order(pos)], res3[order(pos),2],  col = 3, lwd = 2)
lines(pos[order(pos)], res4[order(pos),2],  col = 4, lwd = 2)
abline(h = 0.95, lty = 4, lwd = 2)
legend("bottomright", legend = c("N = 40", "N=100", "N=150", "N=400"), col = 1:4, lty=1, bty="n", cex = 0.8)


plot(pos[order(pos)], res1[order(pos),3], ylim = c(0.7, 1), pch = 20, col = 1, xlab="Genomic Position (t)", ylab="CI Coverage",
     main = expression(paste(beta[2], "(t)")), type="l", lty = 1, lwd = 1.5, xaxt ="n")
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)
lines(pos[order(pos)], res2[order(pos),3],  col = 2, lwd = 2)
lines(pos[order(pos)], res3[order(pos),3],  col = 3, lwd = 2)
lines(pos[order(pos)], res4[order(pos),3],  col = 4, lwd = 2)
abline(h = 0.95, lty = 4, lwd = 2)
legend("bottomright", legend = c("N = 40", "N=100", "N=150", "N=400"), col = 1:4, lty=1, bty="n", cex = 0.8)


plot(pos[order(pos)], res1[order(pos),4], ylim = c(0.7, 1), pch = 20, col = 1, xlab="Genomic Position (t)", ylab="CI Coverage",
     main = expression(paste(beta[3], "(t)")), type="l", lty = 1, lwd = 1.5, xaxt ="n")
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)
lines(pos[order(pos)], res2[order(pos),4],  col = 2, lwd = 2)
lines(pos[order(pos)], res3[order(pos),4],  col = 3, lwd = 2)
lines(pos[order(pos)], res4[order(pos),4],  col = 4, lwd = 2)
abline(h = 0.95, lty = 4, lwd = 2)
legend("bottomright", legend = c("N = 40", "N=100", "N=150", "N=400"), col = 1:4, lty=1, bty="n", cex = 0.8)

dev.off()