

P1 <-  matrix(NA, nrow = 1000, ncol = 4)
colnames(P1) <- paste0("Samp", c(40, 100, 150, 400))
P2 <- P3 <- P1
load("Scenario 1/Exp_No_Error/NoError_Three_ZSamp40Simu1000.RData")
P1[,1] <-sum.est[4,,2]
P2[,1]<-sum.est.BS[,2]
P3[,1] <-sum.est.DM[,2]
load("Scenario 1/Exp_No_Error/NoError_Three_ZSamp100Simu1000.RData")
P1[,2] <-sum.est[4,,2]
P2[,2]<-sum.est.BS[,2]
P3[,2] <-sum.est.DM[,2]
load("Scenario 1/Exp_No_Error/NoError_Three_ZSamp150Simu1000.RData")
P1[,3] <-sum.est[4,,2]
P2[,3]<-sum.est.BS[,2]
P3[,3] <-sum.est.DM[,2]
load("Scenario 1/Exp_No_Error/NoError_Three_ZSamp400Simu1000.RData")
P1[,4] <-sum.est[4,,2]
P2[,4]<-sum.est.BS[,2]
P3[,4] <-sum.est.DM[,2]


#save(P1, P2, P3, file = "TypeOnePvalue-3Z-NoError.RData")


pdf(file = "QQ-plot-3Z-NoError.pdf", width = 7, height = 5)
par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
)

pval.1 <- P1[,1]
p.exp <-  -log10((1:length(pval.1))/(length(pval.1)+1 ))
plot(p.exp, -log10(P1[,1])[(order(P1[,1]))], xlab = "-log10(Expected p values)",
     ylab = "-log10(Observed p values)", pch = 20, col = 1 , main = "N = 40")
points(p.exp, -log10(P2[,1])[(order(P2[,1]))], pch = 20, col ="blue")
points(p.exp, -log10(P3[,1])[(order(P3[,1]))], pch = 20, col ="red")
abline(a = 0, b = 1, col = 1, lty = 4)
legend("topleft", c("smoothed-EM", "BiSeq", "dmrseq"),  fill = c(1, "blue", "red"), bty = "n", cex = 1)

plot(p.exp, -log10(P1[,2])[(order(P1[,2]))], xlab = "-log10(Expected p values)",
     ylab = "-log10(Observed p values)", pch = 20, col = 1 , main = "N = 100")
points(p.exp, -log10(P2[,2])[(order(P2[,2]))], pch = 20, col ="blue")
points(p.exp, -log10(P3[,2])[(order(P3[,2]))], pch = 20, col ="red")
abline(a = 0, b = 1, col = 1, lty = 4)

plot(p.exp, -log10(P1[,3])[(order(P1[,3]))], xlab = "-log10(Expected p values)",
     ylab = "-log10(Observed p values)", pch = 20, col = 1 , main = "N = 150")
points(p.exp, -log10(P2[,3])[(order(P2[,3]))], pch = 20, col ="blue")
points(p.exp, -log10(P3[,3])[(order(P3[,3]))], pch = 20, col ="red")
abline(a = 0, b = 1, col = 1, lty = 4)

plot(p.exp, -log10(P1[,4])[(order(P1[,4]))], xlab = "-log10(Expected p values)",
     ylab = "-log10(Observed p values)", pch = 20, col = 1 , main = "N = 400")
points(p.exp, -log10(P2[,4])[(order(P2[,4]))], pch = 20, col ="blue")
points(p.exp, -log10(P3[,4])[(order(P3[,4]))], pch = 20, col ="red")
abline(a = 0, b = 1, col = 1, lty = 4)
dev.off()