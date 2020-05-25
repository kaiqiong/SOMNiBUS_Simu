

P1 <-  matrix(NA, nrow = 1000, ncol = 4)
colnames(P1) <- paste0("Samp", c(40, 100, 150, 400))
P4 <- P5 <- P2 <- P3 <- P1

setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/SOMNiBUS_Simu")
load("Scenario 1/Exp_No_Error/NoError_Three_ZSamp40Simu1000.RData")
P1[,1] <-sum.est[4,,2]
P2[,1]<-sum.est.BS[,2]
P3[,1] <-sum.est.DM[,2]


load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Scripts/Exp_Three_Z_Null_Z_BSmooth_No_Error/NoError_Three_ZSamp40Simu1000.RData")


P4[,1] <- sum.est.BS[,1]
P5[,1] <- sum.est.SMSC[,1]



load("Scenario 1/Exp_No_Error/NoError_Three_ZSamp100Simu1000.RData")
P1[,2] <-sum.est[4,,2]
P2[,2]<-sum.est.BS[,2]
P3[,2] <-sum.est.DM[,2]


load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Scripts/Exp_Three_Z_Null_Z_BSmooth_No_Error/NoError_Three_ZSamp100Simu1000.RData")

P4[,2] <- sum.est.BS[,1]
P5[,2] <- sum.est.SMSC[,1]



# To be added for sample size 150 and 400

load("Scenario 1/Exp_No_Error/NoError_Three_ZSamp150Simu1000.RData")
P1[,3] <-sum.est[4,,2]
P2[,3]<-sum.est.BS[,2]
P3[,3] <-sum.est.DM[,2]


load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Scripts/Exp_Three_Z_Null_Z_BSmooth_No_Error/NoError_Three_ZSamp150Simu1000.RData")

P4[,3] <- sum.est.BS[,1]
P5[,3] <- sum.est.SMSC[,1]



load("Scenario 1/Exp_No_Error/NoError_Three_ZSamp400Simu1000.RData")
P1[,4] <-sum.est[4,,2]
P2[,4]<-sum.est.BS[,2]
P3[,4] <-sum.est.DM[,2]




load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Scripts/Exp_Three_Z_Null_Z_BSmooth_No_Error/NoError_Three_ZSamp400Simu1000.RData")

P4[,4] <- sum.est.BS[,1]
P5[,4] <- sum.est.SMSC[,1]

#load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Replicate-code/TypeOnePvalue-3Z-NoError.RData")


# Update the results for BiSeq
P2_t <- P2
load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Scripts/Exp_biseq_3Z_noError_correct/NoError_Three_ZSamp40Simu1000.RData")
P2_t[,1] <- sum.est.BS[,2]
load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Scripts/Exp_biseq_3Z_noError_correct/NoError_Three_ZSamp100Simu1000.RData")
P2_t[,2] <- sum.est.BS[,2]
load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Scripts/Exp_biseq_3Z_noError_correct/NoError_Three_ZSamp150Simu1000.RData")
P2_t[,3] <- sum.est.BS[,2]
load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Scripts/Exp_biseq_3Z_noError_correct/NoError_Three_ZSamp400Simu1000.RData")
P2_t[,4] <- sum.est.BS[,2]

setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Replicate-code")
save(P1, P2_t, P3,P4, P5, P2, file = "TypeOnePvalue-3Z-NoError-correct.RData")



#------------------
# Try use QQ plots
#------------------

library(data.table)
library(ggplot2)


myqq <- function(pval.1){
p.exp <-  -log10((1:length(pval.1))/(length(pval.1)+1 ))
p.obs <- -log10(pval.1)[(order(pval.1))]
return(cbind(p.exp, p.obs))
}

P1_plot <- data.frame(cbind( rbind(myqq(P1[,1]), myqq(P1[,2]), myqq(P1[,3]), myqq(P1[,4])), N = rep( c(40, 100, 150, 400), each = 1000)))
P2_plot <- data.frame(cbind( rbind(myqq(P2_t[,1]), myqq(P2_t[,2]), myqq(P2_t[,3]), myqq(P2_t[,4])), N = rep( c(40, 100, 150, 400), each = 1000)))

P3_plot <- data.frame(cbind( rbind(myqq(P3[,1]), myqq(P3[,2]), myqq(P3[,3]), myqq(P3[,4])), N = rep( c(40, 100, 150, 400), each = 1000)))

P4_plot <- data.frame(cbind( rbind(myqq(P4[,1]), myqq(P4[,2]), myqq(P4[,3]), myqq(P4[,4])), N = rep( c(40, 100, 150, 400), each = 1000)))
P5_plot <- data.frame(cbind( rbind(myqq(P5[,1]), myqq(P5[,2]), myqq(P5[,3]), myqq(P5[,4])), N = rep( c(40, 100, 150, 400), each = 1000)))


res_plot <- data.frame(cbind(rbind(P1_plot, P2_plot, P3_plot, P4_plot, P5_plot),
                             type = rep(c("SOMNiBUS", "BiSeq", "dmrseq",  "BSmooth", "SMSC"),
                                        each = 4000)))

res_plot$type = factor(res_plot$type, levels = c("SOMNiBUS", "BiSeq", "dmrseq",  "BSmooth", "SMSC"))

trop = c("darkorange", "dodgerblue", "hotpink","limegreen", "#984ea3")
trop = c("black", "blue", "red", "darkorange", "limegreen")
library(gridExtra)
library(ggpubr)
#grid.arrange(

ggsave(file = "ggQQ-plot-3Z-NoError-trueBiseq.pdf", width = 7, height = 5)

postscript(file = "ggQQ-plot-3Z-NoError-trueBiseq.eps", width = 7, height = 5, horizontal = FALSE, onefile = FALSE,
           paper ="special", pointsize = 12)

ggarrange( 
   ggplot(data = res_plot[res_plot$N %in% c(40, 100),], 
           aes(x = p.exp, y = p.obs,  color = type)) +
        geom_point() + 
        xlab(NULL) +
        facet_grid(cols= vars(N), labeller = label_both)+
        scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))+
     #  scale_y_continuous(limits = c(0,4))+
        geom_abline(linetype = 4)+
        scale_color_manual(values = trop),
    ggplot(data = res_plot[res_plot$N %in% c(150, 400),], 
           aes(x = p.exp, y = p.obs, color = type)) +
        geom_point() +
        facet_grid(cols = vars(N), labeller = label_both)+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))+
        scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))+
        geom_abline(linetype = 4)+
       scale_color_manual(values = trop) ,
    nrow=2, common.legend = TRUE, legend ="bottom")

dev.off()


setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Replicate-code")
load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/Smooth Covariate Effect/Revision_simu/Replicate-code/TypeOnePvalue-3Z-NoError-correct.RData")

pdf(file = "QQ-plot-3Z-NoError-trueBiseq.pdf", width = 7, height = 5)

postscript(file = "QQ-plot-3Z-NoError-trueBiseq.eps", width = 7, height = 5, horizontal = FALSE, onefile = FALSE,
           paper ="special", pointsize = 12)
par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
)

pval.1 <- P1[,1]
p.exp <-  -log10((1:length(pval.1))/(length(pval.1)+1 ))
plot(p.exp, -log10(P1[,1])[(order(P1[,1]))], xlab = "-log10(Expected p values)",
     ylab = "-log10(Observed p values)", pch = 20, col = 1 , main = "N = 40", ylim = c(0,3))
points(p.exp, -log10(P2_t[,1])[(order(P2_t[,1]))], pch = 20, col ="blue")
points(p.exp, -log10(P3[,1])[(order(P3[,1]))], pch = 20, col ="red")

points(p.exp, -log10(P4[,1])[(order(P4[,1]))], pch = 20, col ="darkorange")
points(p.exp, -log10(P5[,1])[(order(P5[,1]))], pch = 20, col ="limegreen")
points(p.exp, -log10(P2[,1])[(order(P2[,1]))], pch = 20, col ="#984ea3")
abline(a = 0, b = 1, col = 1, lty = 4)
legend("topleft", c("SOMNiBUS", "GlobalTest", "dmrseq",  "BSmooth", "SMSC", "BiSeq"),  fill = c(1, "#984ea3", "red","darkorange","limegreen", "blue"), bty = "n", cex = 1)

plot(p.exp, -log10(P1[,2])[(order(P1[,2]))], xlab = "-log10(Expected p values)",
     ylab = "-log10(Observed p values)", pch = 20, col = 1 , main = "N = 100", ylim = c(0,3))
points(p.exp, -log10(P2_t[,2])[(order(P2_t[,2]))], pch = 20, col ="blue")
points(p.exp, -log10(P3[,2])[(order(P3[,2]))], pch = 20, col ="red")


points(p.exp, -log10(P4[,2])[(order(P4[,2]))], pch = 20, col ="darkorange")
points(p.exp, -log10(P5[,2])[(order(P5[,2]))], pch = 20, col ="limegreen")
points(p.exp, -log10(P2[,2])[(order(P2[,2]))], pch = 20, col ="#984ea3")
abline(a = 0, b = 1, col = 1, lty = 4)





plot(p.exp, -log10(P1[,3])[(order(P1[,3]))], xlab = "-log10(Expected p values)",
     ylab = "-log10(Observed p values)", pch = 20, col = 1 , main = "N = 150", ylim = c(0,3))
points(p.exp, -log10(P2_t[,3])[(order(P2_t[,3]))], pch = 20, col ="blue")
points(p.exp, -log10(P3[,3])[(order(P3[,3]))], pch = 20, col ="red")

points(p.exp, -log10(P4[,3])[(order(P4[,3]))], pch = 20, col ="darkorange")
points(p.exp, -log10(P5[,3])[(order(P5[,3]))], pch = 20, col ="limegreen")
points(p.exp, -log10(P2[,3])[(order(P2[,3]))], pch = 20, col ="#984ea3")
abline(a = 0, b = 1, col = 1, lty = 4)

plot(p.exp, -log10(P1[,4])[(order(P1[,4]))], xlab = "-log10(Expected p values)",
     ylab = "-log10(Observed p values)", pch = 20, col = 1 , main = "N = 400", ylim = c(0,3))
points(p.exp, -log10(P2_t[,4])[(order(P2_t[,4]))], pch = 20, col ="blue")
points(p.exp, -log10(P3[,4])[(order(P3[,4]))], pch = 20, col ="red")

points(p.exp, -log10(P4[,4])[(order(P4[,4]))], pch = 20, col ="darkorange")
points(p.exp, -log10(P5[,4])[(order(P5[,4]))], pch = 20, col ="limegreen")
points(p.exp, -log10(P2[,4])[(order(P2[,4]))], pch = 20, col ="#984ea3")
abline(a = 0, b = 1, col = 1, lty = 4)
dev.off()
