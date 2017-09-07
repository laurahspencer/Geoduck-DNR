# OPTIONAL SCRIPT #4 FOR NON-NORMALIZED DATA
# PLOTS ALL SAMPLE DATA FOR TRANSITIONS, COLOR CODED BY SITE/TREATMENT


#### PLOT ABUNDANCE FOR EACH PEPTIDE, FOR EACH SITE, TREATMENT ### 
SRM.proteins <- data.frame(protein=SRM.data[2:116,c(1,3,4)]) #protein name to each transition
SRM.proteins[,1] <- sub(" cds.*", "", SRM.proteins[,1])
SRM.data.mean.ordered <- SRM.data.mean[-1,-1]
SRM.data.mean.ordered <- as.matrix(SRM.data.mean.ordered[ , order(names(SRM.data.mean.ordered))])

CI.B.etc <- cbind(CI.b, Order=c(rep(1, times=nrow(CI.b))), Color=c(rep(2, times=nrow(CI.b))), Symbol=c(rep(20, times=nrow(CI.b))))
CI.E.etc <- cbind(CI.e, Order=c(rep(2, times=nrow(CI.e))), Color=c(rep(2, times=nrow(CI.e))), Symbol=c(rep(100, times=nrow(CI.e))))
PG.B.etc <- cbind(PG.b, Order=c(rep(3, times=nrow(PG.b))), Color=c(rep(6, times=nrow(PG.b))), Symbol=c(rep(20, times=nrow(PG.b))))
PG.E.etc <- cbind(PG.e, Order=c(rep(4, times=nrow(PG.e))), Color=c(rep(6, times=nrow(PG.e))), Symbol=c(rep(100, times=nrow(PG.e))))
WB.B.etc <- cbind(WB.b, Order=c(rep(5, times=nrow(WB.b))), Color=c(rep(3, times=nrow(WB.b))), Symbol=c(rep(20, times=nrow(WB.b))))
WB.E.etc <- cbind(WB.e, Order=c(rep(6, times=nrow(WB.e))), Color=c(rep(3, times=nrow(WB.e))), Symbol=c(rep(100, times=nrow(WB.e))))
FB.B.etc <- cbind(FB.b, Order=c(rep(7, times=nrow(FB.b))), Color=c(rep(4, times=nrow(FB.b))), Symbol=c(rep(20, times=nrow(FB.b))))
FB.E.etc <- cbind(FB.e, Order=c(rep(8, times=nrow(FB.e))), Color=c(rep(4, times=nrow(FB.e))), Symbol=c(rep(100, times=nrow(FB.e))))
samples4plots <- rbind(CI.B.etc, CI.E.etc, PG.B.etc, PG.E.etc, WB.B.etc, WB.E.etc, FB.B.etc, FB.E.etc)
library(plyr)
samples4plots$SAMPLE <- as.character(samples4plots$SAMPLE)
samples4plots <- arrange(samples4plots, samples4plots$SAMPLE)
SRM.data4plots <- rbind.data.frame(SRM.data.mean.ordered, Order=samples4plots$Order, Color=samples4plots$Color, Symbol=samples4plots$Symbol)
SRM.data4plots[which(rownames(SRM.data4plots) == 'Order'), ] <- as.numeric(SRM.data4plots[which(rownames(SRM.data4plots) == 'Order'), ]) #have to make Order row numeric before using it to sort columns
SRM.data4plots.ordered <- as.matrix(SRM.data4plots[, order(SRM.data4plots[which(rownames(SRM.data4plots) == 'Order'), ]) ]) #sort columns (samples) by designated order


#### THE FOLLOWING PLOTS DATA FOR EACH BIOLOGICAL REP (MEAN OF GOOD TECH REPS)
tail(SRM.data4plots.ordered)
SRM.data4plots.ordered
SRM.proteins
barplot(SRM.data4plots.ordered[1,], main=SRM.proteins[1,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# HSP 90
par(mfrow=c(3,3))
barplot(SRM.data4plots.ordered[1,], main=SRM.proteins[1,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[2,], main=SRM.proteins[4,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[3,], main=SRM.proteins[7,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[4,], main=SRM.proteins[1,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[5,], main=SRM.proteins[4,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[6,], main=SRM.proteins[7,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[7,], main=SRM.proteins[1,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[8,], main=SRM.proteins[4,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[9,], main=SRM.proteins[7,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
 

# # HSP 70
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(10:12),], main=SRM.proteins[10,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(13:15),], main=SRM.proteins[13,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(16:18),], main=SRM.proteins[16,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # SUPEROXIDE DISMUTASE 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(19:21),], main=SRM.proteins[19,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(22:24),], main=SRM.proteins[22,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(25:26),], main=SRM.proteins[25,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # GLYCOGEN PHOSPHORYLASE, MUSCLE FORM 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(27:29),], main=SRM.proteins[27,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(30:32),], main=SRM.proteins[30,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(33:35),], main=SRM.proteins[33,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # CYTOCHROME P450 #WTF is happening here???? why can't i sum 3 rows together that are 39 and above?!?!?!
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(36:38),], main=SRM.proteins[36,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(39:41),], main=SRM.proteins[39,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(42:44),], main=SRM.proteins[42,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(45:47),], main=SRM.proteins[45,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(48:50),], main=SRM.proteins[48,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(51:53),], main=SRM.proteins[51,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(54:56),], main=SRM.proteins[54,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(57:59),], main=SRM.proteins[57,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(60:62),], main=SRM.proteins[60,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(63:65),], main=SRM.proteins[63,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(66:68),], main=SRM.proteins[66,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(69:71),], main=SRM.proteins[69,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(72:74),], main=SRM.proteins[72,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(75:77),], main=SRM.proteins[75,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(78:80),], main=SRM.proteins[78,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(81:82),], main=SRM.proteins[81,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(83:85),], main=SRM.proteins[83,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(86:88),], main=SRM.proteins[86,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(89:91),], main=SRM.proteins[89,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(92:94),], main=SRM.proteins[92,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(95:97),], main=SRM.proteins[95,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(98:100),], main=SRM.proteins[98,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(101:103),], main=SRM.proteins[101,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(104:106),], main=SRM.proteins[104,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 
# # 
# par(mfrow=c(3,1))
# barplot(SRM.data4plots.ordered[sum(107:109),], main=SRM.proteins[107,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(110:112),], main=SRM.proteins[110,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# barplot(SRM.data4plots.ordered[sum(113:115),], main=SRM.proteins[113,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
# 

#### PLOTTING ALL TRANSITION DATA SEPARATELY IN BAR PLOTS #### 
par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[10,], main=SRM.proteins[10,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[11,], main=SRM.proteins[11,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[12,], main=SRM.proteins[12,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[13,], main=SRM.proteins[13,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[14,], main=SRM.proteins[14,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[15,], main=SRM.proteins[15,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[16,], main=SRM.proteins[16,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[17,], main=SRM.proteins[17,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[18,], main=SRM.proteins[18,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[19,], main=SRM.proteins[19,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[20,], main=SRM.proteins[20,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[21,], main=SRM.proteins[21,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[22,], main=SRM.proteins[22,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[23,], main=SRM.proteins[23,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[24,], main=SRM.proteins[24,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[25,], main=SRM.proteins[25,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[26,], main=SRM.proteins[26,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[27,], main=SRM.proteins[27,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[28,], main=SRM.proteins[28,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[29,], main=SRM.proteins[29,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[30,], main=SRM.proteins[30,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[31,], main=SRM.proteins[31,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[32,], main=SRM.proteins[32,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[33,], main=SRM.proteins[33,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[34,], main=SRM.proteins[34,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[35,], main=SRM.proteins[35,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[36,], main=SRM.proteins[36,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[37,], main=SRM.proteins[37,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[38,], main=SRM.proteins[38,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[39,], main=SRM.proteins[39,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[40,], main=SRM.proteins[40,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[41,], main=SRM.proteins[41,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[42,], main=SRM.proteins[42,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[43,], main=SRM.proteins[43,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[44,], main=SRM.proteins[44,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[45,], main=SRM.proteins[45,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[46,], main=SRM.proteins[46,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[47,], main=SRM.proteins[47,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[48,], main=SRM.proteins[48,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[49,], main=SRM.proteins[49,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[50,], main=SRM.proteins[50,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[51,], main=SRM.proteins[51,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[52,], main=SRM.proteins[52,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[53,], main=SRM.proteins[53,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[54,], main=SRM.proteins[54,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[55,], main=SRM.proteins[55,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[56,], main=SRM.proteins[56,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[57,], main=SRM.proteins[57,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[58,], main=SRM.proteins[58,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[59,], main=SRM.proteins[59,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[60,], main=SRM.proteins[60,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[61,], main=SRM.proteins[61,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[62,], main=SRM.proteins[62,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[63,], main=SRM.proteins[63,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[64,], main=SRM.proteins[64,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[65,], main=SRM.proteins[65,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[66,], main=SRM.proteins[66,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[67,], main=SRM.proteins[67,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[68,], main=SRM.proteins[68,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[69,], main=SRM.proteins[69,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[70,], main=SRM.proteins[70,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[71,], main=SRM.proteins[71,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[72,], main=SRM.proteins[72,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[73,], main=SRM.proteins[73,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[74,], main=SRM.proteins[74,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[75,], main=SRM.proteins[75,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[76,], main=SRM.proteins[76,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[77,], main=SRM.proteins[77,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[78,], main=SRM.proteins[78,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[79,], main=SRM.proteins[79,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[80,], main=SRM.proteins[80,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[83,], main=SRM.proteins[83,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[84,], main=SRM.proteins[84,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[85,], main=SRM.proteins[85,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[86,], main=SRM.proteins[86,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[87,], main=SRM.proteins[87,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[88,], main=SRM.proteins[88,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[81,], main=SRM.proteins[81,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[82,], main=SRM.proteins[82,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[89,], main=SRM.proteins[89,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[90,], main=SRM.proteins[90,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[91,], main=SRM.proteins[91,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[92,], main=SRM.proteins[92,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[93,], main=SRM.proteins[93,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[94,], main=SRM.proteins[94,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[95,], main=SRM.proteins[95,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[96,], main=SRM.proteins[96,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[97,], main=SRM.proteins[97,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[98,], main=SRM.proteins[98,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[99,], main=SRM.proteins[99,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[100,], main=SRM.proteins[100,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[101,], main=SRM.proteins[101,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[102,], main=SRM.proteins[102,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[103,], main=SRM.proteins[103,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[104,], main=SRM.proteins[104,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[105,], main=SRM.proteins[105,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[106,], main=SRM.proteins[106,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
barplot(SRM.data.mean.ordered[107,], main=SRM.proteins[107,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[108,], main=SRM.proteins[108,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[109,], main=SRM.proteins[109,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[110,], main=SRM.proteins[110,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[111,], main=SRM.proteins[111,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[112,], main=SRM.proteins[112,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[113,], main=SRM.proteins[113,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[114,], main=SRM.proteins[114,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
barplot(SRM.data.mean.ordered[115,], main=SRM.proteins[115,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

par(mfrow=c(3,3))
plot(SRM.data.mean.ordered[107,], main=SRM.proteins[107,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
plot(SRM.data.mean.ordered[108,], main=SRM.proteins[108,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
plot(SRM.data.mean.ordered[109,], main=SRM.proteins[109,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
plot(SRM.data.mean.ordered[110,], main=SRM.proteins[110,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
plot(SRM.data.mean.ordered[111,], main=SRM.proteins[111,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
plot(SRM.data.mean.ordered[112,], main=SRM.proteins[112,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
plot(SRM.data.mean.ordered[113,], main=SRM.proteins[113,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
plot(SRM.data.mean.ordered[114,], main=SRM.proteins[114,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
plot(SRM.data.mean.ordered[115,], main=SRM.proteins[115,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
