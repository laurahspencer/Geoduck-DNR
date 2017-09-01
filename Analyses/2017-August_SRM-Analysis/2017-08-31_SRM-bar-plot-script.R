##### READ IN DATA FROM GITHUB REPO

SRM.data4plots.ordered <- as.matrix(read.csv(url("https://raw.githubusercontent.com/laurahspencer/Geoduck-DNR/master/Data/2017-08-31_SRM.data4plots.ordered.csv"), header=TRUE, stringsAsFactors=FALSE, na.strings = "NA", row.names = 1))
SRM.proteins <- as.matrix(read.csv(url("https://raw.githubusercontent.com/laurahspencer/Geoduck-DNR/master/Data/2017-08-31_SRM.proteins.csv"), header=TRUE, stringsAsFactors=FALSE, na.strings = "NA", row.names = 1))

SRM.data4plots.ordered <- as.data.frame(SRM.data4plots.ordered)

SRM.data4plots.ordered[sum(36:38),] #sum of rows 36-38 works
SRM.data4plots.ordered[sum(39:41),1] #sum of rows 39-41 does NOT work
SRM.data4plots.ordered[mean(39:41),] # #UPDATE: the mean of rows 39-41 produces values, but they are actually not the mean, and in values from row 40
SRM.data4plots.ordered[sum(39:40),] #this works
SRM.data4plots.ordered[sum(40:41),] #this works
SRM.data4plots.ordered[sum(39,41),] #this works
is.numeric(SRM.data4plots.ordered[39:41,]) #all data in all three rows are numeric
any(is.na(SRM.data4plots.ordered[39:41,])) #no NA data in these three rows


#Color = SRM.data4plots.ordered[117,]
# Symbol =SRM.data4plots.ordered[118,]

#### PLOT NON-STRINGENT MEAN DATA FOR EACH BIOLOGICAL REP., organized, & color/pattern coded by site/treatment

# HSP 90
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(1:3),], main=SRM.proteins[1,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(4:6),], main=SRM.proteins[4,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(7:9),], main=SRM.proteins[7,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# HSP 70
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(10:12),], main=SRM.proteins[10,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(13:15),], main=SRM.proteins[13,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(16:18),], main=SRM.proteins[16,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# SUPEROXIDE DISMUTASE 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(19:21),], main=SRM.proteins[19,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(22:24),], main=SRM.proteins[22,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(25:26),], main=SRM.proteins[25,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# GLYCOGEN PHOSPHORYLASE, MUSCLE FORM 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(27:29),], main=SRM.proteins[27,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(30:32),], main=SRM.proteins[30,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(33:35),], main=SRM.proteins[33,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# CYTOCHROME P450 #WTF is happening here???? why can't i sum 3 rows together that are 39 and above?!?!?!
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(36:38),], main=SRM.proteins[36,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(39:41),], main=SRM.proteins[39,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(42:44),], main=SRM.proteins[42,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# trying to troubleshoot / figure out where the problem lies
SRM.data4plots.ordered[sum(36:38),] #sum of rows 36-38 works
SRM.data4plots.ordered[sum(39:41),] #sum of rows 39-41 does NOT work
SRM.data4plots.ordered[mean(39:41),] #but the MEAN of rows 39-41 does work
SRM.data4plots.ordered[sum(39:40),] #this works
SRM.data4plots.ordered[sum(40:41),] #this works
SRM.data4plots.ordered[sum(39,41),] #this works
is.numeric(SRM.data4plots.ordered[39:41,]) #all data in all three rows are numeric
any(is.na(SRM.data4plots.ordered[39:41,])) #no NA data in these three rows

#and this weirdness happens: 
SRM.data4plots.ordered[c(38:40),] #these are rows 38-40
SRM.data4plots.ordered[sum(38:40),] #sum of rows 38-40 is NOT CORRECT 

#### THE REST OF THE BAR PLOTS WON'T WORK, ALL WITH THE SAME ERROR and similar curious behavior 

# 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(45:47),], main=SRM.proteins[45,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(48:50),], main=SRM.proteins[48,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(51:53),], main=SRM.proteins[51,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(54:56),], main=SRM.proteins[54,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(57:59),], main=SRM.proteins[57,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(60:62),], main=SRM.proteins[60,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(63:65),], main=SRM.proteins[63,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(66:68),], main=SRM.proteins[66,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(69:71),], main=SRM.proteins[69,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(72:74),], main=SRM.proteins[72,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(75:77),], main=SRM.proteins[75,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(78:80),], main=SRM.proteins[78,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(81:82),], main=SRM.proteins[81,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(83:85),], main=SRM.proteins[83,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(86:88),], main=SRM.proteins[86,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(89:91),], main=SRM.proteins[89,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(92:94),], main=SRM.proteins[92,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(95:97),], main=SRM.proteins[95,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(98:100),], main=SRM.proteins[98,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(101:103),], main=SRM.proteins[101,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(104:106),], main=SRM.proteins[104,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[sum(107:109),], main=SRM.proteins[107,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(110:112),], main=SRM.proteins[110,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[sum(113:115),], main=SRM.proteins[113,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])