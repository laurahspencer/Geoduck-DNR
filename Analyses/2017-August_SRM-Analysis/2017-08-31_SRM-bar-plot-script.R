##### READ IN DATA FROM GITHUB REPO

SRM.data4plots.ordered <- as.matrix(read.csv(url("https://raw.githubusercontent.com/laurahspencer/Geoduck-DNR/master/Data/2017-08-31_SRM.data4plots.ordered.csv"), header=TRUE, stringsAsFactors=FALSE, na.strings = "NA", row.names = 1))
SRM.proteins <- as.matrix(read.csv(url("https://raw.githubusercontent.com/laurahspencer/Geoduck-DNR/master/Data/2017-08-31_SRM.proteins.csv"), header=TRUE, stringsAsFactors=FALSE, na.strings = "NA", row.names = 1))

#### PLOT NON-STRINGENT MEAN DATA FOR EACH BIOLOGICAL REP., organized, & color/pattern coded by site/treatment

# HSP 90
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[c(1:3),], main=SRM.proteins[1,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[c(4:6),], main=SRM.proteins[4,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[c(7:9),], main=SRM.proteins[7,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# HSP 70
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[10:12,], main=SRM.proteins[10,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[13:15,], main=SRM.proteins[13,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[16:18,], main=SRM.proteins[16,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# SUPEROXIDE DISMUTASE 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[19:21,], main=SRM.proteins[19,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[22:24,], main=SRM.proteins[22,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[25:26,], main=SRM.proteins[25,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# GLYCOGEN PHOSPHORYLASE, MUSCLE FORM 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[27:29,], main=SRM.proteins[27,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[30:32,], main=SRM.proteins[30,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[33:35,], main=SRM.proteins[33,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# CYTOCHROME P450 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[36:38,], main=SRM.proteins[36,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[39:41,], main=SRM.proteins[39,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[42:44,], main=SRM.proteins[42,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# CATALASE 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[45:47,], main=SRM.proteins[45,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[48:50,], main=SRM.proteins[48,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[51:53,], main=SRM.proteins[51,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# PEROXIREDOXIN
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[54:56,], main=SRM.proteins[54,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[57:59,], main=SRM.proteins[57,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[60:62,], main=SRM.proteins[60,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# TRIFUNCTIONAL ENZYME SUBUNIT BETA, MITOCHONDRIAL 
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[63:65,], main=SRM.proteins[63,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[66:68,], main=SRM.proteins[66,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[69:71,], main=SRM.proteins[69,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# RAS-RELATED PROTEIN RAB-11B
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[72:74,], main=SRM.proteins[72,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[75:77,], main=SRM.proteins[75,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[78:80,], main=SRM.proteins[78,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

#  NA/K TRANSPORTING ATPASE SUBUNIT ALPHA-4
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[81:82,], main=SRM.proteins[81,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[83:85,], main=SRM.proteins[83,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[86:88,], main=SRM.proteins[86,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# PUROMYCIN-SENSITIVE AMINOPEPTIDASE
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[89:91,], main=SRM.proteins[89,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[92:94,], main=SRM.proteins[92,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[95:97,], main=SRM.proteins[95,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# PROTEIN DISULFIDE ISOMERASE
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[98:100,], main=SRM.proteins[98,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[101:103,], main=SRM.proteins[101,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[104:106,], main=SRM.proteins[104,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# ARACHIDONATE 5-LIPOXYGENASE
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[107:109,], main=SRM.proteins[107,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[110:112,], main=SRM.proteins[110,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[113:115,], main=SRM.proteins[113,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

# RAS-RELATED PROTEIN RAB-11B
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[72,], main=SRM.proteins[72,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[75,], main=SRM.proteins[75,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[78,], main=SRM.proteins[78,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

#  NA/K TRANSPORTING ATPASE SUBUNIT ALPHA-4
par(mfrow=c(3,1))
barplot(SRM.data4plots.ordered[81,], main=SRM.proteins[81,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[83,], main=SRM.proteins[83,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])
barplot(SRM.data4plots.ordered[86,], main=SRM.proteins[86,c(1:2)], col=SRM.data4plots.ordered[117,], density=SRM.data4plots.ordered[118,])

