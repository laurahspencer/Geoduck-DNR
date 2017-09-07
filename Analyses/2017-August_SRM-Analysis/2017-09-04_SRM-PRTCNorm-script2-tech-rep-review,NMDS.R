# Script #2 of PRTC-Normalized SRM data analysis


#### CREATE NMDS PLOT ########

#Load the source file for the biostats package, biostats.R script must be saved in working directory

source("References/biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
library(vegan)

#Transpose the file so that rows and columns are switched 
SRM.data.t <- t(SRM.data.normalized) # t() function transposes

#Replace NA cells with 0; metaMDS() does not handle NA's
SRM.data.t.noNA <- SRM.data.t
SRM.data.t.noNA[is.na(SRM.data.t.noNA)] <- 0
head(SRM.data.t.noNA)

#Make MDS dissimilarity matrix
#
SRM.nmds <- metaMDS(SRM.data.t.noNA, distance = 'bray', k = 2, trymax = 3000, autotransform = FALSE)
# comm= your data.frame or matrix
# distance= bray, (not sure what this means)
# k= # of dimensions to assess
# trymax = max # iterations to attempt if no solution is reached

# Create Shepard plot, which shows scatter around the regression between the interpoint distances in the final configuration (i.e., the distances between each pair of communities) against their original dissimilarities.
stressplot(SRM.nmds) 

#Make figure

plot(SRM.nmds)
# site (sample) in black circle
# species (variable) in red ticks

# make figure with sample annotations https://stat.ethz.ch/pipermail/r-sig-ecology/2011-September/002371.html
SRM.nmds.samples <- scores(SRM.nmds, display = "sites")
SRM.nmds.transitions <- scores(SRM.nmds, display = "species")
SRM.nmds.samples.sorted <- SRM.nmds.samples[ order(row.names(SRM.nmds.samples)), ]
rownames(SRM.nmds.samples.sorted)

# Generate 50 distint color ID's in a vector for plotting NMDS data
library(RColorBrewer)
n <- 50
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
colors <- sample(col_vector, 50)

### PLOTTING ALL REPS BY SAMPLE NUMBER & TREATMENT ### 
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-1.5,3), ylim=c(-0.5,0.5), asp=NA)
# symbol key
# 15 = eelgrass = filled square
# 21 = bare = open circle

points(SRM.nmds.samples.sorted[c("G001-A", "G001-B"),], col=colors[1], pch=15)
points(SRM.nmds.samples.sorted[c("G002-A", "G002-B", "G002-C"),], col=colors[50], pch=15)
points(SRM.nmds.samples.sorted[c("G003-A", "G003-B"),], col=colors[3], pch=15) #G003-C is very different
points(SRM.nmds.samples.sorted[c("G007-A", "G007-B"),], col=colors[4], pch=15)
points(SRM.nmds.samples.sorted[c("G008-A", "G008-B"),], col=colors[5], pch=15)
points(SRM.nmds.samples.sorted[c("G009-A", "G009-B"),], col=colors[6], pch=15)
points(SRM.nmds.samples.sorted[c("G012-A", "G012-B", "G012-C"),], col=colors[7], pch=21)
points(SRM.nmds.samples.sorted[c("G013-A", "G013-C"),], col=colors[8], pch=21)
points(SRM.nmds.samples.sorted[c("G014-A", "G014-B"),], col=colors[9], pch=21)
points(SRM.nmds.samples.sorted[c("G015-A", "G015-B"),], col=colors[10], pch=21)
points(SRM.nmds.samples.sorted[c("G016-A", "G016-B", "G016-C"),], col=colors[11], pch=21) #not good; three pretty distant, one definitely more distant
points(SRM.nmds.samples.sorted[c("G017-A", "G017-B"),], col=colors[12], pch=21) #not great
points(SRM.nmds.samples.sorted[c("G031-A", "G031-B", "G031-C"),], col=colors[13], pch=15)
points(SRM.nmds.samples.sorted[c("G032-A", "G032-B"),], col=colors[14], pch=15)
points(SRM.nmds.samples.sorted[c("G040-A", "G040-B"),], col=colors[15], pch=21)
points(SRM.nmds.samples.sorted[c("G041-A", "G041-B"),], col=colors[16], pch=21)
points(SRM.nmds.samples.sorted[c("G042-A", "G042-B", "G042-C"),], col=colors[17], pch=21) #G042-C very off
points(SRM.nmds.samples.sorted[c("G043-A", "G043-B"),], col=colors[18], pch=21)
points(SRM.nmds.samples.sorted[c("G045-A", "G045-B"),], col=colors[19], pch=15)
points(SRM.nmds.samples.sorted[c("G047-A", "G047-B"),], col=colors[20], pch=15)
points(SRM.nmds.samples.sorted[c("G049-A", "G049-B"),], col=colors[21], pch=15)
points(SRM.nmds.samples.sorted[c("G053-A", "G053-B", "G053-remake-C", "G053-remake-D"),], col=colors[22], pch=15) #B not good
points(SRM.nmds.samples.sorted[c("G054-A", "G054-B"),], col=colors[23], pch=15)
points(SRM.nmds.samples.sorted[c("G055-A", "G055-B"),], col=colors[24], pch=15) #one is very off
points(SRM.nmds.samples.sorted[c("G057-A", "G057-B", "G057-C"),], col=colors[25], pch=21) #one is off
points(SRM.nmds.samples.sorted[c("G060-A", "G060-B"),], col=colors[26], pch=21) 
points(SRM.nmds.samples.sorted[c("G062-B", "G062-C"),], col=colors[27], pch=21)
points(SRM.nmds.samples.sorted[c("G064-A", "G064-B"),], col=colors[50], pch=21)
points(SRM.nmds.samples.sorted[c("G066-A", "G066-B"),], col=colors[29], pch=21)
points(SRM.nmds.samples.sorted[c("G070-A", "G070-B", "G070-C"),], col=colors[30], pch=21)
points(SRM.nmds.samples.sorted[c("G071-A-A", "G071-A-B"),], col=colors[31], pch=15)
points(SRM.nmds.samples.sorted[c("G071-B-A", "G071-B-B"),], col=colors[32], pch=15)
points(SRM.nmds.samples.sorted[c("G073-A", "G073-B", "G073-C"),], col=colors[33], pch=15) #all three very off
points(SRM.nmds.samples.sorted[c("G074-A", "G074-B"),], col=colors[34], pch=15)
points(SRM.nmds.samples.sorted[c("G079-A", "G079-B"),], col=colors[35], pch=21)
points(SRM.nmds.samples.sorted[c("G081-A", "G081-B"),], col=colors[36], pch=21)
points(SRM.nmds.samples.sorted[c("G104-A", "G104-B", "G104-remake-C", "G104-remake-D"),], col=colors[37], pch=21) #B
points(SRM.nmds.samples.sorted[c("G105-A", "G105-B"),], col=colors[38], pch=21)
points(SRM.nmds.samples.sorted[c("G109-A", "G109-B", "G109-C"),], col=colors[39], pch=15)
points(SRM.nmds.samples.sorted[c("G110-A", "G110-B"),], col=colors[40], pch=15)
points(SRM.nmds.samples.sorted[c("G114-B", "G114-remake-C", "G114-remake-D"),], col=colors[41], pch=21) # check these out!
points(SRM.nmds.samples.sorted[c("G116-A", "G116-B"),], col=colors[42], pch=21)
points(SRM.nmds.samples.sorted[c("G120-A", "G120-B"),], col=colors[43], pch=21)
points(SRM.nmds.samples.sorted[c("G122-A", "G122-B"),], col=colors[44], pch=21)
points(SRM.nmds.samples.sorted[c("G127-A", "G127-B", "G127-C"),], col=colors[45], pch=15) #one is off
points(SRM.nmds.samples.sorted[c("G128-A", "G128-C", "G128-D"),], col=colors[46], pch=15)
points(SRM.nmds.samples.sorted[c("G129-A", "G129-B"),], col=colors[47], pch=15)
points(SRM.nmds.samples.sorted[c("G132-A", "G132-C", "G132-D"),], col=colors[48], pch=15)

#### NEXT, REMOVE SAMPLES THAT DON'T LOOK GOOD, AVERAGE TECH REPS, THEN RE-PLOT BY SITE/TREATMENT #### 

# average sample technical reps.  (there's probably an easier way to do this to not manually enter the tech rep names for each sample, possibly via a loop?); remove reps that were poor quality as per NMDS

G013 <- ave(SRM.data.normalized$'G013-A', SRM.data.normalized$'G013-C')
G120 <- ave(SRM.data.normalized$`G120-A`, SRM.data.normalized$`G120-B`)
G047 <- ave(SRM.data.normalized$`G047-A`, SRM.data.normalized$`G047-B`)
G017 <- ave(SRM.data.normalized$`G017-A`, SRM.data.normalized$`G017-B`)
G079 <- ave(SRM.data.normalized$`G079-A`, SRM.data.normalized$`G079-B`)
G127 <- ave(SRM.data.normalized$`G127-A`, SRM.data.normalized$`G127-C`) #B removed
G060 <- ave(SRM.data.normalized$`G060-A`, SRM.data.normalized$`G060-B`)
G009 <- ave(SRM.data.normalized$`G009-A`, SRM.data.normalized$`G009-B`)
G002 <- ave(SRM.data.normalized$`G002-B`, SRM.data.normalized$`G002-C`) #A removed
G128 <- ave(SRM.data.normalized$`G128-A`, SRM.data.normalized$`G128-C`,SRM.data.normalized$`G128-D`)
G016 <- ave(SRM.data.normalized$`G016-A`, SRM.data.normalized$`G016-B`, SRM.data.normalized$`G016-C`)
G071.A <- ave(SRM.data.normalized$`G071-A-A`, SRM.data.normalized$`G071-A-B`)
G114 <- ave(SRM.data.normalized$`G114-A`, SRM.data.normalized$`G114-B`, SRM.data.normalized$`G114-remake-C`, SRM.data.normalized$`G114-remake-D`)
G045 <- ave(SRM.data.normalized$`G045-A`, SRM.data.normalized$`G045-B`)
G132 <- ave(SRM.data.normalized$`G132-A`, SRM.data.normalized$`G132-C`, SRM.data.normalized$`G132-D`)
G031 <- ave(SRM.data.normalized$`G031-A`, SRM.data.normalized$`G031-B`, SRM.data.normalized$`G031-C`)
G012 <- ave(SRM.data.normalized$`G012-A`, SRM.data.normalized$`G012-B`, SRM.data.normalized$`G012-C`)
G116 <- ave(SRM.data.normalized$`G116-A`, SRM.data.normalized$`G116-B`)
G043 <- ave(SRM.data.normalized$`G043-A`, SRM.data.normalized$`G043-B`)
G015 <- ave(SRM.data.normalized$`G015-A`, SRM.data.normalized$`G015-B`)
G040 <- ave(SRM.data.normalized$`G040-A`, SRM.data.normalized$`G040-B`)
G110 <- ave(SRM.data.normalized$`G110-A`, SRM.data.normalized$`G110-B`)
G008 <- ave(SRM.data.normalized$`G008-A`, SRM.data.normalized$`G008-B`)
G109 <- ave(SRM.data.normalized$`G109-A`, SRM.data.normalized$`G109-B`, SRM.data.normalized$`G109-C`)
G122 <- ave(SRM.data.normalized$`G122-A`, SRM.data.normalized$`G122-B`)
G041 <- ave(SRM.data.normalized$`G041-A`, SRM.data.normalized$`G041-B`)
G066 <- ave(SRM.data.normalized$`G066-A`, SRM.data.normalized$`G066-B`)
G105 <- ave(SRM.data.normalized$`G105-A`, SRM.data.normalized$`G105-B`)
G032 <- ave(SRM.data.normalized$`G032-A`, SRM.data.normalized$`G032-B`)
G129 <- ave(SRM.data.normalized$`G129-A`, SRM.data.normalized$`G129-B`)
G054 <- ave(SRM.data.normalized$`G054-A`, SRM.data.normalized$`G054-B`)
G081 <- ave(SRM.data.normalized$`G081-A`, SRM.data.normalized$`G081-B`)
G003 <- ave(SRM.data.normalized$`G003-A`, SRM.data.normalized$`G003-B`) #C removed
G074 <- ave(SRM.data.normalized$`G074-A`, SRM.data.normalized$`G074-B`)
G014 <- ave(SRM.data.normalized$`G014-A`, SRM.data.normalized$`G014-B`)
G049 <- ave(SRM.data.normalized$`G049-A`, SRM.data.normalized$`G049-B`)
G053 <- ave(SRM.data.normalized$`G053-A`, SRM.data.normalized$`G053-remake-C`, SRM.data.normalized$`G053-remake-D`) #B removed 
G104 <- ave(SRM.data.normalized$`G104-A`, SRM.data.normalized$`G104-remake-C`, SRM.data.normalized$`G104-remake-D`) #B removed
G055 <- ave(SRM.data.normalized$`G055-A`, SRM.data.normalized$`G055-B`, SRM.data.normalized$`G055-C`)
G042 <- ave(SRM.data.normalized$`G042-A`, SRM.data.normalized$`G042-B`) #C removed
G064 <- ave(SRM.data.normalized$`G064-A`, SRM.data.normalized$`G064-B`)
G073 <- ave(SRM.data.normalized$`G073-A`, SRM.data.normalized$`G073-B`, SRM.data.normalized$`G073-C`)
G057 <- ave(SRM.data.normalized$`G057-A`, SRM.data.normalized$`G057-B`, SRM.data.normalized$`G057-C`)
G007 <- ave(SRM.data.normalized$`G007-A`, SRM.data.normalized$`G007-B`)
G070 <- ave(SRM.data.normalized$`G070-A`, SRM.data.normalized$`G070-B`, SRM.data.normalized$`G070-C`)
G001 <- ave(SRM.data.normalized$`G001-A`, SRM.data.normalized$`G001-B`)
G071.B <- ave(SRM.data.normalized$`G071-B-A`, SRM.data.normalized$`G071-B-B`)
G062 <- ave(SRM.data.normalized$`G062-B`, SRM.data.normalized$`G062-C`)

SRM.data.mean <- cbind.data.frame(rownames(SRM.data.normalized), G013, G120, G047, G017, G079, G127, G060, G009, G002, G128, G016, G071.A, G114, G045, G132, G031, G012, G116, G043, G015, G040, G110, G008, G109, G122, G041, G066, G105, G032, G129, G054, G081, G003, G074, G014, G049, G053, G104, G055, G042, G064, G073, G057, G007, G070,  G001, G071.B, G062) # combine all tech. replicate mean vectors into new data frame 
head(SRM.data.mean)

#### CREATE NMDS PLOT, MEAN OF TECH REPS - NOT LOG TRANSFORMED ########

#Load the source file for the biostats package, biostats.R script must be saved in working directory

source("References/biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
library(vegan)

#Transpose the file so that rows and columns are switched 
rownames(SRM.data.mean) <- SRM.data.mean[,1]
SRM.data.mean
SRM.data.mean <- SRM.data.mean[,-1]
SRM.data.mean.t <- t(SRM.data.mean) # t() function transposes
SRM.data.mean.t

#Replace NA cells with 0; metaMDS() does not handle NA's
SRM.data.mean.t.noNA <- SRM.data.mean.t
SRM.data.mean.t.noNA[is.na(SRM.data.mean.t.noNA)] <- 0
SRM.data.mean.t.noNA

#Make MDS dissimilarity matrix
#
SRM.mean.nmds <- metaMDS(SRM.data.mean.t.noNA, distance = 'bray', k = 2, trymax = 3000, autotransform = FALSE)
stressplot(SRM.mean.nmds) 
plot(SRM.mean.nmds)
# site (sample) in black circle
# species (variable) in red ticks

# make figure with sample annotations https://stat.ethz.ch/pipermail/r-sig-ecology/2011-September/002371.html
SRM.nmds.mean.samples <- scores(SRM.mean.nmds, display = "sites")

### Let's plot using ordiplot()

library(RColorBrewer)
marker = c(color = brewer.pal(4, "Set1"))
# png("Analyses/2017-September_SRM-results/2017-09-04_PRTCNorm-SRM-NMDS-plot.png")
ordiplot(SRM.mean.nmds, type="n")
points(SRM.nmds.mean.samples[c(CI.B.samples),], col=marker[1], pch=8)
points(SRM.nmds.mean.samples[c(CI.E.samples),], col=marker[1], pch=15)
points(SRM.nmds.mean.samples[c(PG.B.samples),], col=marker[2], pch=8)
points(SRM.nmds.mean.samples[c(PG.E.samples),], col=marker[2], pch=15)
points(SRM.nmds.mean.samples[c(WB.B.samples),], col=marker[3], pch=8)
points(SRM.nmds.mean.samples[c(WB.E.samples),], col=marker[3], pch=15)
points(SRM.nmds.mean.samples[c(FB.B.samples),], col=marker[4], pch=8)
points(SRM.nmds.mean.samples[c(FB.E.samples),], col=marker[4], pch=15)
legend(1.5,0.7, pch=c(rep(16,4), 8, 15), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Port Gamble", "Bare", "Eelgrass"), col=c(marker[1], marker[2], marker[3], marker[4], "black", "black"))
# dev.off()

#### Create plot with forced aspect ratio to zoom in ### 

# png("Analyses/2017-September_SRM-results/2017-09-04_PRTCNorm-SRM-NMDS-plot-zoomed.png")
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-1,2.5), ylim=c(-0.4,0.2), asp=NA)
points(SRM.nmds.mean.samples[c(CI.B.samples),], col=marker[1], pch=8)
points(SRM.nmds.mean.samples[c(CI.E.samples),], col=marker[1], pch=15)
points(SRM.nmds.mean.samples[c(PG.B.samples),], col=marker[2], pch=8)
points(SRM.nmds.mean.samples[c(PG.E.samples),], col=marker[2], pch=15)
points(SRM.nmds.mean.samples[c(WB.B.samples),], col=marker[3], pch=8)
points(SRM.nmds.mean.samples[c(WB.E.samples),], col=marker[3], pch=15)
points(SRM.nmds.mean.samples[c(FB.B.samples),], col=marker[4], pch=8)
points(SRM.nmds.mean.samples[c(FB.E.samples),], col=marker[4], pch=15)
legend(1.5,0.2, pch=c(rep(16,4), 8, 15), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Port Gamble", "Bare", "Eelgrass"), col=c(marker[1], marker[2], marker[3], marker[4], "black", "black"))
# dev.off()

#### CREATE NMDS PLOT, MEAN OF TECH REPS - LOG TRANSFORMED ########

#Transpose the file so that rows and columns are switched and normalized by log(x+1)
SRM.data.mean.t.log <- SRM.data.mean.t
SRM.data.mean.t.log[is.na(SRM.data.mean.t.log)] <- 0
SRM.data.mean.t.log <- (SRM.data.mean.t.log+1)
SRM.data.mean.t.log <- data.trans(SRM.data.mean.t.log, method = 'log', plot = FALSE)

#Make MDS dissimilarity matrix
SRM.mean.log.nmds <- metaMDS(SRM.data.mean.t.log, distance = 'bray', k = 2, trymax = 3000, autotransform = FALSE)
stressplot(SRM.mean.log.nmds) 
plot(SRM.mean.log.nmds)
# site (sample) in black circle
# species (variable) in red ticks

# make figure with sample annotations https://stat.ethz.ch/pipermail/r-sig-ecology/2011-September/002371.html
SRM.nmds.mean.log.samples <- scores(SRM.mean.log.nmds, display = "sites")
SRM.nmds.mean.log.transitions <- scores(SRM.mean.log.nmds, display = "species")
# this probably isn't necessary

### Let's plot using ordiplot()

library(RColorBrewer)
marker = c(color = brewer.pal(4, "Set1"))

# png("Analyses/2017-September_SRM-results/2017-09-04_PRTCNorm-SRM-NMDS-log-plot.png")
ordiplot(SRM.mean.log.nmds, type="n")
points(SRM.nmds.mean.log.samples[c(CI.B.samples),], col=marker[1], pch=8)
points(SRM.nmds.mean.log.samples[c(CI.E.samples),], col=marker[1], pch=15)
points(SRM.nmds.mean.log.samples[c(PG.B.samples),], col=marker[2], pch=8)
points(SRM.nmds.mean.log.samples[c(PG.E.samples),], col=marker[2], pch=15)
points(SRM.nmds.mean.log.samples[c(WB.B.samples),], col=marker[3], pch=8)
points(SRM.nmds.mean.log.samples[c(WB.E.samples),], col=marker[3], pch=15)
points(SRM.nmds.mean.log.samples[c(FB.B.samples),], col=marker[4], pch=8)
points(SRM.nmds.mean.log.samples[c(FB.E.samples),], col=marker[4], pch=15)
legend(.01,0.01, pch=c(rep(16,4), 8, 15), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Port Gamble", "Bare", "Eelgrass"), col=c(marker[1], marker[2], marker[3], marker[4], "black", "black"))
# dev.off()

#### Create plot with forced aspect ratio to zoom in ### 

# png("Analyses/2017-September_SRM-results/2017-09-04_PRTCNorm-SRM-NMDS-plot-log-zoomed.png")
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-.005,.01), ylim=c(-0.002,0.002), asp=NA)
points(SRM.nmds.mean.log.samples[c(CI.B.samples),], col=marker[1], pch=8)
points(SRM.nmds.mean.log.samples[c(CI.E.samples),], col=marker[1], pch=15)
points(SRM.nmds.mean.log.samples[c(PG.B.samples),], col=marker[2], pch=8)
points(SRM.nmds.mean.log.samples[c(PG.E.samples),], col=marker[2], pch=15)
points(SRM.nmds.mean.log.samples[c(WB.B.samples),], col=marker[3], pch=8)
points(SRM.nmds.mean.log.samples[c(WB.E.samples),], col=marker[3], pch=15)
points(SRM.nmds.mean.log.samples[c(FB.B.samples),], col=marker[4], pch=8)
points(SRM.nmds.mean.log.samples[c(FB.E.samples),], col=marker[4], pch=15)
legend(1.5,0.2, pch=c(rep(16,4), 8, 15), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Port Gamble", "Bare", "Eelgrass"), col=c(marker[1], marker[2], marker[3], marker[4], "black", "black"))
# dev.off()

#### MORE STRINGENT REMOVAL OF POOR QUALITY POINTS (VERY, VERY STRINGENT)
# average sample technical reps.  (there's probably an easier way to do this to not manually enter the tech rep names for each sample, possibly via a loop?); remove reps that were poor quality as per NMDS

G013.s <- ave(SRM.data.normalized$'G013-A', SRM.data.normalized$'G013-C')
G120.s <- ave(SRM.data.normalized$`G120-A`, SRM.data.normalized$`G120-B`)
G047.s <- ave(SRM.data.normalized$`G047-A`, SRM.data.normalized$`G047-B`)
G017.s <- ave(SRM.data.normalized$`G017-A`, SRM.data.normalized$`G017-B`)
G079.s <- ave(SRM.data.normalized$`G079-A`, SRM.data.normalized$`G079-B`)
G127.s <- ave(SRM.data.normalized$`G127-A`, SRM.data.normalized$`G127-C`) #B removed
G060.s <- ave(SRM.data.normalized$`G060-A`, SRM.data.normalized$`G060-B`)
G009.s <- ave(SRM.data.normalized$`G009-A`, SRM.data.normalized$`G009-B`)
G002.s <- ave(SRM.data.normalized$`G002-B`, SRM.data.normalized$`G002-C`) #A removed
G128.s <- ave(SRM.data.normalized$`G128-C`,SRM.data.normalized$`G128-D`)
G114.s <- ave(SRM.data.normalized$`G114-A`, SRM.data.normalized$`G114-B`, SRM.data.normalized$`G114-remake-C`, SRM.data.normalized$`G114-remake-D`)
G045.s <- ave(SRM.data.normalized$`G045-A`, SRM.data.normalized$`G045-B`)
G132.s <- ave(SRM.data.normalized$`G132-A`, SRM.data.normalized$`G132-C`, SRM.data.normalized$`G132-D`)
G031.s <- ave(SRM.data.normalized$`G031-A`, SRM.data.normalized$`G031-B`, SRM.data.normalized$`G031-C`)
G012.s <- ave(SRM.data.normalized$`G012-A`, SRM.data.normalized$`G012-B`, SRM.data.normalized$`G012-C`)
G015.s <- ave(SRM.data.normalized$`G015-A`, SRM.data.normalized$`G015-B`)
G040.s <- ave(SRM.data.normalized$`G040-A`, SRM.data.normalized$`G040-B`)
G110.s <- ave(SRM.data.normalized$`G110-A`, SRM.data.normalized$`G110-B`)
G008.s <- ave(SRM.data.normalized$`G008-A`, SRM.data.normalized$`G008-B`)
G122.s <- ave(SRM.data.normalized$`G122-A`, SRM.data.normalized$`G122-B`)
G066.s <- ave(SRM.data.normalized$`G066-A`, SRM.data.normalized$`G066-B`)
G105.s <- ave(SRM.data.normalized$`G105-A`, SRM.data.normalized$`G105-B`)
G032.s <- ave(SRM.data.normalized$`G032-A`, SRM.data.normalized$`G032-B`)
G054.s <- ave(SRM.data.normalized$`G054-A`, SRM.data.normalized$`G054-B`)
G003.s <- ave(SRM.data.normalized$`G003-A`, SRM.data.normalized$`G003-B`) #C removed
G074.s <- ave(SRM.data.normalized$`G074-A`, SRM.data.normalized$`G074-B`)
G014.s <- ave(SRM.data.normalized$`G014-A`, SRM.data.normalized$`G014-B`)
G049.s <- ave(SRM.data.normalized$`G049-A`, SRM.data.normalized$`G049-B`)
G053.s <- ave(SRM.data.normalized$`G053-A`, SRM.data.normalized$`G053-remake-C`, SRM.data.normalized$`G053-remake-D`) #B removed 
G104.s <- ave(SRM.data.normalized$`G104-A`, SRM.data.normalized$`G104-remake-C`, SRM.data.normalized$`G104-remake-D`) #B removed
G007.s <- ave(SRM.data.normalized$`G007-A`, SRM.data.normalized$`G007-B`)
G070.s <- ave(SRM.data.normalized$`G070-A`, SRM.data.normalized$`G070-B`, SRM.data.normalized$`G070-C`)
G071.B.s <- ave(SRM.data.normalized$`G071-B-A`, SRM.data.normalized$`G071-B-B`)
G062.s <- ave(SRM.data.normalized$`G062-B`, SRM.data.normalized$`G062-C`)

SRM.data.mean.stringent <- cbind.data.frame(rownames(SRM.data.normalized), G013.s, G120.s, G047.s, G017.s, G079.s, G127.s, G060.s, G009.s, G002.s, G128.s, G114.s, G045.s, G132.s, G031.s, G040.s, G110.s, G122.s, G066.s, G105.s, G032.s, G054.s, G074.s, G014.s, G049.s, G053.s, G104.s, G007.s, G070.s, G071.B.s, G062.s) # combine all tech. replicate mean vectors into new data frame 

#### CREATE NMDS PLOT, MEAN OF TECH REPS, STRINGENT SELECTION - NOT LOG TRANSFORMED ########

#Transpose the file so that rows and columns are switched 
rownames(SRM.data.mean.stringent) <- SRM.data.mean.stringent[,1]
SRM.data.mean.stringent <- SRM.data.mean.stringent[,-1]
SRM.data.mean.stringent.t <- t(SRM.data.mean.stringent) # t() function transposes
rownames(SRM.data.mean.stringent.t)

#Replace NA cells with 0; metaMDS() does not handle NA's
SRM.data.mean.s.t.noNA <- SRM.data.mean.stringent.t
SRM.data.mean.s.t.noNA[is.na(SRM.data.mean.s.t.noNA)] <- 0
SRM.data.mean.s.t.noNA

#Make MDS dissimilarity matrix

SRM.mean.s.nmds <- metaMDS(SRM.data.mean.s.t.noNA, distance = 'bray', k = 2, trymax = 3000, autotransform = FALSE)
stressplot(SRM.mean.s.nmds) 
plot(SRM.mean.s.nmds)

# Eigenvectors
eigen.s <- envfit(SRM.mean.s.nmds$points, SRM.data.mean.s.t.noNA, perm=1000)
eigen.s

# make figure with sample annotations https://stat.ethz.ch/pipermail/r-sig-ecology/2011-September/002371.html
SRM.nmds.mean.s.samples <- scores(SRM.mean.s.nmds, display = "sites")

### Let's plot using ordiplot()

library(RColorBrewer)
marker = c(color = brewer.pal(4, "Set1"))

CI.B.samples.s <- c("G013.s", "G014.s", "G017.s")
CI.E.samples.s <- c("G002.s", "G007.s", "G009.s")
PG.B.samples.s <- c("G062.s", "G066.s", "G070.s", "G079.s")
PG.E.samples.s <- c("G031.s", "G032.s", "G074.s", "G071.B.s")
WB.B.samples.s <- c("G104.s", "G105.s", "G114.s", "G120.s", "G122.s")
WB.E.samples.s <- c("G110.s", "G127.s", "G128.s", "G132.s")
FB.B.samples.s <- c("G040.s", "G060.s")
FB.E.samples.s <- c("G045.s", "G047.s", "G049.s", "G053.s", "G054.s")

# png("Analyses/2017-September_SRM-results/2017-09-04_PRTCNorm-SRM-NMDS-stringent-plot.png")
ordiplot(SRM.mean.s.nmds, type="n")
points(SRM.nmds.mean.s.samples[c(CI.B.samples.s),], col=marker[1], pch=8)
points(SRM.nmds.mean.s.samples[c(CI.E.samples.s),], col=marker[1], pch=15)
points(SRM.nmds.mean.s.samples[c(PG.B.samples.s),], col=marker[2], pch=8)
points(SRM.nmds.mean.s.samples[c(PG.E.samples.s),], col=marker[2], pch=15)
points(SRM.nmds.mean.s.samples[c(WB.B.samples.s),], col=marker[3], pch=8)
points(SRM.nmds.mean.s.samples[c(WB.E.samples.s),], col=marker[3], pch=15)
points(SRM.nmds.mean.s.samples[c(FB.B.samples.s),], col=marker[4], pch=8)
points(SRM.nmds.mean.s.samples[c(FB.E.samples.s),], col=marker[4], pch=15)
# dev.off()

legend(1.5,0.7, pch=c(rep(16,4), 8, 15), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Port Gamble", "Bare", "Eelgrass"), col=c(marker[1], marker[2], marker[3], marker[4], "black", "black"))

#### Create plot with forced aspect ratio to zoom in ### 

# png("Analyses/2017-September_SRM-results/2017-09-04_PRTCNorm-SRM-NMDS-stringent-plot-zoomed.png")
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-1,1), ylim=c(-0.2,0.2), asp=NA)
points(SRM.nmds.mean.s.samples[c(CI.B.samples.s),], col=marker[1], pch=8)
points(SRM.nmds.mean.s.samples[c(CI.E.samples.s),], col=marker[1], pch=15)
points(SRM.nmds.mean.s.samples[c(PG.B.samples.s),], col=marker[2], pch=8)
points(SRM.nmds.mean.s.samples[c(PG.E.samples.s),], col=marker[2], pch=15)
points(SRM.nmds.mean.s.samples[c(WB.B.samples.s),], col=marker[3], pch=8)
points(SRM.nmds.mean.s.samples[c(WB.E.samples.s),], col=marker[3], pch=15)
points(SRM.nmds.mean.s.samples[c(FB.B.samples.s),], col=marker[4], pch=8)
points(SRM.nmds.mean.s.samples[c(FB.E.samples.s),], col=marker[4], pch=15)
legend(0.5,-0.06, pch=c(rep(16,4), 8, 15), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Port Gamble", "Bare", "Eelgrass"), col=c(marker[1], marker[2], marker[3], marker[4], "black", "black"))
# dev.off()

#### CREATE NMDS PLOT, MEAN OF TECH REPS, STRINGENT SELECTION -  LOG TRANSFORMED ########

#Transform by log(x+1)
SRM.data.mean.stringent.t.log <- SRM.data.mean.stringent.t
SRM.data.mean.stringent.t.log[is.na(SRM.data.mean.stringent.t.log)] <- 0 #change NA's to zero
SRM.data.mean.stringent.t.log <- (SRM.data.mean.stringent.t.log+1)
SRM.data.mean.stringent.t.log <- data.trans(SRM.data.mean.stringent.t.log, method = 'log', plot = FALSE)

#Make MDS dissimilarity matrix
SRM.mean.s.log.nmds <- metaMDS(SRM.data.mean.stringent.t.log, distance = 'bray', k = 2, trymax = 3000, autotransform = FALSE)
stressplot(SRM.mean.s.log.nmds) 
plot(SRM.mean.s.log.nmds)

#eigenvector
eigen.s.log <- envfit(SRM.mean.s.log.nmds$points, SRM.data.mean.s.t.noNA, perm=1000)
eigen.s.log

# make figure with sample annotations https://stat.ethz.ch/pipermail/r-sig-ecology/2011-September/002371.html
SRM.nmds.mean.s.log.samples <- scores(SRM.mean.s.log.nmds, display = "sites")
SRM.nmds.mean.s.log.transitions <- scores(SRM.mean.s.log.nmds, display = "species")
# this probably isn't necessary

### Let's plot using ordiplot()

# png("Analyses/2017-September_SRM-results/2017-09-04_PRTCNorm-SRM-NMDS-stringent-plot-log.png")
ordiplot(SRM.mean.s.log.nmds, type="n")
points(SRM.nmds.mean.s.log.samples[c(CI.B.samples.s),], col=marker[1], pch=8)
points(SRM.nmds.mean.s.log.samples[c(CI.E.samples.s),], col=marker[1], pch=15)
points(SRM.nmds.mean.s.log.samples[c(PG.B.samples.s),], col=marker[2], pch=8)
points(SRM.nmds.mean.s.log.samples[c(PG.E.samples.s),], col=marker[2], pch=15)
points(SRM.nmds.mean.s.log.samples[c(WB.B.samples.s),], col=marker[3], pch=8)
points(SRM.nmds.mean.s.log.samples[c(WB.E.samples.s),], col=marker[3], pch=15)
points(SRM.nmds.mean.s.log.samples[c(FB.B.samples.s),], col=marker[4], pch=8)
points(SRM.nmds.mean.s.log.samples[c(FB.E.samples.s),], col=marker[4], pch=15)
legend(0.01,0.01, pch=c(rep(16,4), 8, 15), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Port Gamble", "Bare", "Eelgrass"), col=c(marker[1], marker[2], marker[3], marker[4], "black", "black"))
# dev.off()

#### Create plot with forced aspect ratio to zoom in ### 

# png("Analyses/2017-September_SRM-results/2017-09-04_PRTCNorm-SRM-NMDS-stringent-plot-log-zoomed.png")
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-.005,.010), ylim=c(-0.0025,0.0025), asp=NA)
points(SRM.nmds.mean.s.log.samples[c(CI.B.samples.s),], col=marker[1], pch=8)
points(SRM.nmds.mean.s.log.samples[c(CI.E.samples.s),], col=marker[1], pch=15)
points(SRM.nmds.mean.s.log.samples[c(PG.B.samples.s),], col=marker[2], pch=8)
points(SRM.nmds.mean.s.log.samples[c(PG.E.samples.s),], col=marker[2], pch=15)
points(SRM.nmds.mean.s.log.samples[c(WB.B.samples.s),], col=marker[3], pch=8)
points(SRM.nmds.mean.s.log.samples[c(WB.E.samples.s),], col=marker[3], pch=15)
points(SRM.nmds.mean.s.log.samples[c(FB.B.samples.s),], col=marker[4], pch=8)
points(SRM.nmds.mean.s.log.samples[c(FB.E.samples.s),], col=marker[4], pch=15)
legend(0.006, 0.0023, pch=c(rep(16,4), 8, 15), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Port Gamble", "Bare", "Eelgrass"), col=c(marker[1], marker[2], marker[3], marker[4], "black", "black"))
# dev.off()