### Note: this script works with datasets downloaded from Skyline.  Skyline data is exported as report, and this script works with data with the following metrics: Protein Name, Transitions, Peptide Sequence, Fragment Ion, Peptide Retention Time, Area

# Script #1 in data processing for NOT-NORMALIZED DATA

############# IMPORT DATASETS ########################################################################################

setwd("~/Documents/Roberts Lab/Geoduck-DNR/") 
SRMreport <- read.csv("Data/2017-08-11_Transition Results_LHS modified-noRT-pivoted.csv", header=FALSE, na.strings = "#N/A", stringsAsFactors = FALSE) # import local file
SRMsequence <- read.csv("Data/2017-07-28_SRM-Sequence-final.csv", header=TRUE, stringsAsFactors = FALSE)
sample.key <- read.csv("Data/2017-08-14-Geoduck-samples.csv", header=TRUE, stringsAsFactors = FALSE)
dilution.curve <- read.csv("Data/2017-09-05_Dilution-Curve-Results.csv", header=TRUE, stringsAsFactors = FALSE)
OutplantData <- read.csv("Data/Outplant-Temp-Data.csv", header=TRUE, stringsAsFactors =TRUE)
SRMsamples <- noquote(as.character(c("G013", "G120", "G047", "G017", "G079", "G127", "G060", "G009", "G002", "G128", "G016", "G071-A", "G114", "G045", "G132", "G031", "G012", "G116", "G043", "G015", "G040", "G110", "G008", "G109", "G122", "G041", "G066", "G105", "G032", "G129", "G054", "G081", "G003", "G074", "G014", "G049", "G053", "G104", "G055", "G042", "G064", "G073", "G057", "G007", "G070", "G001", "G071-B", "G062")))

############ REPLACE REP NAMES WITH SAMPLE NAMES ###################################################################
# I could also probably use the function merge() for this task

rep.names <- SRMreport[1,] # create vector of replicate names
rep.names.short <- noquote(gsub(' Area', '', rep.names)) # remove Area from rep name, and don't include quotes 
rep.names.short <- noquote(gsub('2017_July_10_bivalves_', '', rep.names.short)) #remove the extra long rep name that is a residual from the .raw file name
repsTOsamples <- as.data.frame(SRMsequence[,c(2,3,5)])
library(dplyr)
repsTOsamples.filtered <- filter(repsTOsamples, repsTOsamples[,1] %in% rep.names.short)
samples <- as.character(repsTOsamples.filtered$Sample...rep.name)
other.headers <- as.character(rep.names.short[1:4])
samples.vector <- noquote(c(other.headers, samples, stringsAsFactors = FALSE))
samples.vector <- samples.vector[-121]
SRM.data <- SRMreport
SRM.data[1,] <- samples.vector
colnames(SRM.data) <- SRM.data[1,] #make first row column names

############ ANNOTATE SAMPLE NAMES WITH SITE & TREATMENT ##################################################################################

repsTOsamples.filtered.annotated <- filter(sample.key[,c(8,9)], sample.key$PRVial %in% repsTOsamples.filtered$Comment) #pull site & treatment from sample key
length(SRMsamples) == nrow(repsTOsamples.filtered.annotated) # should equal TRUE if I haven't lost any sample data
repsTOsamples.filtered.annotated # ISSUE IDENTIFIED: missing 71-A & 71-B, need to fix 

# Add G071-A & G071-B coding to the annotated key
s71.A <- data.frame(matrix(0, ncol=2, nrow=1))
s71.A[1,1] <- "G071.A"
s71.A[1,2] <- "PG-E"
colnames(s71.A) <- colnames(repsTOsamples.filtered.annotated)
s71.B <- data.frame(matrix(0, ncol=2, nrow=1))
s71.B[1,1] <- "G071.B"
s71.B[1,2] <- "PG-E"
colnames(s71.B) <- colnames(repsTOsamples.filtered.annotated)
sample.key.annotated <- rbind(repsTOsamples.filtered.annotated, s71.A, s71.B) # row bind annotated key w/ 71 info

# Subset sample names for site & treatment combos
CI.E <- sample.key.annotated[c(sample.key.annotated$Sample.Shorthand == "CI-E"),]
CI.B <- sample.key.annotated[c(sample.key.annotated$Sample.Shorthand == "CI-B"),]
PG.E <- sample.key.annotated[c(sample.key.annotated$Sample.Shorthand == "PG-E"),]
PG.B <- sample.key.annotated[c(sample.key.annotated$Sample.Shorthand == "PG-B"),]
WB.E <- sample.key.annotated[c(sample.key.annotated$Sample.Shorthand == "WB-E"),]
WB.B <- sample.key.annotated[c(sample.key.annotated$Sample.Shorthand == "WB-B"),]
FB.E <- sample.key.annotated[c(sample.key.annotated$Sample.Shorthand == "FB-E"),]
FB.B <- sample.key.annotated[c(sample.key.annotated$Sample.Shorthand == "FB-B"),]

# Isolate just sample names for each site/treatment combo
CI.E.samples <- CI.E$PRVial
CI.B.samples <- CI.B$PRVial
PG.E.samples <- PG.E$PRVial
PG.B.samples <- PG.B$PRVial
WB.E.samples <- WB.E$PRVial
WB.B.samples <- WB.B$PRVial
FB.E.samples <- FB.E$PRVial
FB.B.samples <- FB.B$PRVial

# Isolate eelgrass and bare group sample names
Eelgrass.samples <- c(CI.E.samples, PG.E.samples, WB.E.samples, FB.E.samples)
Bare.samples <- c(CI.B.samples, PG.B.samples, WB.B.samples, FB.B.samples)

############ CONVERT AREA DATA TO NUMERIC FORMAT #######################################################################################

SRM.data.numeric <- SRM.data[-1,] # First, remove row #1 (excess info) & change all area values to numeric, so I can average, etc. I know that my area data is from column 5 to 120 
SRM.data.numeric[,5:120] <- as.numeric( 
  as.character(
    unlist(
      SRM.data.numeric[,5:120])
  )
)
is.numeric(SRM.data.numeric[5,20]) # confirm area data is numeric, using a random cell. Should equal TRUE.

############ NAME EACH ROW WITH A UNIQUE TRANSITION ID ##############################################################
nTransitions <- length(SRM.data.numeric$Transition) # How many transitions are there
Transition.ID <- vector(length=nTransitions) # create empty vector with length= number of transitions
for (i in 1:nTransitions) {  
  Transition.ID[i] <- paste(SRM.data.numeric[i,3], SRM.data.numeric[i,4])}  # loop that fills empty vector with unique transition ID, built from the peptide sequence (column 3) and the fragment ion (columm 4)
Transition.ID # confirm correctly named transition IDs
length(SRM.data.numeric$Transition) == length(Transition.ID) # confirm that I didn't lose any transitions
row.names(SRM.data.numeric) <- Transition.ID # assign newly created transition IDs as row names
write.csv(SRM.data.numeric, file="Analyses/2017-September_SRM-results/2017-09-04_SRM-data-NotNORM-annotated.csv") #write this file out for safe keeping

########### REMOVE POOR QUALITY PEPTIDES IDENTIFIED VIA SKYLINE & DILUTION CURVE RESULTS #############################
# Poor quality, determined via Skyline due to lack of consistent signal as compared to other peptides in the protein: 
  # Superoxide Dismutase: THGAPTDEER
  # PDI: NNKPSDYQGGR
  # Na/K: MVTGDNVNTAR
# Poor quality, determined via Dilution curve
  # HSP70: TTPSYVAFNDTER
  # Peroxiredoxin: LVQAFQFTDK
  # Ras-related Rab: QITMNDLPVGR & VVLVGDSGVGK
  # Na/K: AQLWDTAGQER & MVTGDNVNTAR

nrow(SRM.data.screened.noPRTC) / nrow(SRM.data[2:116,])


SRM.data.screened <- SRM.data.numeric[!grepl(c("THGAPTDEER|NNKPSDYQGGR|MVTGDNVNTAR|TTPSYVAFNDTER|LVQAFQFTDK|QITMNDLPVGR|VVLVGDSGVGK|AQLWDTAGQER"), SRM.data.numeric$`Peptide Sequence`),]
SRM.data.screened.noPRTC <- SRM.data.screened[!grepl("PRTC peptides", SRM.data.screened$`Protein Name`),]
write.csv(SRM.data.screened.noPRTC, file="Analyses/2017-September_SRM-results/2017-09-04_SRM-data-NotNORM-screenednoPRTC.csv")

############ CREATE NMDS PLOT ########################################################################################

#Load the source file for the biostats package, biostats.R script must be saved in working directory

source("References/biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
library(vegan)

#Transpose the file so that rows and columns are switched 
SRM.data.t <- t(SRM.data.screened.noPRTC[, -1:-4]) # t() function transposes, removes PRTC transitions, extraneous info

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

###### Make figure
# site (aka sample) in black circle
# species (aka transition) in red ticks
plot(SRM.nmds)

# Principal Component Analysis
SRM.nmds.pca <- rda(SRM.data.t.noNA, scale = TRUE)
summary(SRM.nmds.pca)
plot(SRM.nmds.pca, scaling = 3)
dim(SRM.data.t.noNA)
biplot(SRM.nmds.pca, scaling = -1)
SRM.nmds.ca <- cca(SRM.data.t.noNA)
plot(SRM.nmds.ca)
#inertia is the sum of all variance in transitions; eigenvalues sum to total inertia, aka each eigenvalue "explains" a certain proportion of the total variance. Percent that each eigenvalue is responsible for total variance is: eigenvalue/total inertia. For example, PC1/total inertia = 83%


# make figure with sample annotations https://stat.ethz.ch/pipermail/r-sig-ecology/2011-September/002371.html
SRM.nmds.samples <- scores(SRM.nmds, display = "sites")
SRM.nmds.samples.sorted <- SRM.nmds.samples[ order(row.names(SRM.nmds.samples)), ]
rownames(SRM.nmds.samples.sorted)
colors <- colorRampPalette(brewer.pal(8,"Dark2"))(48)

### PLOTTING ALL REPS WITH SAMPLE NUMBER ID'S ### 
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-1,3), ylim=c(-0.5,0.5), asp=NA, main= "NMDS of SRM data for technical rep QA")
text(SRM.nmds.samples.sorted[c("G001-A", "G001-B"),], labels=c("1A", "1B"), col=colors[1])
text(SRM.nmds.samples.sorted[c("G002-A", "G002-B", "G002-C"),], labels=c("2A", "2B", "2C"), col=colors[2])
text(SRM.nmds.samples.sorted[c("G003-A", "G003-B", "G003-C"),], labels=c("3A", "3B", "3C"), col=colors[3]) #G003-C is very different
text(SRM.nmds.samples.sorted[c("G007-A", "G007-B"),], labels=c("7A", "7B"), col=colors[4])
text(SRM.nmds.samples.sorted[c("G008-A", "G008-B"),], labels=c("8A", "8B"), col=colors[5])
text(SRM.nmds.samples.sorted[c("G009-A", "G009-B"),], labels=c("9A", "9B"), col=colors[6])
text(SRM.nmds.samples.sorted[c("G012-A", "G012-B", "G012-C"),], labels=c("12A", "12B", "12C"), col=colors[7])
text(SRM.nmds.samples.sorted[c("G013-A", "G013-C"),], labels=c("13A", "13B"), col=colors[8])
text(SRM.nmds.samples.sorted[c("G014-A", "G014-B"),], labels=c("14A", "14B"), col=colors[9])
text(SRM.nmds.samples.sorted[c("G015-A", "G015-B"),], labels=c("15A", "15B"), col=colors[10])
text(SRM.nmds.samples.sorted[c("G016-A", "G016-B", "G016-C"),], labels=c("16A", "16B", "16C"), col=colors[11]) # not perfect 
text(SRM.nmds.samples.sorted[c("G017-A", "G017-B"),], labels=c("17A", "17B"), col=colors[12]) # not great
text(SRM.nmds.samples.sorted[c("G031-A", "G031-B", "G031-C"),], labels=c("31A", "31B", "31C"), col=colors[13]) #31- not great
text(SRM.nmds.samples.sorted[c("G032-A", "G032-B"),], labels=c("32A", "32B"), col=colors[14])
text(SRM.nmds.samples.sorted[c("G040-A", "G040-B"),], labels=c("40A", "40B"), col=colors[15])
text(SRM.nmds.samples.sorted[c("G041-A", "G041-B"),], labels=c("41A", "41B"), col=colors[16])
text(SRM.nmds.samples.sorted[c("G042-A", "G042-B", "G042-C"),], labels=c("42A", "42B", "42C"), col=colors[17]) #G042-C very off
text(SRM.nmds.samples.sorted[c("G043-A", "G043-B"),], labels=c("43A", "43B"), col=colors[18])
text(SRM.nmds.samples.sorted[c("G045-A", "G045-B"),], labels=c("45A", "45B"), col=colors[19])
text(SRM.nmds.samples.sorted[c("G047-A", "G047-B"),], labels=c("47A", "47B"), col=colors[20])
text(SRM.nmds.samples.sorted[c("G049-A", "G049-B"),], labels=c("49A", "49B"), col=colors[21])
text(SRM.nmds.samples.sorted[c("G053-A", "G053-B", "G053-remake-C", "G053-remake-D"),], labels=c("53A", "53B", "53C-redo", "53D-redo"), col=colors[22]) #53-B bad
text(SRM.nmds.samples.sorted[c("G054-A", "G054-B"),], labels=c("54A", "54B"), col=colors[23])
text(SRM.nmds.samples.sorted[c("G055-A", "G055-B"),], labels=c("55A", "55B"), col=colors[24]) #one is very off
text(SRM.nmds.samples.sorted[c("G057-A", "G057-B", "G057-C"),], labels=c("57A", "57B", "57C"), col=colors[25]) #57-B BAD
text(SRM.nmds.samples.sorted[c("G060-A", "G060-B"),], labels=c("60A", "60B"), col=colors[26]) 
text(SRM.nmds.samples.sorted[c("G062-B", "G062-C"),], labels=c("62B", "62C"), col=colors[27])
text(SRM.nmds.samples.sorted[c("G064-A", "G064-B"),], labels=c("64A", "64B"), col=colors[28])
text(SRM.nmds.samples.sorted[c("G066-A", "G066-B"),], labels=c("66A", "66B"), col=colors[29])
text(SRM.nmds.samples.sorted[c("G070-A", "G070-B", "G070-C"),], labels=c("70A", "70B", "70C"), col=colors[30])
text(SRM.nmds.samples.sorted[c("G071-A-A", "G071-A-B"),], labels=c("71aA", "71aB"), col=colors[31])
text(SRM.nmds.samples.sorted[c("G071-B-A", "G071-B-B"),], labels=c("71bA", "71bB"), col=colors[32])
text(SRM.nmds.samples.sorted[c("G073-A", "G073-B", "G073-C"),], labels=c("73A", "73B", "73C"), col=colors[33]) #all three very off
text(SRM.nmds.samples.sorted[c("G074-A", "G074-B"),], labels=c("74A", "74B"), col=colors[34])
text(SRM.nmds.samples.sorted[c("G079-A", "G079-B"),], labels=c("79A", "79B"), col=colors[35])
text(SRM.nmds.samples.sorted[c("G081-A", "G081-B"),], labels=c("81A", "81B"), col=colors[36])
text(SRM.nmds.samples.sorted[c("G104-A", "G104-B", "G104-remake-C", "G104-remake-D"),], labels=c("104A", "104B", "104C-redo", "104D-redo"), col=colors[37]) #B
text(SRM.nmds.samples.sorted[c("G105-A", "G105-B"),], labels=c("105A", "105B"), col=colors[38])
text(SRM.nmds.samples.sorted[c("G109-A", "G109-B", "G109-C"),], labels=c("109A", "109B", "109C"), col=colors[39]) #109-B BAD
text(SRM.nmds.samples.sorted[c("G110-A", "G110-B"),], labels=c("110A", "110B"), col=colors[40])
text(SRM.nmds.samples.sorted[c("G114-B", "G114-remake-C", "G114-remake-D"),], labels=c("114B", "114C-redo", "114D-redo" ), col=colors[41]) # check these out!
text(SRM.nmds.samples.sorted[c("G116-A", "G116-B"),], labels=c("116A", "116B"), col=colors[42])
text(SRM.nmds.samples.sorted[c("G120-A", "G120-B"),], labels=c("120A", "120B"), col=colors[43])
text(SRM.nmds.samples.sorted[c("G122-A", "G122-B"),], labels=c("122A", "122B"), col=colors[44])
text(SRM.nmds.samples.sorted[c("G127-A", "G127-B", "G127-C"),], labels=c("127A", "127B"), col=colors[45]) #127 B BAD
text(SRM.nmds.samples.sorted[c("G128-A", "G128-C", "G128-D"),], labels=c("128A", "128C", "128D"), col=colors[46])
text(SRM.nmds.samples.sorted[c("G129-A", "G129-B"),], labels=c("129A", "129B"), col=colors[47])
text(SRM.nmds.samples.sorted[c("G132-A", "G132-C", "G132-D"),], labels=c("132A", "132C", "132D"), col=colors[48])


### PLOTTING ALL REPS COLOR CODED AND WITH TREATMENT SYMBOL ### 
# symbol key
# 15 = eelgrass = filled square
# 21 = bare = open circle
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-1,3), ylim=c(-0.5,0.5), asp=NA)
points(SRM.nmds.samples.sorted[c("G001-A", "G001-B"),], col=colors[1], pch=15)
points(SRM.nmds.samples.sorted[c("G002-A", "G002-B", "G002-C"),], col=colors[2], pch=15) #GOO2-B very diff
points(SRM.nmds.samples.sorted[c("G003-A", "G003-B", "G003-C"),], col=colors[3], pch=15) #G003-C is very different
points(SRM.nmds.samples.sorted[c("G007-A", "G007-B"),], col=colors[4], pch=15)
points(SRM.nmds.samples.sorted[c("G008-A", "G008-B"),], col=colors[5], pch=15)
points(SRM.nmds.samples.sorted[c("G009-A", "G009-B"),], col=colors[6], pch=15)
points(SRM.nmds.samples.sorted[c("G012-A", "G012-B", "G012-C"),], col=colors[7], pch=21)
points(SRM.nmds.samples.sorted[c("G013-A", "G013-C"),], col=colors[8], pch=21)
points(SRM.nmds.samples.sorted[c("G014-A", "G014-B"),], col=colors[9], pch=21)
points(SRM.nmds.samples.sorted[c("G015-A", "G015-B"),], col=colors[10], pch=21)
points(SRM.nmds.samples.sorted[c("G016-A", "G016-B", "G016-C"),], col=colors[11], pch=21) #GO16-A BAD 
points(SRM.nmds.samples.sorted[c("G017-A", "G017-B"),], col=colors[12], pch=21) #BAD
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
points(SRM.nmds.samples.sorted[c("G057-A", "G057-B", "G057-C"),], col=colors[25], pch=21) #57-B BAD
points(SRM.nmds.samples.sorted[c("G060-A", "G060-B"),], col=colors[26], pch=21) 
points(SRM.nmds.samples.sorted[c("G062-B", "G062-C"),], col=colors[27], pch=21)
points(SRM.nmds.samples.sorted[c("G064-A", "G064-B"),], col=colors[28], pch=21)
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
points(SRM.nmds.samples.sorted[c("G109-A", "G109-B", "G109-C"),], col=colors[39], pch=15) #109-B BAD
points(SRM.nmds.samples.sorted[c("G110-A", "G110-B"),], col=colors[40], pch=15)
points(SRM.nmds.samples.sorted[c("G114-B", "G114-remake-C", "G114-remake-D"),], col=colors[41], pch=21) # check these out!
points(SRM.nmds.samples.sorted[c("G116-A", "G116-B"),], col=colors[42], pch=21)
points(SRM.nmds.samples.sorted[c("G120-A", "G120-B"),], col=colors[43], pch=21)
points(SRM.nmds.samples.sorted[c("G122-A", "G122-B"),], col=colors[44], pch=21)
points(SRM.nmds.samples.sorted[c("G127-A", "G127-B", "G127-C"),], col=colors[45], pch=15) #127 B BAD
points(SRM.nmds.samples.sorted[c("G128-A", "G128-C", "G128-D"),], col=colors[46], pch=15)
points(SRM.nmds.samples.sorted[c("G129-A", "G129-B"),], col=colors[47], pch=15)
points(SRM.nmds.samples.sorted[c("G132-A", "G132-C", "G132-D"),], col=colors[48], pch=15)


#### NEXT, REMOVE SAMPLES THAT DON'T LOOK GOOD, AVERAGE TECH REPS, THEN RE-PLOT BY SITE/TREATMENT #### 

# average sample technical reps.  (there's probably an easier way to do this to not manually enter the tech rep names for each sample, possibly via a loop?); remove reps that were poor quality as per NMDS

# Mean of tech reps
G001 <- ave(SRM.data.screened.noPRTC$`G001-A`, SRM.data.screened.noPRTC$`G001-B`)
G002 <- ave(SRM.data.screened.noPRTC$`G002-A`, SRM.data.screened.noPRTC$`G002-B`, SRM.data.screened.noPRTC$`G002-C`)
G003 <- ave(SRM.data.screened.noPRTC$`G003-A`, SRM.data.screened.noPRTC$`G003-B`) #C removed
G007 <- ave(SRM.data.screened.noPRTC$`G007-A`, SRM.data.screened.noPRTC$`G007-B`)
G008 <- ave(SRM.data.screened.noPRTC$`G008-A`, SRM.data.screened.noPRTC$`G008-B`)
G009 <- ave(SRM.data.screened.noPRTC$`G009-A`, SRM.data.screened.noPRTC$`G009-B`)
G110 <- ave(SRM.data.screened.noPRTC$`G110-A`, SRM.data.screened.noPRTC$`G110-B`)
G012 <- ave(SRM.data.screened.noPRTC$`G012-A`, SRM.data.screened.noPRTC$`G012-B`, SRM.data.screened.noPRTC$`G012-C`)
G013 <- ave(SRM.data.screened.noPRTC$'G013-A', SRM.data.screened.noPRTC$'G013-C')
G014 <- ave(SRM.data.screened.noPRTC$`G014-A`, SRM.data.screened.noPRTC$`G014-B`)
G015 <- ave(SRM.data.screened.noPRTC$`G015-A`, SRM.data.screened.noPRTC$`G015-B`)
G016 <- ave(SRM.data.screened.noPRTC$`G016-A`, SRM.data.screened.noPRTC$`G016-B`, SRM.data.screened.noPRTC$`G016-C`)
G017 <- ave(SRM.data.screened.noPRTC$`G017-A`, SRM.data.screened.noPRTC$`G017-B`)
G031 <- ave(SRM.data.screened.noPRTC$`G031-A`, SRM.data.screened.noPRTC$`G031-B`, SRM.data.screened.noPRTC$`G031-C`)
G032 <- ave(SRM.data.screened.noPRTC$`G032-A`, SRM.data.screened.noPRTC$`G032-B`)
G040 <- ave(SRM.data.screened.noPRTC$`G040-A`, SRM.data.screened.noPRTC$`G040-B`)
G041 <- ave(SRM.data.screened.noPRTC$`G041-A`, SRM.data.screened.noPRTC$`G041-B`)
G042 <- ave(SRM.data.screened.noPRTC$`G042-A`, SRM.data.screened.noPRTC$`G042-B`) #C removed
G043 <- ave(SRM.data.screened.noPRTC$`G043-A`, SRM.data.screened.noPRTC$`G043-B`)
G045 <- ave(SRM.data.screened.noPRTC$`G045-A`, SRM.data.screened.noPRTC$`G045-B`)
G047 <- ave(SRM.data.screened.noPRTC$`G047-A`, SRM.data.screened.noPRTC$`G047-B`)
G049 <- ave(SRM.data.screened.noPRTC$`G049-A`, SRM.data.screened.noPRTC$`G049-B`)
G053 <- ave(SRM.data.screened.noPRTC$`G053-A`, SRM.data.screened.noPRTC$`G053-remake-C`, SRM.data.screened.noPRTC$`G053-remake-D`) #B removed 
G054 <- ave(SRM.data.screened.noPRTC$`G054-A`, SRM.data.screened.noPRTC$`G054-B`)
G055 <- ave(SRM.data.screened.noPRTC$`G055-A`, SRM.data.screened.noPRTC$`G055-B`, SRM.data.screened.noPRTC$`G055-C`)
G057 <- ave(SRM.data.screened.noPRTC$`G057-A`, SRM.data.screened.noPRTC$`G057-C`) #B removed
G060 <- ave(SRM.data.screened.noPRTC$`G060-A`, SRM.data.screened.noPRTC$`G060-B`)
G062 <- ave(SRM.data.screened.noPRTC$`G062-B`, SRM.data.screened.noPRTC$`G062-C`)
G064 <- ave(SRM.data.screened.noPRTC$`G064-A`, SRM.data.screened.noPRTC$`G064-B`)
G066 <- ave(SRM.data.screened.noPRTC$`G066-A`, SRM.data.screened.noPRTC$`G066-B`)
G070 <- ave(SRM.data.screened.noPRTC$`G070-A`, SRM.data.screened.noPRTC$`G070-B`, SRM.data.screened.noPRTC$`G070-C`)
G071.A <- ave(SRM.data.screened.noPRTC$`G071-A-A`, SRM.data.screened.noPRTC$`G071-A-B`)
G071.B <- ave(SRM.data.screened.noPRTC$`G071-B-A`, SRM.data.screened.noPRTC$`G071-B-B`)
G073 <- ave(SRM.data.screened.noPRTC$`G073-A`, SRM.data.screened.noPRTC$`G073-B`, SRM.data.numeric$`G073-C`)
G074 <- ave(SRM.data.screened.noPRTC$`G074-A`, SRM.data.screened.noPRTC$`G074-B`)
G079 <- ave(SRM.data.screened.noPRTC$`G079-A`, SRM.data.screened.noPRTC$`G079-B`)
G081 <- ave(SRM.data.screened.noPRTC$`G081-A`, SRM.data.screened.noPRTC$`G081-B`)
G104 <- ave(SRM.data.screened.noPRTC$`G104-A`, SRM.data.screened.noPRTC$`G104-remake-C`, SRM.data.screened.noPRTC$`G104-remake-D`) #B removed
G105 <- ave(SRM.data.screened.noPRTC$`G105-A`, SRM.data.screened.noPRTC$`G105-B`)
G109 <- ave(SRM.data.screened.noPRTC$`G109-A`, SRM.data.screened.noPRTC$`G109-C`)
G114 <- ave(SRM.data.screened.noPRTC$`G114-A`, SRM.data.screened.noPRTC$`G114-B`, SRM.data.screened.noPRTC$`G114-remake-C`, SRM.data.screened.noPRTC$`G114-remake-D`)
G116 <- ave(SRM.data.screened.noPRTC$`G116-A`, SRM.data.screened.noPRTC$`G116-B`)
G120 <- ave(SRM.data.screened.noPRTC$`G120-A`, SRM.data.screened.noPRTC$`G120-B`)
G122 <- ave(SRM.data.screened.noPRTC$`G122-A`, SRM.data.screened.noPRTC$`G122-B`)
G127 <- ave(SRM.data.screened.noPRTC$`G127-A`, SRM.data.screened.noPRTC$`G127-C`) #B removed
G128 <- ave(SRM.data.screened.noPRTC$`G128-A`, SRM.data.screened.noPRTC$`G128-C`,SRM.data.screened.noPRTC$`G128-D`)
G129 <- ave(SRM.data.screened.noPRTC$`G129-A`, SRM.data.screened.noPRTC$`G129-B`)
G132 <- ave(SRM.data.screened.noPRTC$`G132-A`, SRM.data.screened.noPRTC$`G132-C`, SRM.data.screened.noPRTC$`G132-D`)
# Tech reps removed: 3C, 42C, 53B, 57B, 104B, 127B

# Entire sample removed: 73
PG.E.samples <- PG.E.samples[!PG.E.samples %in% "G073"] #revised PG.E.sample list 

SRM.data.mean <- cbind.data.frame(rownames(SRM.data.screened.noPRTC), G001,G002,G003,G007,G008,G009,G110,G012,G013,G014,G015,G016,G017,G031,G032,G040,G041,G042,G043,G045,G047,G049,G053,G054,G055,G057,G060,G062,G064,G066,G070,G071.A,G071.B,G074,G079,G081,G104,G105,G109,G114,G116,G120,G122,G127,G128,G129,G132)
SRM.data.mean <- data.frame(SRM.data.mean[,-1], row.names=SRM.data.mean[,1]) #make first column row names, and delete first column
write.csv(SRM.data.mean, file="Analyses/2017-September_SRM-results/2017-09-04_SRM-data-meanBYsample.csv")
View(SRM.data.mean)

require(plotrix)
# Standard error for tech reps ###FYI THIS IS NOT WORKING !!!!!!!!!!!!
G001.err <- std.error(c(SRM.data.screened.noPRTC$`G001-A`, SRM.data.screened.noPRTC$`G001-B`))
G002.err <- std.error(c(SRM.data.screened.noPRTC$`G002-A`, SRM.data.screened.noPRTC$`G002-B`, SRM.data.screened.noPRTC$`G002-C`))
G003.err <- std.error(c(SRM.data.screened.noPRTC$`G003-A`, SRM.data.screened.noPRTC$`G003-B`)) #C removed
G007.err <- std.error(c(SRM.data.screened.noPRTC$`G007-A`, SRM.data.screened.noPRTC$`G007-B`))
G008.err <- std.error(c(SRM.data.screened.noPRTC$`G008-A`, SRM.data.screened.noPRTC$`G008-B`))
G009.err <- std.error(c(SRM.data.screened.noPRTC$`G009-A`, SRM.data.screened.noPRTC$`G009-B`))
G110.err <- std.error(c(SRM.data.screened.noPRTC$`G110-A`, SRM.data.screened.noPRTC$`G110-B`))
G012.err <- std.error(c(SRM.data.screened.noPRTC$`G012-A`, SRM.data.screened.noPRTC$`G012-B`, SRM.data.screened.noPRTC$`G012-C`))
G013.err <- std.error(c(SRM.data.screened.noPRTC$'G013-A', SRM.data.screened.noPRTC$'G013-C'))
G014.err <- std.error(c(SRM.data.screened.noPRTC$`G014-A`, SRM.data.screened.noPRTC$`G014-B`))
G015.err <- std.error(c(SRM.data.screened.noPRTC$`G015-A`, SRM.data.screened.noPRTC$`G015-B`))
G016.err <- std.error(c(SRM.data.screened.noPRTC$`G016-A`, SRM.data.screened.noPRTC$`G016-B`, SRM.data.screened.noPRTC$`G016-C`))
G017.err <- std.error(c(SRM.data.screened.noPRTC$`G017-A`, SRM.data.screened.noPRTC$`G017-B`))
G031.err <- std.error(c(SRM.data.screened.noPRTC$`G031-A`, SRM.data.screened.noPRTC$`G031-B`, SRM.data.screened.noPRTC$`G031-C`))
G032.err <- std.error(c(SRM.data.screened.noPRTC$`G032-A`, SRM.data.screened.noPRTC$`G032-B`))
G040.err <- std.error(c(SRM.data.screened.noPRTC$`G040-A`, SRM.data.screened.noPRTC$`G040-B`))
G041.err <- std.error(c(SRM.data.screened.noPRTC$`G041-A`, SRM.data.screened.noPRTC$`G041-B`))
G042.err <- std.error(c(SRM.data.screened.noPRTC$`G042-A`, SRM.data.screened.noPRTC$`G042-B`)) #C removed
G043.err <- std.error(c(SRM.data.screened.noPRTC$`G043-A`, SRM.data.screened.noPRTC$`G043-B`))
G045.err <- std.error(c(SRM.data.screened.noPRTC$`G045-A`, SRM.data.screened.noPRTC$`G045-B`))
G047.err <- std.error(c(SRM.data.screened.noPRTC$`G047-A`, SRM.data.screened.noPRTC$`G047-B`))
G049.err <- std.error(c(SRM.data.screened.noPRTC$`G049-A`, SRM.data.screened.noPRTC$`G049-B`))
G053.err <- std.error(c(SRM.data.screened.noPRTC$`G053-A`, SRM.data.screened.noPRTC$`G053-remake-C`, SRM.data.screened.noPRTC$`G053-remake-D`)) #B removed 
G054.err <- std.error(c(SRM.data.screened.noPRTC$`G054-A`, SRM.data.screened.noPRTC$`G054-B`))
G055.err <- std.error(c(SRM.data.screened.noPRTC$`G055-A`, SRM.data.screened.noPRTC$`G055-B`, SRM.data.screened.noPRTC$`G055-C`))
G057.err <- std.error(c(SRM.data.screened.noPRTC$`G057-A`, SRM.data.screened.noPRTC$`G057-C`)) #B removed
G060.err <- std.error(c(SRM.data.screened.noPRTC$`G060-A`, SRM.data.screened.noPRTC$`G060-B`))
G062.err <- std.error(c(SRM.data.screened.noPRTC$`G062-B`, SRM.data.screened.noPRTC$`G062-C`))
G064.err <- std.error(c(SRM.data.screened.noPRTC$`G064-A`, SRM.data.screened.noPRTC$`G064-B`))
G066.err <- std.error(c(SRM.data.screened.noPRTC$`G066-A`, SRM.data.screened.noPRTC$`G066-B`))
G070.err <- std.error(c(SRM.data.screened.noPRTC$`G070-A`, SRM.data.screened.noPRTC$`G070-B`, SRM.data.screened.noPRTC$`G070-C`))
G071.A.err <- std.error(cor(SRM.data.screened.noPRTC$`G071-A-A`, SRM.data.screened.noPRTC$`G071-A-B`))
G071.B.err <- std.error(cor(SRM.data.screened.noPRTC$`G071-B-A`, SRM.data.screened.noPRTC$`G071-B-B`))
G073.err <- std.error(c(SRM.data.screened.noPRTC$`G073-A`, SRM.data.screened.noPRTC$`G073-B`, SRM.data.numeric$`G073-C`))
G074.err <- std.error(c(SRM.data.screened.noPRTC$`G074-A`, SRM.data.screened.noPRTC$`G074-B`))
G079.err <- std.error(c(SRM.data.screened.noPRTC$`G079-A`, SRM.data.screened.noPRTC$`G079-B`))
G081.err <- std.error(c(SRM.data.screened.noPRTC$`G081-A`, SRM.data.screened.noPRTC$`G081-B`))
G104.err <- std.error(c(SRM.data.screened.noPRTC$`G104-A`, SRM.data.screened.noPRTC$`G104-remake-C`, SRM.data.screened.noPRTC$`G104-remake-D`)) #B removed
G105.err <- std.error(c(SRM.data.screened.noPRTC$`G105-A`, SRM.data.screened.noPRTC$`G105-B`))
G109.err <- std.error(c(SRM.data.screened.noPRTC$`G109-A`, SRM.data.screened.noPRTC$`G109-C`))
G114.err <- std.error(c(SRM.data.screened.noPRTC$`G114-A`, SRM.data.screened.noPRTC$`G114-B`, SRM.data.screened.noPRTC$`G114-remake-C`, SRM.data.screened.noPRTC$`G114-remake-D`))
G116.err <- std.error(c(SRM.data.screened.noPRTC$`G116-A`, SRM.data.screened.noPRTC$`G116-B`))
G120.err <- std.error(c(SRM.data.screened.noPRTC$`G120-A`, SRM.data.screened.noPRTC$`G120-B`))
G122.err <- std.error(c(SRM.data.screened.noPRTC$`G122-A`, SRM.data.screened.noPRTC$`G122-B`))
G127.err <- std.error(c(SRM.data.screened.noPRTC$`G127-A`, SRM.data.screened.noPRTC$`G127-C`)) #B removed
G128.err <- std.error(c(SRM.data.screened.noPRTC$`G128-A`, SRM.data.screened.noPRTC$`G128-C`,SRM.data.screened.noPRTC$`G128-D`))
G129.err <- std.error(c(SRM.data.screened.noPRTC$`G129-A`, SRM.data.screened.noPRTC$`G129-B`))
G132.err <- std.error(c(SRM.data.screened.noPRTC$`G132-A`, SRM.data.screened.noPRTC$`G132-C`, SRM.data.screened.noPRTC$`G132-D`))

SRM.data.stderr <- cbind.data.frame(rownames(SRM.data.screened.noPRTC), G001.err,G002.err,G003.err,G007.err,G008.err,G009.err,G110.err,G012.err,G013.err,G014.err,G015.err,G016.err,G017.err,G031.err,G032.err,G040.err,G041.err,G042.err,G043.err,G045.err,G047.err,G049.err,G053.err,G054.err,G055.err,G057.err,G060.err,G062.err,G064.err,G066.err,G070.err,G071.A.err,G071.B.err,G073.err, G074.err,G079.err,G081.err,G104.err,G105.err,G109.err,G114.err,G116.err,G120.err,G122.err,G127.err,G128.err,G129.err,G132.err)
