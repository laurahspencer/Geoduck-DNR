### IMPORT DATASET ### 
### Note: this script works with datasets downloaded from Skyline.  Skyline data is exported as report, and this script works with data with the following metrics: Protein Name, Transitions, Peptide Sequence, Fragment Ion, Peptide Retention Time, Area

setwd("~/Documents/Roberts Lab/Geoduck-DNR/") 
SRMreport <- read.csv("Data/2017-08-11_Transition Results_LHS modified-noRT-pivoted.csv", header=FALSE, na.strings = "#N/A", stringsAsFactors = FALSE) # import local file
head(SRMreport)
SRMsequence <- read.csv("Data/2017-07-28_SRM-Sequence-final.csv", header=TRUE, stringsAsFactors = FALSE)
head(SRMsequence)
sample.key <- read.csv("Data/2017-08-14-Geoduck-samples.csv", header=TRUE, stringsAsFactors = FALSE)
head(sample.key)
SRMsamples <- noquote(as.character(c("G013", "G120", "G047", "G017", "G079", "G127", "G060", "G009", "G002", "G128", "G016", "G071-A", "G114", "G045", "G132", "G031", "G012", "G116", "G043", "G015", "G040", "G110", "G008", "G109", "G122", "G041", "G066", "G105", "G032", "G129", "G054", "G081", "G003", "G074", "G014", "G049", "G053", "G104", "G055", "G042", "G064", "G073", "G057", "G007", "G070", "G001", "G071-B", "G062")))

### REPLACE REP NAMES WITH SAMPLE NAMES ###
# I could also probably use the function merge() for this task

SRMreport[1,] # view replicate names
length(SRMreport[1,]) # Number of replicates I ran on mass spec
rep.names <- SRMreport[1,] # create vector of replicate names
rep.names.short <- noquote(gsub(' Area', '', rep.names)) # remove Area from rep name, and don't include quotes 
rep.names.short # check to confirm correct names
rep.names.short <- noquote(gsub('2017_July_10_bivalves_', '', rep.names.short)) #remove the extra long rep name that is a residual from the .raw file name
length(rep.names.short)
SRMsequence$Sample...rep.name
noquote(as.character(SRMsequence$Sample...rep.name))
repsTOsamples <- as.data.frame(SRMsequence[,c(2,3,5)])
library(dplyr)
repsTOsamples.filtered <- filter(repsTOsamples, repsTOsamples[,1] %in% rep.names.short)
samples <- as.character(repsTOsamples.filtered$Sample...rep.name)
other.headers <- as.character(rep.names.short[1:4])
samples.vector <- noquote(c(other.headers, samples, stringsAsFactors = FALSE))
samples.vector <- samples.vector[-121]
SRM.data <- SRMreport
SRM.data[1,] <- samples.vector
ncol(SRM.data) #confirm still have the correct # columns
colnames(SRM.data) <- SRM.data[1,] #make first row column names

### ANNOTATE SAMPLE NAMES WITH SITE & TREATMENT ########## needs work

head(sample.key)
sample.key[,c(8,9)]
repsTOsamples.filtered.annotated <- filter(sample.key[,c(8,9)], sample.key$PRVial %in% repsTOsamples.filtered$Comment) #pull site & treatment from sample key
length(SRMsamples) # check # samples there should be total
nrow(repsTOsamples.filtered.annotated) # check # samples is appropriate
repsTOsamples.filtered.annotated # NOTE: missing 71-A & 71-B, need to fix 

# Add G071-A & G071-B coding to the annotated key
s71.A <- data.frame(matrix(0, ncol=2, nrow=1))
s71.A[1,1] <- "G071.A"
s71.A[1,2] <- "PG-E"
colnames(s71.A) <- colnames(repsTOsamples.filtered.annotated)
s71.B <- data.frame(matrix(0, ncol=2, nrow=1))
s71.B[1,1] <- "G071.B"
s71.B[1,2] <- "PG-E"
colnames(s71.B) <- colnames(repsTOsamples.filtered.annotated)

# row bind annotated key w/ 71 info
sample.key.annotated <- rbind(repsTOsamples.filtered.annotated, s71.A, s71.B)

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

### CONVERT AREA DATA TO NUMERIC FORMAT ### 

SRM.data.numeric <- SRM.data # First, change all area values to numeric, so I can average, etc. I know that my area data is from column 5 to 120 
SRM.data.numeric[,5:120] <- as.numeric( 
  as.character(
    unlist(
      SRM.data.numeric[,5:120])
  )
)
is.numeric(SRM.data.numeric[2,20]) # confirm area data is numeric, using a random cell. Should equal TRUE.

#### NAME EACH ROW WITH A UNIQUE TRANSITION ID ###
nTransitions <- length(SRM.data.numeric$Transition) # How many transitions are there
Transition.ID <- vector(length=nTransitions) # create empty vector with length= number of transitions
for (i in 1:nTransitions) {  
  Transition.ID[i] <- paste(SRM.data.numeric[i,3], SRM.data.numeric[i,4])}  # loop that fills empty vector with unique transition ID, built from the peptide sequence (column 3) and the fragment ion (columm 4)
Transition.ID # confirm correctly named transition IDs
length(SRM.data.numeric$Transition) == length(Transition.ID) # confirm that I didn't lose any transitions
row.names(SRM.data.numeric) <- Transition.ID # assign newly created transition IDs as row names
head(SRM.data.numeric) # confirm changes
# write.csv(SRM.data.numeric, file="Data/2017-08-19_SRM-data-R1.csv") #write this file out for safe keeping

### NORMALIZE BASED ON PRTC ABUNDANCE ### 

# Extract PRTC peptides 
SRM.data.numeric[115:142,1:5] #check out which rows pertain to PRTC peptides, it's 117->142
SRM.PRTC <- SRM.data.numeric[117:142,] #pull out PRTC data for each sample
head(SRM.PRTC)
ncol(SRM.PRTC[,-1:-4])

# Discard PRTC peptides what eluted early and late for this analysis since their quantities are less stable. These were identified via Skyline.
# Good quality are: "LTILEELR", "GLILVGGYGTR", "ELGQSGVDTYLQTK", "SAAGAFGPELSR", "TASEFDSAIAQDK"
# Poor Quality are: "IGDYAGIK", "DIPVPKPK", "HVLTSIGEK", "GISNEGQNASIK", "SSAAPPPPPR"
rownames(SRM.PRTC)
rownames(SRM.PRTC[1:12,])
rownames(SRM.PRTC[13:27,])
SRM.PRTC.good <- SRM.PRTC[1:12,]
head(SRM.PRTC.good)
ncol(SRM.PRTC.good[,-1:-4])

###### Assign PRTC batches based on which PRTC mix samples were made with ######

SRM.PRTC.a <- cbind(
  SRM.PRTC.good[,1:4], 
  SRM.PRTC.good[, grepl("G013|G120|G047|G017|G079|G127|G060|G009|G002|G128|G016|G071-A|G114|G045|G132|G031|G012|G116|G043|G015|G040", names(SRM.PRTC)) & !grepl("remake", names(SRM.PRTC))]
)

SRM.PRTC.b <- cbind(
  SRM.PRTC.good[,1:4],
  SRM.PRTC.good[, grepl("G110|G008|G109|G122", names(SRM.PRTC))]
)

SRM.PRTC.c <- cbind(
  SRM.PRTC.good[,1:4],
  SRM.PRTC.good[, grepl("G041|G066|G105|G032|G129|G054|G081|G003|G074|G014|G049|G053|G104|G055|G042|G064|G073|G057|G007|G070|G001|G071-B|G062", names(SRM.PRTC)) & !grepl("remake", names(SRM.PRTC))]
)

SRM.PRTC.d <- cbind(
  SRM.PRTC.good[,1:4],
  SRM.PRTC.good[, grepl("G114.remake|G053.remake|G104.remake", names(SRM.PRTC))]
)

############# Calculate mean abundance for each transition within PRTC batches ##############

# PRTC batch A
SRM.PRTC.a.mean <- 1:nrow(SRM.PRTC.a)# Create vector to be filled by loop with mean PRTC area
for (i in SRM.PRTC.a.mean) { 
  SRM.PRTC.a.mean[i] <- rowMeans(SRM.PRTC.a[i,-1:-4], na.rm=TRUE)  # create a sample, each entry represents a sample with mean abundance for all PRTC transitions in that sample
}
print(SRM.PRTC.a.mean)

# PRTC batch B
SRM.PRTC.b.mean <- 1:nrow(SRM.PRTC.b)# Create vector to be filled by loop with mean PRTC area
for (i in SRM.PRTC.b.mean) { 
  SRM.PRTC.b.mean[i] <- rowMeans(SRM.PRTC.b[i,-1:-4], na.rm=TRUE)  # create a sample, each entry represents a sample with mean abundance for all PRTC transitions in that sample
}
print(SRM.PRTC.b.mean)

# PRTC batch C
SRM.PRTC.c.mean <- 1:nrow(SRM.PRTC.c)# Create vector to be filled by loop with mean PRTC area
for (i in SRM.PRTC.c.mean) { 
  SRM.PRTC.c.mean[i] <- rowMeans(SRM.PRTC.c[i,-1:-4], na.rm=TRUE)  # create a sample, each entry represents a sample with mean abundance for all PRTC transitions in that sample
}
print(SRM.PRTC.c.mean)

# PRTC batch D 
SRM.PRTC.d.mean <- 1:nrow(SRM.PRTC.d)# Create vector to be filled by loop with mean PRTC area
for (i in SRM.PRTC.d.mean) { 
  SRM.PRTC.d.mean[i] <- rowMeans(SRM.PRTC.d[i,-1:-4], na.rm=TRUE)  # create a sample, each entry represents a sample with mean abundance for all PRTC transitions in that sample
}
print(SRM.PRTC.d.mean)

# Combine batch means into new data frame
PRTC.batch.means <- cbind(SRM.PRTC.good[,1:4], SRM.PRTC.a.mean, SRM.PRTC.b.mean, SRM.PRTC.c.mean, SRM.PRTC.d.mean)
# write.csv(PRTC.batch.means, file="Data/2017-08-22_SRM-PRTC-batch-means.csv") #write this file out for safe keeping

# divide batches a, b & c by b, since it has the highest abundance
PRTC.a.ratio <- mean(PRTC.batch.means$SRM.PRTC.a.mean/PRTC.batch.means$SRM.PRTC.b.mean, na.rm = TRUE )
PRTC.b.ratio <- mean(PRTC.batch.means$SRM.PRTC.b.mean/PRTC.batch.means$SRM.PRTC.b.mean, na.rm = TRUE )
PRTC.c.ratio <- mean(PRTC.batch.means$SRM.PRTC.c.mean/PRTC.batch.means$SRM.PRTC.b.mean, na.rm = TRUE )
PRTC.d.ratio <- mean(PRTC.batch.means$SRM.PRTC.d.mean/PRTC.batch.means$SRM.PRTC.b.mean, na.rm = TRUE )

# The following is average relative abundance to batch B
PRTC.a.ratio
PRTC.b.ratio
PRTC.c.ratio
PRTC.d.ratio

####  ADJUST PRTC abundance based on batch ratios ###

SRM.PRTC.adjusted <- cbind.data.frame(
  SRM.PRTC.good[, grepl("G013|G120|G047|G017|G079|G127|G060|G009|G002|G128|G016|G071-A|G114|G045|G132|G031|G012|G116|G043|G015|G040", names(SRM.PRTC.good)) & !grepl("remake", names(SRM.PRTC.good))]/PRTC.a.ratio,
  SRM.PRTC.good[,grepl(c("G110|G008|G109|G122"), names(SRM.PRTC.good))& !grepl("remake", names(SRM.PRTC.good))]/PRTC.b.ratio,
  SRM.PRTC.good[,grepl("G041|G066|G105|G032|G129|G054|G081|G003|G074|G014|G049|G053|G104|G055|G042|G064|G073|G057|G007|G070|G001|G071-B|G062", names(SRM.PRTC.good))& !grepl("remake", names(SRM.PRTC.good))]/PRTC.c.ratio,
  SRM.PRTC.good[, grepl("G114.remake|G053.remake|G104.remake", names(SRM.PRTC.good))]/PRTC.d.ratio
)
ncol(SRM.PRTC.adjusted) == ncol(SRM.PRTC.good[,-1:-4])

# Plot transition abundance pre- and post- adjustment
png("Analyses/2017-August_SRM-Analysis/PRTC-pre-post-adjusted-1.png", width = 1000, height = 1200, units = "px")
par(mfrow=c(6,2))
barplot(as.matrix(SRM.PRTC.adjusted[1,]), main=rownames(SRM.PRTC.adjusted[1,])) 
barplot(as.matrix(SRM.PRTC.good[1,-1:-4]))
barplot(as.matrix(SRM.PRTC.adjusted[2,]), main=rownames(SRM.PRTC.adjusted[2,]))
barplot(as.matrix(SRM.PRTC.good[2,-1:-4])) 
barplot(as.matrix(SRM.PRTC.adjusted[3,]), main=rownames(SRM.PRTC.adjusted[3,]))
barplot(as.matrix(SRM.PRTC.good[3,-1:-4])) 
barplot(as.matrix(SRM.PRTC.adjusted[4,]), main=rownames(SRM.PRTC.adjusted[4,])) 
barplot(as.matrix(SRM.PRTC.good[4,-1:-4])) 
barplot(as.matrix(SRM.PRTC.adjusted[5,]), main=rownames(SRM.PRTC.adjusted[5,])) 
barplot(as.matrix(SRM.PRTC.good[5,-1:-4])) 
barplot(as.matrix(SRM.PRTC.adjusted[6,]), main=rownames(SRM.PRTC.adjusted[6,])) 
barplot(as.matrix(SRM.PRTC.good[6,-1:-4])) 
dev.off()

png("Analyses/2017-August_SRM-Analysis/PRTC-pre-post-adjusted-2.png", width = 1000, height = 1200, units = "px")
par(mfrow=c(6,2))
barplot(as.matrix(SRM.PRTC.adjusted[7,]), main=rownames(SRM.PRTC.adjusted[7,])) 
barplot(as.matrix(SRM.PRTC.good[7,-1:-4]))
barplot(as.matrix(SRM.PRTC.adjusted[8,]), main=rownames(SRM.PRTC.adjusted[8,]))
barplot(as.matrix(SRM.PRTC.good[8,-1:-4])) 
barplot(as.matrix(SRM.PRTC.adjusted[9,]), main=rownames(SRM.PRTC.adjusted[9,]))
barplot(as.matrix(SRM.PRTC.good[9,-1:-4])) 
barplot(as.matrix(SRM.PRTC.adjusted[10,]), main=rownames(SRM.PRTC.adjusted[10,])) 
barplot(as.matrix(SRM.PRTC.good[10,-1:-4])) 
barplot(as.matrix(SRM.PRTC.adjusted[11,]), main=rownames(SRM.PRTC.adjusted[11,])) 
barplot(as.matrix(SRM.PRTC.good[11,-1:-4])) 
barplot(as.matrix(SRM.PRTC.adjusted[12,]), main=rownames(SRM.PRTC.adjusted[12,])) 
barplot(as.matrix(SRM.PRTC.good[12,-1:-4])) 
dev.off()

#### COMPARE PRE- AND POST- BATCH ADJUSTMENT COEFFICIENTO F VARIATION FOR PRTC TRANSITIONS

# Transpose PRTC peptide data to run stats
SRM.PRTC.t <- t(SRM.PRTC.good[,-1:-4]) #All samples
SRM.PRTC.a.t <- t(SRM.PRTC.a[,-1:-4]) #Batch A
SRM.PRTC.b.t <- t(SRM.PRTC.b[,-1:-4]) #Batch B
SRM.PRTC.c.t <- t(SRM.PRTC.c[,-1:-4]) #Batch C
SRM.PRTC.d.t <- t(SRM.PRTC.d[,-1:-4]) #Batch D
SRM.PRTC.adjusted.t <- t(SRM.PRTC.adjusted) #All samples, adjusted

# Loop prepares a vector with Coefficient of Variance for each PRTC transition
SRM.PRTC.t.cv <- 1:ncol(SRM.PRTC.t) #All samples together
for (i in SRM.PRTC.t.cv) { 
  SRM.PRTC.t.cv[i] <- sd(SRM.PRTC.t[,i], na.rm = TRUE)/mean(SRM.PRTC.t[,i], na.rm = TRUE)*100
}
SRM.PRTC.a.t.cv <- 1:ncol(SRM.PRTC.a.t) #Batch A
for (i in SRM.PRTC.a.t.cv) { 
  SRM.PRTC.a.t.cv[i] <- sd(SRM.PRTC.a.t[,i], na.rm = TRUE)/mean(SRM.PRTC.a.t[,i], na.rm = TRUE)*100
}
SRM.PRTC.b.t.cv <- 1:ncol(SRM.PRTC.b.t) #Batch B
for (i in SRM.PRTC.b.t.cv) { 
  SRM.PRTC.b.t.cv[i] <- sd(SRM.PRTC.b.t[,i], na.rm = TRUE)/mean(SRM.PRTC.b.t[,i], na.rm = TRUE)*100
}
SRM.PRTC.c.t.cv <- 1:ncol(SRM.PRTC.c.t) #Batch C
for (i in SRM.PRTC.c.t.cv) { 
  SRM.PRTC.c.t.cv[i] <- sd(SRM.PRTC.c.t[,i], na.rm = TRUE)/mean(SRM.PRTC.c.t[,i], na.rm = TRUE)*100
}
SRM.PRTC.d.t.cv <- 1:ncol(SRM.PRTC.d.t) #Batch D
for (i in SRM.PRTC.d.t.cv) { 
  SRM.PRTC.d.t.cv[i] <- sd(SRM.PRTC.d.t[,i], na.rm = TRUE)/mean(SRM.PRTC.d.t[,i], na.rm = TRUE)*100
}
SRM.PRTC.adjusted.t.cv <- 1:ncol(SRM.PRTC.adjusted.t) #All samples together
for (i in SRM.PRTC.adjusted.t.cv) { 
  SRM.PRTC.adjusted.t.cv[i] <- sd(SRM.PRTC.adjusted.t[,i], na.rm = TRUE)/mean(SRM.PRTC.adjusted.t[,i], na.rm = TRUE)*100
}

#Create dataframe with coefficient of variation for each transition, for all samples then by batches (a,b,c,d)
PRTC.cv.comparison <- cbind.data.frame(transitions=colnames(SRM.PRTC.t), Not.Adjusted=as.integer(SRM.PRTC.t.cv), A=as.integer(SRM.PRTC.a.t.cv), B=as.integer(SRM.PRTC.b.t.cv), C=as.integer(SRM.PRTC.c.t.cv), D=as.integer(SRM.PRTC.d.t.cv), Adjusted=as.integer(SRM.PRTC.adjusted.t.cv))
PRTC.cv.comparison #this is the result

### NORMALIZE ASSAY DATA BASED ON PRTC ABUNDANCE 

# calculate mean abundance fo all PRTC transitions within samples
SRM.PRTC.adjusted.mean <- data.frame(colMeans(SRM.PRTC.adjusted, na.rm=TRUE))

# Now, normalize all sample abundance data based on PRTC mean abundances
SRM.data.numeric.1 <- SRM.data.numeric[c(-1,-117:-142),-1:-4] #remove non-abundance-data  columns, remove PRTC transitions from assay data
PRTC.norm.vector <- SRM.PRTC.adjusted.mean[,1] #create vector of mean PRTC abundances for each sample
SRM.PRTC.adjusted.mean
length(PRTC.norm.vector) == ncol(SRM.data.numeric.1) #confirm PRTC normalization vector length equals # samples in srm data
SRM.data.normalized <- sweep(SRM.data.numeric.1, 2, PRTC.norm.vector, "/") #normalize srm data (averaged by tech. rep) by mean PRTC abundance for that sample

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
# comm= your data.frame or matrix
# distance= bray, (not sure what this means)
# k= # of dimensions to assess
# trymax = max # iterations to attempt if no solution is reached

# Create Shepard plot, which shows scatter around the regression between the interpoint distances in the final configuration (i.e., the distances between each pair of communities) against their original dissimilarities.
stressplot(SRM.mean.nmds) 

#Make figure

plot(SRM.mean.nmds)
# site (sample) in black circle
# species (variable) in red ticks

# make figure with sample annotations https://stat.ethz.ch/pipermail/r-sig-ecology/2011-September/002371.html
SRM.nmds.mean.samples <- scores(SRM.mean.nmds, display = "sites")
SRM.nmds.mean.transitions <- scores(SRM.nmds, display = "species")
# this probably isn't necessary

### Let's plot using ordiplot()

library(RColorBrewer)
marker = c(color = brewer.pal(4, "Set1"))

png("Analyses/2017-August_SRM-Analysis/2017-08-24_SRM-NMDS-plot.png")
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
dev.off()

#### Create plot with forced aspect ratio to zoom in ### 

png("Analyses/2017-August_SRM-Analysis/2017-08-24_SRM-NMDS-plot-zoomed.png")
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
dev.off()

#### CREATE NMDS PLOT, MEAN OF TECH REPS - LOG TRANSFORMED ########

#Transpose the file so that rows and columns are switched and normalized by log(x+1)
SRM.data.mean.t.log <- SRM.data.mean.t
SRM.data.mean.t.log[is.na(SRM.data.mean.t.log)] <- 0
SRM.data.mean.t.log <- (SRM.data.mean.t.log+1)
SRM.data.mean.t.log <- data.trans(SRM.data.mean.t.log, method = 'log', plot = FALSE)

#Make MDS dissimilarity matrix
#
SRM.mean.log.nmds <- metaMDS(SRM.data.mean.t.log, distance = 'bray', k = 2, trymax = 3000, autotransform = FALSE)
# comm= your data.frame or matrix
# distance= bray, (not sure what this means)
# k= # of dimensions to assess
# trymax = max # iterations to attempt if no solution is reached

# Create Shepard plot, which shows scatter around the regression between the interpoint distances in the final configuration (i.e., the distances between each pair of communities) against their original dissimilarities.
stressplot(SRM.mean.log.nmds) 

#Make figure

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

png("Analyses/2017-August_SRM-Analysis/2017-08-24_SRM-NMDS-log-plot.png")
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
dev.off()

#### Create plot with forced aspect ratio to zoom in ### 

png("Analyses/2017-August_SRM-Analysis/2017-08-24_SRM-NMDS-plot-log-zoomed.png")
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
dev.off()

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
head(SRM.data.mean.stringent)

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

png("Analyses/2017-August_SRM-Analysis/2017-08-24_SRM-NMDS-stringent-plot.png")
ordiplot(SRM.mean.s.nmds, type="n")
points(SRM.nmds.mean.s.samples[c(CI.B.samples.s),], col=marker[1], pch=8)
points(SRM.nmds.mean.s.samples[c(CI.E.samples.s),], col=marker[1], pch=15)
points(SRM.nmds.mean.s.samples[c(PG.B.samples.s),], col=marker[2], pch=8)
points(SRM.nmds.mean.s.samples[c(PG.E.samples.s),], col=marker[2], pch=15)
points(SRM.nmds.mean.s.samples[c(WB.B.samples.s),], col=marker[3], pch=8)
points(SRM.nmds.mean.s.samples[c(WB.E.samples.s),], col=marker[3], pch=15)
points(SRM.nmds.mean.s.samples[c(FB.B.samples.s),], col=marker[4], pch=8)
points(SRM.nmds.mean.s.samples[c(FB.E.samples.s),], col=marker[4], pch=15)
dev.off()

legend(1.5,0.7, pch=c(rep(16,4), 8, 15), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Port Gamble", "Bare", "Eelgrass"), col=c(marker[1], marker[2], marker[3], marker[4], "black", "black"))

#### Create plot with forced aspect ratio to zoom in ### 

png("Analyses/2017-August_SRM-Analysis/2017-08-24_SRM-NMDS-stringent-plot-zoomed.png")
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
dev.off()

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

#ANOSIM
rownames(SRM.data.mean.stringent.t) <- sub(".s", "", rownames(SRM.data.mean.stringent.t))

CI.sb <- data.frame(SAMPLE=CI.B.samples.s, SITE=rep("CI", times=length(CI.B.samples.s)), TREATMENT=rep("Bare", times=length(CI.B.samples.s)))
CI.se <- data.frame(SAMPLE=CI.E.samples.s, SITE=rep("CI", times=length(CI.E.samples.s)), TREATMENT=rep("Eelgrass", times=length(CI.E.samples.s)))
PG.sb <- data.frame(SAMPLE=PG.B.samples.s, SITE=rep("PG", times=length(PG.B.samples.s)), TREATMENT=rep("Bare", times=length(PG.B.samples.s)))
PG.se <- data.frame(SAMPLE=PG.E.samples.s, SITE=rep("PG", times=length(PG.E.samples.s)), TREATMENT=rep("Eelgrass", times=length(PG.E.samples.s)))
WB.sb <- data.frame(SAMPLE=WB.B.samples.s, SITE=rep("WB", times=length(WB.B.samples.s)), TREATMENT=rep("Bare", times=length(WB.B.samples.s)))
WB.se <- data.frame(SAMPLE=WB.E.samples.s, SITE=rep("WB", times=length(WB.E.samples.s)), TREATMENT=rep("Eelgrass", times=length(WB.E.samples.s)))
FB.sb <- data.frame(SAMPLE=FB.B.samples.s, SITE=rep("FB", times=length(FB.B.samples.s)), TREATMENT=rep("Bare", times=length(FB.B.samples.s)))
FB.se <- data.frame(SAMPLE=FB.E.samples.s, SITE=rep("FB", times=length(FB.E.samples.s)), TREATMENT=rep("Eelgrass", times=length(FB.E.samples.s)))
samples4anosim.s <- rbind.data.frame(CI.sb, CI.se, PG.sb, PG.se, WB.sb, WB.se, FB.se, FB.sb, stringsAsFactors = TRUE)
samples4anosim.s[,1] <- sub(".s", "", samples4anosim.s[,1])
prtc.batches <- c("A", "C", "A", "A", "C", "A", "A", "C", "A", "A", "A", "C", "C", "C", "A", "C", "C", "C", "C", "C", "A", "D", "C", "B", "D", "A", "B", "A", "A", "A")
data4anosim.s <- cbind.data.frame(SRM.data.mean.stringent.t[order(rownames(SRM.data.mean.stringent.t)),], samples4anosim.s[order(samples4anosim.s$SAMPLE),], PRTC.BATCH=prtc.batches)
str(data4anosim.s)
data4anosim.s$SITE <- as.factor(data4anosim.s$SITE)
data4anosim.s$TREATMENT <- as.factor(data4anosim.s$TREATMENT)
data4anosim.s$PRTC.BATCH <- as.factor(data4anosim.s$PRTC.BATCH)
data4anosim.s[,115]

### start

# A G013|G120|G047|G017|G079|G127|G060|G009|G002|G128|G016|G071-A|G114|G045|G132|G031|G012|G116|G043|G015|G040 
# B G110|G008|G109|G122
# C G041|G066|G105|G032|G129|G054|G081|G003|G074|G014|G049|G053|G104|G055|G042|G064|G073|G057|G007|G070|G001|G071-B|G062
# D G114.remake|G053.remake|G104.remake



### stop

# ANOSIM between sites
sdms.vegdist <- vegdist(data4anosim.s[,-116:-119], 'bray', na.rm=TRUE)
site.anos<-anosim(sdms.vegdist, grouping=data4anosim.s$SITE)
summary(site.anos)
plot(site.anos)

# ANOSIM between treatments
treatment.anos<-anosim(sdms.vegdist, grouping=data4anosim.s$TREATMENT)
summary(treatment.anos)
plot(treatment.anos)

# ANOSIM between PRTC batches
PRTC.anos<-anosim(sdms.vegdist, grouping=data4anosim.s$PRTC.BATCH)
summary(PRTC.anos)
plot(PRTC.anos)



# make figure with sample annotations https://stat.ethz.ch/pipermail/r-sig-ecology/2011-September/002371.html
SRM.nmds.mean.s.log.samples <- scores(SRM.mean.s.log.nmds, display = "sites")
SRM.nmds.mean.s.log.transitions <- scores(SRM.mean.s.log.nmds, display = "species")
# this probably isn't necessary

### Let's plot using ordiplot()

png("Analyses/2017-August_SRM-Analysis/2017-08-24_SRM-NMDS-stringent-plot-log.png")
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
dev.off()

#### Create plot with forced aspect ratio to zoom in ### 

# png("Analyses/2017-August_SRM-Analysis/2017-08-24_SRM-NMDS-stringent-plot-log-zoomed.png")
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

#### PLOT ABUNDANCE FOR EACH PEPTIDE, FOR EACH SITE, TREATMENT ### 
SRM.proteins <- data.frame(protein=SRM.data[2:116,c(1,3,4)]) #protein name to each transition
SRM.proteins[,1] <- sub(" cds.*", "", SRM.proteins[,1])
SRM.data.mean.ordered <- as.matrix(SRM.data.mean[ , order(names(SRM.data.mean))])

CI.B.etc <- cbind(CI.B, Order=c(rep(1, times=nrow(CI.B))), Color=c(rep(2, times=nrow(CI.B))), Symbol=c(rep(20, times=nrow(CI.B))))
CI.E.etc <- cbind(CI.E, Order=c(rep(2, times=nrow(CI.E))), Color=c(rep(2, times=nrow(CI.E))), Symbol=c(rep(100, times=nrow(CI.E))))
PG.B.etc <- cbind(PG.B, Order=c(rep(3, times=nrow(PG.B))), Color=c(rep(6, times=nrow(PG.B))), Symbol=c(rep(20, times=nrow(PG.B))))
PG.E.etc <- cbind(PG.E, Order=c(rep(4, times=nrow(PG.E))), Color=c(rep(6, times=nrow(PG.E))), Symbol=c(rep(100, times=nrow(PG.E))))
WB.B.etc <- cbind(WB.B, Order=c(rep(5, times=nrow(WB.B))), Color=c(rep(3, times=nrow(WB.B))), Symbol=c(rep(20, times=nrow(WB.B))))
WB.E.etc <- cbind(WB.E, Order=c(rep(6, times=nrow(WB.E))), Color=c(rep(3, times=nrow(WB.E))), Symbol=c(rep(100, times=nrow(WB.E))))
FB.B.etc <- cbind(FB.B, Order=c(rep(7, times=nrow(FB.B))), Color=c(rep(4, times=nrow(FB.B))), Symbol=c(rep(20, times=nrow(FB.B))))
FB.E.etc <- cbind(FB.E, Order=c(rep(8, times=nrow(FB.E))), Color=c(rep(4, times=nrow(FB.E))), Symbol=c(rep(100, times=nrow(FB.E))))
samples4plots <- rbind(CI.B.etc, CI.E.etc, PG.B.etc, PG.E.etc, WB.B.etc, WB.E.etc, FB.B.etc, FB.E.etc)
library(plyr)
samples4plots <- arrange(samples4plots, samples4plots$PRVial)

SRM.data4plots <- rbind.data.frame(SRM.data.mean.ordered, Order=samples4plots$Order, Color=samples4plots$Color, Symbol=samples4plots$Symbol)
SRM.data4plots.ordered <- as.matrix(SRM.data4plots[ ,order(SRM.data4plots[which(rownames(SRM.data4plots) == 'Order'), ]) ])

# Write out data files for folks to help me troubleshoot
write.csv(SRM.data4plots.ordered, file="Data/2017-08-31_SRM.data4plots.ordered.csv") #write this file out for safe keeping
write.csv(SRM.proteins, file="Data/2017-08-31_SRM.proteins.csv") #write this file out for safe keeping


# Where to find vectors of plot characteristics in SRM.data4plots.ordered 
#Order = SRM.data4plots.ordered[116,]
#Color = SRM.data4plots.ordered[117,]
# Symbol =SRM.data4plots.ordered[118,]

#### THE FOLLOWING PLOTS NON-STRINGENT MEAN DATA FOR EACH BIOLOGICAL REP. NEED TO 1) ALSO PLOT STRINGENT MEANS, AND 2) FIGURE OUT WHAT THE HECK IS GOING ON WITH MOST OF MY BAR PLOTS. WTF.

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


#### PLOTTING ALL TRANSITION DATA SEPARATELY IN BAR PLOTS #### 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[10,], main=SRM.proteins[10,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[11,], main=SRM.proteins[11,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[12,], main=SRM.proteins[12,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[13,], main=SRM.proteins[13,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[14,], main=SRM.proteins[14,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[15,], main=SRM.proteins[15,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[16,], main=SRM.proteins[16,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[17,], main=SRM.proteins[17,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[18,], main=SRM.proteins[18,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[19,], main=SRM.proteins[19,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[20,], main=SRM.proteins[20,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[21,], main=SRM.proteins[21,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[22,], main=SRM.proteins[22,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[23,], main=SRM.proteins[23,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[24,], main=SRM.proteins[24,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[25,], main=SRM.proteins[25,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[26,], main=SRM.proteins[26,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[27,], main=SRM.proteins[27,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[28,], main=SRM.proteins[28,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[29,], main=SRM.proteins[29,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[30,], main=SRM.proteins[30,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[31,], main=SRM.proteins[31,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[32,], main=SRM.proteins[32,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[33,], main=SRM.proteins[33,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[34,], main=SRM.proteins[34,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[35,], main=SRM.proteins[35,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[36,], main=SRM.proteins[36,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[37,], main=SRM.proteins[37,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[38,], main=SRM.proteins[38,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[39,], main=SRM.proteins[39,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[40,], main=SRM.proteins[40,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[41,], main=SRM.proteins[41,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[42,], main=SRM.proteins[42,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[43,], main=SRM.proteins[43,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[44,], main=SRM.proteins[44,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[45,], main=SRM.proteins[45,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[46,], main=SRM.proteins[46,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[47,], main=SRM.proteins[47,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[48,], main=SRM.proteins[48,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[49,], main=SRM.proteins[49,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[50,], main=SRM.proteins[50,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[51,], main=SRM.proteins[51,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[52,], main=SRM.proteins[52,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[53,], main=SRM.proteins[53,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[54,], main=SRM.proteins[54,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[55,], main=SRM.proteins[55,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[56,], main=SRM.proteins[56,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[57,], main=SRM.proteins[57,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[58,], main=SRM.proteins[58,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[59,], main=SRM.proteins[59,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[60,], main=SRM.proteins[60,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[61,], main=SRM.proteins[61,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[62,], main=SRM.proteins[62,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[63,], main=SRM.proteins[63,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[64,], main=SRM.proteins[64,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[65,], main=SRM.proteins[65,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[66,], main=SRM.proteins[66,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[67,], main=SRM.proteins[67,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[68,], main=SRM.proteins[68,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[69,], main=SRM.proteins[69,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[70,], main=SRM.proteins[70,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[71,], main=SRM.proteins[71,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[72,], main=SRM.proteins[72,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[73,], main=SRM.proteins[73,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[74,], main=SRM.proteins[74,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[75,], main=SRM.proteins[75,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[76,], main=SRM.proteins[76,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[77,], main=SRM.proteins[77,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[78,], main=SRM.proteins[78,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[79,], main=SRM.proteins[79,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[80,], main=SRM.proteins[80,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[83,], main=SRM.proteins[83,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[84,], main=SRM.proteins[84,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[85,], main=SRM.proteins[85,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[86,], main=SRM.proteins[86,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[87,], main=SRM.proteins[87,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[88,], main=SRM.proteins[88,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[81,], main=SRM.proteins[81,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[82,], main=SRM.proteins[82,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[89,], main=SRM.proteins[89,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[90,], main=SRM.proteins[90,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[91,], main=SRM.proteins[91,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[92,], main=SRM.proteins[92,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[93,], main=SRM.proteins[93,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[94,], main=SRM.proteins[94,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[95,], main=SRM.proteins[95,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[96,], main=SRM.proteins[96,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[97,], main=SRM.proteins[97,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[98,], main=SRM.proteins[98,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[99,], main=SRM.proteins[99,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[100,], main=SRM.proteins[100,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[101,], main=SRM.proteins[101,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[102,], main=SRM.proteins[102,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[103,], main=SRM.proteins[103,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[104,], main=SRM.proteins[104,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[105,], main=SRM.proteins[105,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[106,], main=SRM.proteins[106,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# barplot(SRM.data.mean.ordered[107,], main=SRM.proteins[107,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[108,], main=SRM.proteins[108,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[109,], main=SRM.proteins[109,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[110,], main=SRM.proteins[110,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[111,], main=SRM.proteins[111,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[112,], main=SRM.proteins[112,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[113,], main=SRM.proteins[113,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[114,], main=SRM.proteins[114,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# barplot(SRM.data.mean.ordered[115,], main=SRM.proteins[115,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# 
# par(mfrow=c(3,3))
# plot(SRM.data.mean.ordered[107,], main=SRM.proteins[107,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# plot(SRM.data.mean.ordered[108,], main=SRM.proteins[108,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# plot(SRM.data.mean.ordered[109,], main=SRM.proteins[109,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# plot(SRM.data.mean.ordered[110,], main=SRM.proteins[110,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# plot(SRM.data.mean.ordered[111,], main=SRM.proteins[111,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# plot(SRM.data.mean.ordered[112,], main=SRM.proteins[112,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# plot(SRM.data.mean.ordered[113,], main=SRM.proteins[113,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# plot(SRM.data.mean.ordered[114,], main=SRM.proteins[114,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)
# plot(SRM.data.mean.ordered[115,], main=SRM.proteins[115,c(1:3)], col=samples4plots$Color, density=samples4plots$Symbol)

#### BONEYARD AND TROUBLESHOOTING ####

# somehow I have more columns than I should; let's figure out why 
# ncol(SRM.PRTC.a[,-1:-4]) + ncol(SRM.PRTC.b[,-1:-4]) + ncol(SRM.PRTC.c[,-1:-4]) + ncol(SRM.PRTC.d[,-1:-4])
# names <- c(names(SRM.PRTC.a[,-1:-4]), names(SRM.PRTC.b[,-1:-4]), names(SRM.PRTC.c[,-1:-4]), names(SRM.PRTC.d[,-1:-4]))
# length(names) #122
# length(unique(names)) #116; somehow I've picked up 6 extra columns! 
# names[duplicated(names)] #ah, I've included the remade samples in batches a-c, when they should only be in batch d. This is b/c grepl must only include the string

# figuring out why i'm missing data in my NMDS.mean, which won't let me plot bare points
# length(rownames(SRM.nmds.mean.samples)) == length(Eelgrass.samples) + length(Bare.samples) #what's happening?
# names <- c(rownames(SRM.nmds.mean.samples), Bare.samples, Eelgrass.samples)
# length(names)
# length(rownames(SRM.nmds.mean.samples))
# length(Eelgrass.samples) + length(Bare.samples) #what's happening?
# sort(names) 
# sort(rownames(SRM.nmds.mean.samples)) #found that I'm missing G062 in NMDS data


### SUM TRANSITION AREA BY PROTEIN ###
# The "Area" values (which have already been normalized by TIC) are peak area for each transition. Sum them to determine total area for each protein

# NormProtAreaAgg <- aggregate(cbind(FB.E.1, CI.E.1, PG.B.1, SK.E.1, FB.B.1, WB.B.1, SK.B.1, CI.B.1, PG.E.1, WB.B.2, PG.E.2, FB.E.2, FB.B.2, CI.B.2, SK.E.2, PG.B.2, SK.B.2, CI.E.2) ~ Protein.Name, FUN = sum, data = NormProtArea, na.rm = TRUE, na.action = NULL)
# head(NormProtAreaAgg) #confirming that proteins are now summed

### Useful functions from R class notes

# str(SRM.data) # examine a data frame specs
# names() # extract names from data frame
# head(data.frame, n=x) where x=# of rows
# ls() # see what is stored in working directory 

#### testing figures
# ordiplot(SRM.nmds,type="n")
# orditorp(SRM.nmds,display="sites",cex=1.25,air=0.1)

# fig.nmds <- ordiplot(SRM.nmds, choices=c(1,2), 'sites', type='none', xlab='Axis 1', ylab='Axis 2', cex=1)
# ?ordiplot
# points(SRM.nmds, 'sites', col=c("red"), pch=16)
# legend(SRM.nmds)
# NMDS.plot.colors


# My attempt at writing a function with a loop inside it, which calculates the CV for each column in a dataframe (in this example, CV for each PRTC transition)
# CV.vector <- function(SRM.DF) {
#  SRM.DF.t <- t(SRM.DF[,-1:-4])
#  for (i in SRM.DF.t) { 
#    SRM.DF.t[i] <- sd(SRM.DF[,i], na.rm = TRUE)/mean(SRM.PRTC.t[,i], na.rm = TRUE)*100
#  }
#  print(SRM.DF.t)
#}
# CV.vector(SRM.PRTC.a)
