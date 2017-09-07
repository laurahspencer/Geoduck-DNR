# Script #1 for PRTC-Normalized SRM Data

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

### ANNOTATE SAMPLE NAMES WITH SITE & TREATMENT ########## 

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
# png("Analyses/2017-August_SRM-Analysis/PRTC-pre-post-adjusted-1.png", width = 1000, height = 1200, units = "px")
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

# png("Analyses/2017-August_SRM-Analysis/PRTC-pre-post-adjusted-2.png", width = 1000, height = 1200, units = "px")
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
View(PRTC.cv.comparison) #this is the result

### NORMALIZE ASSAY DATA BASED ON PRTC ABUNDANCE 

# calculate mean abundance fo all PRTC transitions within samples
SRM.PRTC.adjusted.mean <- data.frame(colMeans(SRM.PRTC.adjusted, na.rm=TRUE))

# Now, normalize all sample abundance data based on PRTC mean abundances
SRM.data.numeric.1 <- SRM.data.numeric[c(-1,-117:-142),-1:-4] #remove non-abundance-data  columns, remove PRTC transitions from assay data
PRTC.norm.vector <- SRM.PRTC.adjusted.mean[,1] #create vector of mean PRTC abundances for each sample
SRM.PRTC.adjusted.mean
length(PRTC.norm.vector) == ncol(SRM.data.numeric.1) #confirm PRTC normalization vector length equals # samples in srm data
SRM.data.normalized <- sweep(SRM.data.numeric.1, 2, PRTC.norm.vector, "/") #normalize srm data (averaged by tech. rep) by mean PRTC abundance for that sample
# write.csv(SRM.data.normalized, file="Analyses/2017-September_SRM-results/2017-09-24_PRTCNorm-dataBYtechrep-annotaed.csv)
