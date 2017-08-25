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
CI.E <- repsTOsamples.filtered.annotated[c(repsTOsamples.filtered.annotated$Sample.Shorthand == "CI-E"),]
CI.B <- repsTOsamples.filtered.annotated[c(repsTOsamples.filtered.annotated$Sample.Shorthand == "CI-B"),]
PG.E <- repsTOsamples.filtered.annotated[c(repsTOsamples.filtered.annotated$Sample.Shorthand == "PG-E"),]
PG.B <- repsTOsamples.filtered.annotated[c(repsTOsamples.filtered.annotated$Sample.Shorthand == "PG-B"),]
WB.E <- repsTOsamples.filtered.annotated[c(repsTOsamples.filtered.annotated$Sample.Shorthand == "WB-E"),]
WB.B <- repsTOsamples.filtered.annotated[c(repsTOsamples.filtered.annotated$Sample.Shorthand == "WB-B"),]
FB.E <- repsTOsamples.filtered.annotated[c(repsTOsamples.filtered.annotated$Sample.Shorthand == "FB-E"),]
FB.B <- repsTOsamples.filtered.annotated[c(repsTOsamples.filtered.annotated$Sample.Shorthand == "FB-B"),]

CI.E.samples <- CI.E$PRVial
CI.B.samples <- CI.B$PRVial
PG.E.samples <- PG.E$PRVial
PG.B.samples <- PG.B$PRVial
WB.E.samples <- WB.E$PRVial
WB.B.samples <- WB.B$PRVial
FB.E.samples <- FB.E$PRVial
FB.B.samples <- FB.B$PRVial

### CONVERT AREA DATA TO NUMERIC FORMAT ### 

SRM.data.numeric <- SRM.data # First, change all area values to numeric, so I can average, etc. I know that my area data is from column 5 to 120 
SRM.data.numeric[,5:120] <- as.numeric( 
  as.character(
    unlist(
      SRM.data.numeric[,5:120])
  )
)
is.numeric(SRM.data.numeric[2,20]) # confirm area data is numeric, using a random cell. Should equal TRUE.

# Name each row with a unique transition ID
nTransitions <- length(SRM.data.numeric$Transition) # How many transitions are there
Transition.ID <- vector(length=nTransitions) # create empty vector with length= number of transitions
for (i in 1:nTransitions) {  
  Transition.ID[i] <- paste(SRM.data.numeric[i,3], SRM.data.numeric[i,4])}  # loop that fills empty vector with unique transition ID, built from the peptide sequence (column 3) and the fragment ion (columm 4)
Transition.ID # confirm correctly named transition IDs
length(SRM.data.numeric$Transition) == length(Transition.ID) # confirm that I didn't lose any transitions
row.names(SRM.data.numeric) <- Transition.ID # assign newly created transition IDs as row names
head(SRM.data.numeric) # confirm changes
# write.csv(SRM.data.numeric, file="Data/2017-08-19_SRM-data-R1.csv") #write this file out for safe keeping

### NORMALIZE ### 
# Discard PRTC peptides that elute early and late for this analysis since their quantities are less stable. These were identified via Skyline, and are the following: 
### DID I ALREADY DO THIS? NEED TO CHECK IT OUT ON MY LAB NOTEBOOK ### 

# Normalize peak area abundances by averaged PRTC peptide abundances for each sample to normalize for protein loading. PRTC concentrations varied between sample prep; I factor this in by normalizing again by X, Y, Z... [TBD]

SRM.data.numeric[115:142,1:5] #check out which rows pertain to PRTC peptides, it's 117->142
SRM.PRTC <- SRM.data.numeric[116:142,] #pull out PRTC data for each sample

# Adjust PRTC abundance for remake samples, since that was made at 10% the concentration of the others PRTC batches

SRM.PRTC.adjusted <- SRM.PRTC
SRM.PRTC.adjusted
SRM.PRTC.adjusted$`G114-remake-C` <- SRM.PRTC.adjusted$`G114-remake-C`*10
SRM.PRTC.adjusted$`G053-remake-C` <- SRM.PRTC.adjusted$`G053-remake-C`*10
SRM.PRTC.adjusted$`G104-remake-C` <- SRM.PRTC.adjusted$`G104-remake-C`*10
SRM.PRTC.adjusted$`G104-remake-D` <- SRM.PRTC.adjusted$`G104-remake-D`*10
SRM.PRTC.adjusted$`G053-remake-D` <- SRM.PRTC.adjusted$`G053-remake-D`*10
SRM.PRTC.adjusted$`G114-remake-D` <- SRM.PRTC.adjusted$`G114-remake-D`*10

# calculate mean abundance fo all PRTC transitions within samples
SRM.PRTC.adjusted.mean <- data.frame(colMeans(SRM.PRTC.adjusted[,-1:-4], na.rm=TRUE))
SRM.PRTC.adjusted.mean

# Now, normalize all sample abundance data based on PRTC mean abundances
head(SRM.data.numeric.1)
SRM.data.numeric.1 <- SRM.data.numeric[c(-1,-116:-142),-1:-4] #remove non-data related columns, remove PRTC transitions
PRTC.norm.vector <- SRM.PRTC.adjusted.mean[,1] #create vector of mean PRTC abundances for each sample
length(PRTC.norm.vector) == ncol(SRM.data.numeric.1) #confirm PRTC normalization vector length equals # samples in srm data
SRM.data.normalized <- sweep(SRM.data.numeric.1, 2, PRTC.norm.vector, "/") #normalize srm data (averaged by tech. rep) by mean PRTC abundance for that sample
head(SRM.data.normalized)

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
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-2,3), ylim=c(-0.5,0.5), asp=NA)
# symbol key
# 15 = eelgrass = filled square
# 21 = bare = open circle

points(SRM.nmds.samples.sorted[c("G001-A", "G001-B"),], col=colors[1], pch=15)
points(SRM.nmds.samples.sorted[c("G002-A", "G002-B", "G002-C"),], col=colors[50], pch=15)
points(SRM.nmds.samples.sorted[c("G003-A", "G003-B", "G003-C"),], col=colors[3], pch=15) #G003-C is very different
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
points(SRM.nmds.samples.sorted[c("G053-A", "G053-B", "G053-remake-C", "G053-remake-D"),], col=colors[22], pch=15) #remakes look good, A&B do not
points(SRM.nmds.samples.sorted[c("G054-A", "G054-B"),], col=colors[23], pch=15)
points(SRM.nmds.samples.sorted[c("G055-A", "G055-B"),], col=colors[24], pch=15) #one is very off
points(SRM.nmds.samples.sorted[c("G057-A", "G057-B", "G057-C"),], col=colors[25], pch=21) #one is off
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
points(SRM.nmds.samples.sorted[c("G104-A", "G104-B", "G104-remake-C", "G104-remake-D"),], col=colors[37], pch=21) #check these out too
points(SRM.nmds.samples.sorted[c("G105-A", "G105-B"),], col=colors[38], pch=21)
points(SRM.nmds.samples.sorted[c("G109-A", "G109-B", "G109-C"),], col=colors[39], pch=15)
points(SRM.nmds.samples.sorted[c("G110-A", "G110-B"),], col=colors[40], pch=15)
points(SRM.nmds.samples.sorted[c("G114-A", "G114-B", "G114-remake-C", "G114-remake-D"),], col=colors[41], pch=21) # check these out!
points(SRM.nmds.samples.sorted[c("G116-A", "G116-B"),], col=colors[42], pch=21)
points(SRM.nmds.samples.sorted[c("G120-A", "G120-B"),], col=colors[43], pch=21)
points(SRM.nmds.samples.sorted[c("G122-A", "G122-B"),], col=colors[44], pch=21)
points(SRM.nmds.samples.sorted[c("G127-A", "G127-B", "G127-C"),], col=colors[45], pch=15) #one is off
points(SRM.nmds.samples.sorted[c("G128-A", "G128-C", "G128-D"),], col=colors[46], pch=15)
points(SRM.nmds.samples.sorted[c("G129-A", "G129-B"),], col=colors[47], pch=15)
points(SRM.nmds.samples.sorted[c("G132-A", "G132-C", "G132-D"),], col=colors[48], pch=15)

#### NEXT, REMOVE SAMPLES THAT DON'T LOOK GOOD, AVERAGE TECH REPS, THEN RE-PLOT BY SITE/TREATMENT #### 

### Plotting by site and treatment ####
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-2,3), ylim=c(-0.5,0.5), asp=NA)
# symbol key
# 15 = eelgrass = filled square
# 21 = bare = open circle
# red = PG
# blue = CI
# green = WB
# black = FB

CI.E.samples
CI.B.samples
PG.E.samples
PG.B.samples
WB.E.samples
WB.B.samples
FB.E.samples
FB.B.samples

#### testing figures
ordiplot(SRM.nmds,type="n")
orditorp(SRM.nmds,display="sites",cex=1.25,air=0.1)

fig.nmds <- ordiplot(SRM.nmds, choices=c(1,2), 'sites', type='none', xlab='Axis 1', ylab='Axis 2', cex=1)
?ordiplot
points(SRM.nmds, 'sites', col=c("red"), pch=16)
legend(SRM.nmds)
NMDS.plot.colors

# average sample technical reps.  (there's probably an easier way to do this to not manually enter the tech rep names for each sample, possibly via a loop?)

G013 <- ave(SRM.data.normalized$'G013-A', SRM.data.normalized$'G013-C')
G120 <- ave(SRM.data.normalized$`G120-A`, SRM.data.normalized$`G120-B`)
G047 <- ave(SRM.data.normalized$`G047-A`, SRM.data.normalized$`G047-B`)
G017 <- ave(SRM.data.normalized$`G017-A`, SRM.data.normalized$`G017-B`)
G079 <- ave(SRM.data.normalized$`G079-A`, SRM.data.normalized$`G079-B`)
G127 <- ave(SRM.data.normalized$`G127-A`, SRM.data.normalized$`G127-B`, SRM.data.normalized$`G127-C`)
G060 <- ave(SRM.data.normalized$`G060-A`, SRM.data.normalized$`G060-B`)
G009 <- ave(SRM.data.normalized$`G009-A`, SRM.data.normalized$`G009-B`)
G002 <- ave(SRM.data.normalized$`G002-A`, SRM.data.normalized$`G002-B`, SRM.data.normalized$`G002-C`)
G128 <- ave(SRM.data.normalized$`G128-A`, SRM.data.normalized$`G128-C`,SRM.data.normalized$`G128-D`)
G016 <- ave(SRM.data.normalized$`G016-A`, SRM.data.normalized$`G016-B`, SRM.data.normalized$`G016-C`)
G071.A <- ave(SRM.data.normalized$`G071-A-A`, SRM.data.normalized$`G071-A-B`)
G114 <- ave(SRM.data.normalized$`G114-A`, SRM.data.normalized$`G114-B`)
G114.remake <- ave(SRM.data.normalized$`G114-remake-C`, SRM.data.normalized$`G114-remake-D`)
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
G003 <- ave(SRM.data.normalized$`G003-A`, SRM.data.normalized$`G003-B`, SRM.data.normalized$`G003-C`)
G074 <- ave(SRM.data.normalized$`G074-A`, SRM.data.normalized$`G074-B`)
G014 <- ave(SRM.data.normalized$`G014-A`, SRM.data.normalized$`G014-B`)
G049 <- ave(SRM.data.normalized$`G049-A`, SRM.data.normalized$`G049-B`)
G053 <- ave(SRM.data.normalized$`G053-A`, SRM.data.normalized$`G053-B`) 
G053.remake <- ave(SRM.data.normalized$`G053-remake-C`, SRM.data.normalized$`G053-remake-D`)
G104 <- ave(SRM.data.normalized$`G104-A`, SRM.data.normalized$`G104-B`)
G104.remake <- ave(SRM.data.normalized$`G104-remake-C`, SRM.data.normalized$`G104-remake-D`)
G055 <- ave(SRM.data.normalized$`G055-A`, SRM.data.normalized$`G055-B`, SRM.data.normalized$`G055-C`)
G042 <- ave(SRM.data.normalized$`G042-A`, SRM.data.normalized$`G042-B`, SRM.data.normalized$`G042-C`)
G064 <- ave(SRM.data.normalized$`G064-A`, SRM.data.normalized$`G064-B`)
G073 <- ave(SRM.data.normalized$`G073-A`, SRM.data.normalized$`G073-B`, SRM.data.normalized$`G073-C`)
G057 <- ave(SRM.data.normalized$`G057-A`, SRM.data.normalized$`G057-B`, SRM.data.normalized$`G057-C`)
G007 <- ave(SRM.data.normalized$`G007-A`, SRM.data.normalized$`G007-B`)
G070 <- ave(SRM.data.normalized$`G070-A`, SRM.data.normalized$`G070-B`, SRM.data.normalized$`G070-C`)
G001 <- ave(SRM.data.normalized$`G001-A`, SRM.data.normalized$`G001-B`)
G071.B <- ave(SRM.data.normalized$`G071-B-A`, SRM.data.normalized$`G071-B-B`)

rownames(SRM.data.normalized)
SRM.data.mean <- noquote(cbind(rownames(SRM.data.normalized), G013, G120, G047, G017, G079, G127, G060, G009, G002, G128, G016, G071.A, G114, G114.remake, G045, G132, G031, G012, G116, G043, G015, G040, G110, G008, G109, G122, G041, G066, G105, G032, G129, G054, G081, G003, G074, G014, G049, G053, G053.remake, G104, G104.remake, G055, G042, G064, G073, G057, G007, G070,  G001, G071.B)) # combine all tech. replicate mean vectors into new data frame 
SRM.data.mean
