### Note: this script works with datasets downloaded from Skyline.  Skyline data is exported as report, and this script works with data with the following metrics: Protein Name, Transitions, Peptide Sequence, Fragment Ion, Peptide Retention Time, Area

# Script #1 in data processing for TOTAL-SAMPLE-ABUNDANCE-NORMALIZED DATA

############# IMPORT DATASETS #######################################################################################

setwd("~/Documents/Roberts Lab/Geoduck-DNR/") 
SRMreport <- read.csv("Data/2017-08-11_Transition Results_LHS modified-noRT-pivoted.csv", header=FALSE, na.strings = "#N/A", stringsAsFactors = FALSE) # import local file
head(SRMreport)
SRMsequence <- read.csv("Data/2017-07-28_SRM-Sequence-final.csv", header=TRUE, stringsAsFactors = FALSE)
head(SRMsequence)
sample.key <- read.csv("Data/2017-08-14-Geoduck-samples.csv", header=TRUE, stringsAsFactors = FALSE)
head(sample.key)
SRMsamples <- noquote(as.character(c("G013", "G120", "G047", "G017", "G079", "G127", "G060", "G009", "G002", "G128", "G016", "G071-A", "G114", "G045", "G132", "G031", "G012", "G116", "G043", "G015", "G040", "G110", "G008", "G109", "G122", "G041", "G066", "G105", "G032", "G129", "G054", "G081", "G003", "G074", "G014", "G049", "G053", "G104", "G055", "G042", "G064", "G073", "G057", "G007", "G070", "G001", "G071-B", "G062")))

############ REPLACE REP NAMES WITH SAMPLE NAMES ###################################################################
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

############ ANNOTATE SAMPLE NAMES WITH SITE & TREATMENT ##################################################################################

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

############ CONVERT AREA DATA TO NUMERIC FORMAT #######################################################################################

SRM.data.numeric <- SRM.data # First, change all area values to numeric, so I can average, etc. I know that my area data is from column 5 to 120 
SRM.data.numeric[,5:120] <- as.numeric( 
  as.character(
    unlist(
      SRM.data.numeric[,5:120])
  )
)
is.numeric(SRM.data.numeric[2,20]) # confirm area data is numeric, using a random cell. Should equal TRUE.

############ NAME EACH ROW WITH A UNIQUE TRANSITION ID ###############################################################################
nTransitions <- length(SRM.data.numeric$Transition) # How many transitions are there
Transition.ID <- vector(length=nTransitions) # create empty vector with length= number of transitions
for (i in 1:nTransitions) {  
  Transition.ID[i] <- paste(SRM.data.numeric[i,3], SRM.data.numeric[i,4])}  # loop that fills empty vector with unique transition ID, built from the peptide sequence (column 3) and the fragment ion (columm 4)
Transition.ID # confirm correctly named transition IDs
length(SRM.data.numeric$Transition) == length(Transition.ID) # confirm that I didn't lose any transitions
row.names(SRM.data.numeric) <- Transition.ID # assign newly created transition IDs as row names
tail(SRM.data.numeric) # confirm changes
# write.csv(SRM.data.numeric, file="Analyses/2017-September_SRM-results/2017-09-04_SRM-data-TotAbundNORM-annotated.csv") #write this file out for safe keeping

############ NORMALIZE BY TOTAL SAMPLE ABUNDANCE ##################################################################
SRM.data.numeric[is.na(SRM.data.numeric)] <- 0
SRM.data.norm <- ((SRM.data.numeric[-1,-1:-4]/colSums(SRM.data.numeric[-1,-1:-4]))*100)
# write.csv(SRM.data.norm, file="Analyses/2017-September_SRM-results/2017-09-04_SRM-data-TotAbundNORM-data-normalized.csv")
View(SRM.data.norm)

############ CREATE NMDS PLOT ######################################################################################

#Load the source file for the biostats package, biostats.R script must be saved in working directory

source("References/biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
library(vegan)

#Transpose the file so that rows and columns are switched 
SRM.data.t <- t(SRM.data.norm) # t() function transposes, removes PRTC transitions, extraneous info

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
SRM.nmds.pca
#inertia is the sum of all variance in transitions; eigenvalues sum to total inertia, aka each eigenvalue "explains" a certain proportion of the total variance. Percent that each eigenvalue is responsible for total variance is: eigenvalue/total inertia. For example, PC1/total inertia = 83%

summary(SRM.nmds.pca)
plot(SRM.nmds.pca, scaling = 3)
dim(SRM.data.t.noNA)
biplot(SRM.nmds.pca, scaling = -1)
SRM.nmds.ca <- cca(SRM.data.t.noNA)
SRM.nmds.ca


# make figure with sample annotations https://stat.ethz.ch/pipermail/r-sig-ecology/2011-September/002371.html
SRM.nmds.samples <- scores(SRM.nmds, display = "sites")
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
plot.default(x=NULL, y=NULL, type="n", xlab="NMDS axis 1", ylab="NMDS axis 2", xlim=c(-1.5,1.3), ylim=c(-1.5,1), asp=NA)
# symbol key
# 15 = eelgrass = filled square
# 21 = bare = open circle

points(SRM.nmds.samples.sorted[c("G001-A", "G001-B"),], col=colors[1], pch=15)
points(SRM.nmds.samples.sorted[c("G002-A", "G002-B", "G002-C"),], col=colors[50], pch=15) #GOO2-B very diff
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

######## NEED TO TO MODIFY THESE TECH REP AVERAGES BASED ON NMDS

G013 <- ave(SRM.data.numeric$'G013-A', SRM.data.numeric$'G013-C')
G120 <- ave(SRM.data.numeric$`G120-A`, SRM.data.numeric$`G120-B`)
G047 <- ave(SRM.data.numeric$`G047-A`, SRM.data.numeric$`G047-B`)
G017 <- ave(SRM.data.numeric$`G017-A`, SRM.data.numeric$`G017-B`)
G079 <- ave(SRM.data.numeric$`G079-A`, SRM.data.numeric$`G079-B`)
G127 <- ave(SRM.data.numeric$`G127-A`, SRM.data.numeric$`G127-C`) #B removed
G060 <- ave(SRM.data.numeric$`G060-A`, SRM.data.numeric$`G060-B`)
G009 <- ave(SRM.data.numeric$`G009-A`, SRM.data.numeric$`G009-B`)
G002 <- ave(SRM.data.numeric$`G002-B`, SRM.data.numeric$`G002-C`) #A removed
G128 <- ave(SRM.data.numeric$`G128-A`, SRM.data.numeric$`G128-C`,SRM.data.numeric$`G128-D`)
G071.A <- ave(SRM.data.numeric$`G071-A-A`, SRM.data.numeric$`G071-A-B`)
G114 <- ave(SRM.data.numeric$`G114-A`, SRM.data.numeric$`G114-B`, SRM.data.numeric$`G114-remake-C`, SRM.data.numeric$`G114-remake-D`)
G045 <- ave(SRM.data.numeric$`G045-A`, SRM.data.numeric$`G045-B`)
G132 <- ave(SRM.data.numeric$`G132-A`, SRM.data.numeric$`G132-C`, SRM.data.numeric$`G132-D`)
G031 <- ave(SRM.data.numeric$`G031-A`, SRM.data.numeric$`G031-B`, SRM.data.numeric$`G031-C`)
G012 <- ave(SRM.data.numeric$`G012-A`, SRM.data.numeric$`G012-B`, SRM.data.numeric$`G012-C`)
G116 <- ave(SRM.data.numeric$`G116-A`, SRM.data.numeric$`G116-B`)
G043 <- ave(SRM.data.numeric$`G043-A`, SRM.data.numeric$`G043-B`)
G015 <- ave(SRM.data.numeric$`G015-A`, SRM.data.numeric$`G015-B`)
G040 <- ave(SRM.data.numeric$`G040-A`, SRM.data.numeric$`G040-B`)
G110 <- ave(SRM.data.numeric$`G110-A`, SRM.data.numeric$`G110-B`)
G008 <- ave(SRM.data.numeric$`G008-A`, SRM.data.numeric$`G008-B`)
G109 <- ave(SRM.data.numeric$`G109-A`, SRM.data.numeric$`G109-C`)
G122 <- ave(SRM.data.numeric$`G122-A`, SRM.data.numeric$`G122-B`)
G041 <- ave(SRM.data.numeric$`G041-A`, SRM.data.numeric$`G041-B`)
G066 <- ave(SRM.data.numeric$`G066-A`, SRM.data.numeric$`G066-B`)
G105 <- ave(SRM.data.numeric$`G105-A`, SRM.data.numeric$`G105-B`)
G032 <- ave(SRM.data.numeric$`G032-A`, SRM.data.numeric$`G032-B`)
G129 <- ave(SRM.data.numeric$`G129-A`, SRM.data.numeric$`G129-B`)
G054 <- ave(SRM.data.numeric$`G054-A`, SRM.data.numeric$`G054-B`)
G081 <- ave(SRM.data.numeric$`G081-A`, SRM.data.numeric$`G081-B`)
G003 <- ave(SRM.data.numeric$`G003-A`, SRM.data.numeric$`G003-B`) #C removed
G074 <- ave(SRM.data.numeric$`G074-A`, SRM.data.numeric$`G074-B`)
G014 <- ave(SRM.data.numeric$`G014-A`, SRM.data.numeric$`G014-B`)
G049 <- ave(SRM.data.numeric$`G049-A`, SRM.data.numeric$`G049-B`)
G053 <- ave(SRM.data.numeric$`G053-A`, SRM.data.numeric$`G053-remake-C`, SRM.data.numeric$`G053-remake-D`) #B removed 
G104 <- ave(SRM.data.numeric$`G104-A`, SRM.data.numeric$`G104-remake-C`, SRM.data.numeric$`G104-remake-D`) #B removed
G055 <- ave(SRM.data.numeric$`G055-A`, SRM.data.numeric$`G055-B`, SRM.data.numeric$`G055-C`)
G042 <- ave(SRM.data.numeric$`G042-A`, SRM.data.numeric$`G042-B`) #C removed
G064 <- ave(SRM.data.numeric$`G064-A`, SRM.data.numeric$`G064-B`)
G057 <- ave(SRM.data.numeric$`G057-A`, SRM.data.numeric$`G057-C`)
G007 <- ave(SRM.data.numeric$`G007-A`, SRM.data.numeric$`G007-B`)
G070 <- ave(SRM.data.numeric$`G070-A`, SRM.data.numeric$`G070-B`, SRM.data.numeric$`G070-C`)
G001 <- ave(SRM.data.numeric$`G001-A`, SRM.data.numeric$`G001-B`)
G071.B <- ave(SRM.data.numeric$`G071-B-A`, SRM.data.numeric$`G071-B-B`)
G062 <- ave(SRM.data.numeric$`G062-B`, SRM.data.numeric$`G062-C`)
# removed G016, G073

SRM.data.mean <- cbind.data.frame(rownames(SRM.data.numeric), G013, G120, G047, G017, G079, G127, G060, G009, G002, G128, G071.A, G114, G045, G132, G031, G012, G116, G043, G015, G040, G110, G008, G109, G122, G041, G066, G105, G032, G129, G054, G081, G003, G074, G014, G049, G053, G104, G055, G042, G064, G057, G007, G070,  G001, G071.B, G062) # combine all tech. replicate mean vectors into new data frame 
# write.csv(SRM.data.mean, file="Analyses/Analyses/2017-September_SRM-results/2017-09-04_SRM-data-meanBYsample.csv")
