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
write.csv(SRM.data.numeric, file="Data/2017-08-19_SRM-data-R1.csv") #write this file out for safe keeping

### NORMALIZE ### 
# Discard PRTC peptides that elute early and late for this analysis since their quantities are less stable. These were identified via Skyline, and are the following: 
### DID I ALREADY DO THIS? NEED TO CHECK IT OUT ON MY LAB NOTEBOOK ### 

# Normalize peak area abundances by averaged PRTC peptide abundances for each sample to normalize for protein loading. PRTC concentrations varied between sample prep; I factor this in by normalizing again by X, Y, Z... [TBD]

SRM.data.numeric[115:142,1:5] #check out which rows pertain to PRTC peptides
PRTC.mean <- 1:ncol(SRM.data.numeric) # Create vector to be filled by loop with mean PRTC area
length(PRTC.mean) == ncol(SRM.data.numeric) #confirm vector has 
for (i in PRTC.mean) { 
  if (i <= 5) {   # skip first 5 iterations since they correspond to non-numeric/area info.
    next
  }
  PRTC.mean[i] <- mean(SRM.data.numeric[117:142,i], na.rm=TRUE)  # create a sample, each entry represents a sample with mean abundance for all PRTC transitions in that sample
}
print(PRTC.mean)

# average sample technical reps.  (there's probably an easier way to do this to not manually enter the tech rep names for each sample, possibly via a loop?)

G013 <- ave(SRM.data.numeric$'G013-A', SRM.data.numeric$'G013-C')
G120 <- ave(SRM.data.numeric$`G120-A`, SRM.data.numeric$`G120-B`)
G047 <- ave(SRM.data.numeric$`G047-A`, SRM.data.numeric$`G047-B`)
G017 <- ave(SRM.data.numeric$`G017-A`, SRM.data.numeric$`G017-B`)
G079 <- ave(SRM.data.numeric$`G079-A`, SRM.data.numeric$`G079-B`)
G127 <- ave(SRM.data.numeric$`G127-A`, SRM.data.numeric$`G127-B`, SRM.data.numeric$`G127-C`)
G060 <- ave(SRM.data.numeric$`G060-A`, SRM.data.numeric$`G060-B`)
G009 <- ave(SRM.data.numeric$`G009-A`, SRM.data.numeric$`G009-B`)
G002 <- ave(SRM.data.numeric$`G002-A`, SRM.data.numeric$`G002-B`, SRM.data.numeric$`G002-C`)
G128 <- ave(SRM.data.numeric$`G128-A`, SRM.data.numeric$`G128-C`,SRM.data.numeric$`G128-D`)
G016 <- ave(SRM.data.numeric$`G016-A`, SRM.data.numeric$`G016-B`, SRM.data.numeric$`G016-C`)
G071.A <- ave(SRM.data.numeric$`G071-A-A`, SRM.data.numeric$`G071-A-B`)
G114 <- ave(SRM.data.numeric$`G114-A`, SRM.data.numeric$`G114-B`)
G114.remake <- ave(SRM.data.numeric$`G114-remake-C`, SRM.data.numeric$`G114-remake-D`)
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
G109 <- ave(SRM.data.numeric$`G109-A`, SRM.data.numeric$`G109-B`, SRM.data$`G109-C`)
G122 <- ave(SRM.data.numeric$`G122-A`, SRM.data.numeric$`G122-B`)
G041 <- ave(SRM.data.numeric$`G041-A`, SRM.data.numeric$`G041-B`)
G066 <- ave(SRM.data.numeric$`G066-A`, SRM.data.numeric$`G066-B`)
G105 <- ave(SRM.data.numeric$`G105-A`, SRM.data.numeric$`G105-B`)
G032 <- ave(SRM.data.numeric$`G032-A`, SRM.data.numeric$`G032-B`)
G129 <- ave(SRM.data.numeric$`G129-A`, SRM.data.numeric$`G129-B`)
G054 <- ave(SRM.data.numeric$`G054-A`, SRM.data.numeric$`G054-B`)
G081 <- ave(SRM.data.numeric$`G081-A`, SRM.data.numeric$`G081-B`)
G003 <- ave(SRM.data.numeric$`G003-A`, SRM.data.numeric$`G003-B`, SRM.data.numeric$`G003-C`)
G074 <- ave(SRM.data.numeric$`G074-A`, SRM.data.numeric$`G074-B`)
G014 <- ave(SRM.data.numeric$`G014-A`, SRM.data.numeric$`G014-B`)
G049 <- ave(SRM.data.numeric$`G049-A`, SRM.data.numeric$`G049-B`)
G053 <- ave(SRM.data.numeric$`G053-A`, SRM.data.numeric$`G053-B`) 
G053.remake <- ave(SRM.data.numeric$`G053-remake-C`, SRM.data.numeric$`G053-remake-D`)
G104 <- ave(SRM.data.numeric$`G104-A`, SRM.data.numeric$`G104-B`)
G104.remake <- ave(SRM.data.numeric$`G104-remake-C`, SRM.data.numeric$`G104-remake-D`)
G055 <- ave(SRM.data.numeric$`G055-A`, SRM.data.numeric$`G055-B`, SRM.data.numeric$`G055-C`)
G042 <- ave(SRM.data.numeric$`G042-A`, SRM.data.numeric$`G042-B`, SRM.data.numeric$`G042-C`)
G064 <- ave(SRM.data.numeric$`G064-A`, SRM.data.numeric$`G064-B`)
G073 <- ave(SRM.data.numeric$`G073-A`, SRM.data.numeric$`G073-B`, SRM.data.numeric$`G073-C`)
G057 <- ave(SRM.data.numeric$`G057-A`, SRM.data.numeric$`G057-B`, SRM.data.numeric$`G057-C`)
G007 <- ave(SRM.data.numeric$`G007-A`, SRM.data.numeric$`G007-B`)
G070 <- ave(SRM.data.numeric$`G070-A`, SRM.data.numeric$`G070-B`, SRM.data.numeric$`G070-C`)
G001 <- ave(SRM.data.numeric$`G001-A`, SRM.data.numeric$`G001-B`)
G071.B <- ave(SRM.data.numeric$`G071-B-A`, SRM.data.numeric$`G071-B-B`)

head(SRM.data.numeric) 
SRM.data.mean <- cbind(SRM.data.numeric[,1:4], G013, G120, G047, G017, G079, G127, G060, G009, G002, G128, G016, G071.A, G114, G114.remake, G045, G132, G031, G012, G116, G043, G015, G040, G110, G008, G109, G122, G041, G066, G105, G032, G129, G054, G081, G003, G074, G014, G049, G053, G053.remake, G104, G104.remake, G055, G042, G064, G073, G057, G007, G070,  G001, G071.B) # combine all tech. replicate mean vectors into new data frame 
SRM.data.mean <- SRM.data.mean[-1,] # remove first row, since it is full of NA's (it contains column names)
SRM.PRTC <- SRM.data.mean[116:142,] # 

# First, try adjusting PRTC abundance for just remake samples, since that was made at 10% the concentration of the others PRTC batches

SRM.PRTC.adjusted <- SRM.PRTC
SRM.PRTC.adjusted$G114.remake <- SRM.PRTC.adjusted$G114.remake*10
SRM.PRTC.adjusted$G053.remake <- SRM.PRTC.adjusted$G053.remake*10
SRM.PRTC.adjusted$G104.remake <- SRM.PRTC.adjusted$G104.remake*10

# calculate mean abundance fo all PRTC transitions within samples
SRM.PRTC.adjusted.mean <- data.frame(colMeans(SRM.PRTC.adjusted[,-1:-4], na.rm=TRUE))

# Now, normalize all sample abundance data based on PRTC mean abundances
SRM.data.mean.1 <- SRM.data.mean[-116:-142,-1:-4] #remove non-data related columns, remove PRTC transitions
PRTC.norm.vector <- SRM.PRTC.adjusted.mean[,1] #create vector of mean PRTC abundances for each sample
SRM.data.mean.normalized <- sweep(SRM.data.mean.1, 2, PRTC.norm.vector, "/") #normalize srm data (averaged by tech. rep) by mean PRTC abundance for that sample

#### CREATE NMDS PLOT ########

#Load the source file for the biostats package, biostats.R script must be saved in working directory

source("References/biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
library(vegan)

#Transpose the file so that rows and columns are switched 
SRM.data.t <- t(SRM.data.mean.normalized) # t() function transposes

#Replace NA cells with 0; metaMDS() does not handle NA's
SRM.data.t.noNA <- SRM.data.t
SRM.data.t.noNA[is.na(SRM.data.t.noNA)] <- 0
head(SRM.data.t.noNA)

#Make MDS dissimilarity matrix
#
SRM.nmds <- metaMDS(SRM.data.t.noNA, distance = 'bray', k = 2, trymax = 1000, autotransform = FALSE)
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

#### testing figures
ordiplot(SRM.nmds,type="n")
orditorp(SRM.nmds,display="sites",cex=1.25,air=0.01)

fig.nmds <- ordiplot(SRM.nmds, choices=c(1,2), 'sites', type='none', xlab='Axis 1', ylab='Axis 2', cex=1)
?ordiplot
points(SRM.nmds, 'sites', col=c("red", "blue"), pch=c(rep(16,5), rep(17,5)))
legend(SRM.nmds)
NMDS.plot.colors

##### THIS NEEDS TO BE UPDATED ###### 
#bare=circle
#eelgrass=triangle
#case=red
#fidalgo=blue
#willapa=black
#skokomish=green
#gamble=magenta

?points
#how manny points do I have?
nrow(SRM.data.t.noNA) 
points(fig.nmds, 'species')


NMDS.plot.colors <- 1:nrow(SRM.data.t.noNA) # create empty vector with length= number of samples
print(NMDS.plot.colors)

legend(0.185,0.1, pch=c(rep(16,5), 1, 2), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'))

# Notes on stress level from https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/ - A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions, < 0.1 is great, < 0.2 is good/ok, and stress < 0.3 provides a poor representation.** To reiterate: high stress is bad, low stress is good! NOTE: metaMDS has automatically applied a square root transformation and calculated the Bray-Curtis distances for our community-by-site matrix.

### ANNOTATE SAMPLE NAMES WITH SITE & TREATMENT ########## needs work

head(sample.key)
sample.key[,c(8,9)]
repsTOsamples.filtered.annotated <- filter(sample.key[,c(8,9)], sample.key$PRVial %in% repsTOsamples.filtered$Comment) #pull site & treatment from sample key
length(SRMsamples) # check # samples there should be total
nrow(repsTOsamples.filtered.annotated) # check # samples is appropriate
repsTOsamples.filtered.annotated # NOTE: missing 71-A & 71-B, need to fix 


#### AVERAGE REPLICATES #### 
# there's probably an easier way to do this to not manually enter the tech rep names for each sample, possibly via a loop?

G013 <- ave(SRM.data.numeric$'G013-A', SRM.data.numeric$'G013-C')
G120 <- ave(SRM.data.numeric$`G120-A`, SRM.data.numeric$`G120-B`)
G047 <- ave(SRM.data.numeric$`G047-A`, SRM.data.numeric$`G047-B`)
G017 <- ave(SRM.data.numeric$`G017-A`, SRM.data.numeric$`G017-B`)
G079 <- ave(SRM.data.numeric$`G079-A`, SRM.data.numeric$`G079-B`)
G127 <- ave(SRM.data.numeric$`G127-A`, SRM.data.numeric$`G127-B`, SRM.data.numeric$`G127-C`)
G060 <- ave(SRM.data.numeric$`G060-A`, SRM.data.numeric$`G060-B`)
G009 <- ave(SRM.data.numeric$`G009-A`, SRM.data.numeric$`G009-B`)
G002 <- ave(SRM.data.numeric$`G002-A`, SRM.data.numeric$`G002-B`, SRM.data.numeric$`G002-C`)
G128 <- ave(SRM.data.numeric$`G128-A`, SRM.data.numeric$`G128-C`,SRM.data.numeric$`G128-D`)
G016 <- ave(SRM.data.numeric$`G016-A`, SRM.data.numeric$`G016-B`, SRM.data.numeric$`G016-C`)
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
G109 <- ave(SRM.data.numeric$`G109-A`, SRM.data.numeric$`G109-B`, SRM.data$`G109-C`)
G122 <- ave(SRM.data.numeric$`G122-A`, SRM.data.numeric$`G122-B`)
G041 <- ave(SRM.data.numeric$`G041-A`, SRM.data.numeric$`G041-B`)
G066 <- ave(SRM.data.numeric$`G066-A`, SRM.data.numeric$`G066-B`)
G105 <- ave(SRM.data.numeric$`G105-A`, SRM.data.numeric$`G105-B`)
G032 <- ave(SRM.data.numeric$`G032-A`, SRM.data.numeric$`G032-B`)
G129 <- ave(SRM.data.numeric$`G129-A`, SRM.data.numeric$`G129-B`)
G054 <- ave(SRM.data.numeric$`G054-A`, SRM.data.numeric$`G054-B`)
G081 <- ave(SRM.data.numeric$`G081-A`, SRM.data.numeric$`G081-B`)
G003 <- ave(SRM.data.numeric$`G003-A`, SRM.data.numeric$`G003-B`, SRM.data.numeric$`G003-C`)
G074 <- ave(SRM.data.numeric$`G074-A`, SRM.data.numeric$`G074-B`)
G014 <- ave(SRM.data.numeric$`G014-A`, SRM.data.numeric$`G014-B`)
G049 <- ave(SRM.data.numeric$`G049-A`, SRM.data.numeric$`G049-B`)
G053 <- ave(SRM.data.numeric$`G053-A`, SRM.data.numeric$`G053-B`, SRM.data.numeric$`G053-remake-C`, SRM.data.numeric$`G053-remake-D`)
G104 <- ave(SRM.data.numeric$`G104-A`, SRM.data.numeric$`G104-B`, SRM.data.numeric$`G104-remake-C`, SRM.data.numeric$`G104-remake-D`)
G055 <- ave(SRM.data.numeric$`G055-A`, SRM.data.numeric$`G055-B`, SRM.data.numeric$`G055-C`)
G042 <- ave(SRM.data.numeric$`G042-A`, SRM.data.numeric$`G042-B`, SRM.data.numeric$`G042-C`)
G064 <- ave(SRM.data.numeric$`G064-A`, SRM.data.numeric$`G064-B`)
G073 <- ave(SRM.data.numeric$`G073-A`, SRM.data.numeric$`G073-B`, SRM.data.numeric$`G073-C`)
G057 <- ave(SRM.data.numeric$`G057-A`, SRM.data.numeric$`G057-B`, SRM.data.numeric$`G057-C`)
G007 <- ave(SRM.data.numeric$`G007-A`, SRM.data.numeric$`G007-B`)
G070 <- ave(SRM.data.numeric$`G070-A`, SRM.data.numeric$`G070-B`, SRM.data.numeric$`G070-C`)
G001 <- ave(SRM.data.numeric$`G001-A`, SRM.data.numeric$`G001-B`)
G071.B <- ave(SRM.data.numeric$`G071-B-A`, SRM.data.numeric$`G071-B-B`)

SRM.data.mean <- cbind(SRM.data[,1:4], G013, G120, G047, G017, G079, G127, G060, G009, G002, G128, G016, G071.A, G114, G045, G132, G031, G012, G116, G043, G015, G040, G110, G008, G109, G122, G041, G066, G105, G032, G129, G054, G081, G003, G074, G014, G049, G053, G104, G055, G042, G064, G073, G057, G007, G070,  G001, G071.B)
head(SRM.data.mean, n=3) # SRM.data.mean is a data frame with with area averaged over all technical replicates for each sample. 



###################### i'm thinking that this might be necessary; TBD ###################
###### Assign PRTC batches based on which PRTC mix samples were made with ######

SRM.PRTC.a <- cbind(
  SRM.PRTC[,1:4], 
  SRM.PRTC[, grepl("G013|G120|G047|G017|G079|G127|G060|G009|G002|G128|G016|G071-A|G114|G045|G132|G031|G012|G116|G043|G015|G040", names(SRM.PRTC))]
)

SRM.PRTC.b <- cbind(
  SRM.PRTC[,1:4],
  SRM.PRTC[, grepl("G110|G008|G109|G122", names(SRM.PRTC))]
)

SRM.PRTC.c <- cbind(
  SRM.PRTC[,1:4],
  SRM.PRTC[, grepl("G041|G066|G105|G032|G129|G054|G081|G003|G074|G014|G049|G053|G104|G055|G042|G064|G073|G057|G007|G070|G001|G071-B|G062", names(SRM.PRTC))]
)

SRM.PRTC.d <- cbind(
  SRM.PRTC[,1:4],
  SRM.PRTC[, grepl("G114.remake|G053.remake|G104.remake", names(SRM.PRTC))]
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
PRTC.batch.means <- cbind(SRM.PRTC[,1:4], SRM.PRTC.a.mean, SRM.PRTC.b.mean, SRM.PRTC.c.mean, SRM.PRTC.d.mean)
write.csv(PRTC.batch.means, file="Data/2017-08-22_SRM-PRTC-batch-means.csv") #write this file out for safe keeping

# divide batches a, b & c by d, since d should have ~10% the PRTC compared to the others
PRTC.a.ratio <- mean(PRTC.batch.means$SRM.PRTC.a.mean/PRTC.batch.means$SRM.PRTC.d.mean, na.rm = TRUE )
PRTC.b.ratio <- mean(PRTC.batch.means$SRM.PRTC.b.mean/PRTC.batch.means$SRM.PRTC.d.mean, na.rm = TRUE )
PRTC.c.ratio <- mean(PRTC.batch.means$SRM.PRTC.c.mean/PRTC.batch.means$SRM.PRTC.d.mean, na.rm = TRUE )
PRTC.d.ratio <- mean(PRTC.batch.means$SRM.PRTC.d.mean/PRTC.batch.means$SRM.PRTC.d.mean, na.rm = TRUE )

PRTC.a.ratio
PRTC.b.ratio
PRTC.c.ratio
PRTC.d.ratio

######### plot PRTC data using NMDS ##########

#Load the source file for the biostats package, biostats.R script must be saved in working directory

source("References/biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
library(vegan)

# Extract PRTC transitions 
### NEED TO USE NORMALIZED DATA SET INSTEAD, THEN DO NMDS WITH NORMALIZED AND PRTC REMOVED DATA SET ### 
SRM.data.numeric[115:142,1:5] #check out which rows pertain to PRTC peptides
SRM.data.PRTC <- SRM.data.numeric[c(117:142),]

# Name each row with a unique transition ID
nTransitions.PRTC <- length(SRM.data.PRTC$Transition) # How many transitions are there
Transition.ID.PRTC <- vector(length=nTransitions) # create empty vector with length= number of transitions
for (i in 1:nTransitions.PRTC) {  
  Transition.ID.PRTC[i] <- paste(SRM.data.PRTC[i,3], SRM.data.PRTC[i,4])}  # loop that fills empty vector with unique transition ID, built from the peptide sequence (column 3) and the fragment ion (columm 4)
Transition.ID.PRTC # confirm correctly named transition IDs
row.names(SRM.data.PRTC) <- Transition.ID.PRTC # assign newly created transition IDs as row names
head(SRM.data.PRTC) # confirm changes
write.csv(SRM.data.PRTC, file="Data/2017-08-23_SRM-PRTC.csv") #write this file out for safe keeping

#Transpose the file so that rows and columns are switched 
SRM.data.PRTC.t <- t(SRM.data.PRTC[,-1:-4]) # t() function transposes, and I removed first 4 columns while transposing since these describe transitions, and I have that info via row names
nrow(SRM.data.PRTC) == ncol(SRM.data.PRTC.t) # confirm same # transitios
ncol(SRM.data.PRTC[,-1:-4]) == nrow(SRM.data.PRTC.t) # confirm same # replicates 

#Replace NA cells with 0; metaMDS() does not handle NA's
SRM.data.PRTC.t.noNA <- SRM.data.PRTC.t
SRM.data.PRTC.t.noNA[is.na(SRM.data.PRTC.t.noNA)] <- 0
head(SRM.data.PRTC.t.noNA)

SRM.data.PRTC.t.noNA

#Make MDS dissimilarity matrix
#
SRM.PRTC.nmds <- metaMDS(SRM.data.PRTC.t.noNA, distance = 'bray', k = 2, trymax = 1000, autotransform = FALSE)
# comm= your data.frame or matrix
# distance= bray, (not sure what this means)
# k= # of dimensions to assess
# trymax = max # iterations to attempt if no solution is reached

# Create Shepard plot, which shows scatter around the regression between the interpoint distances in the final configuration (i.e., the distances between each pair of communities) against their original dissimilarities.
stressplot(SRM.PRTC.nmds) 

#Make figure
plot(SRM.PRTC.nmds)
# site (sample) in black circle
# species (variable) in red ticks






### OLD SCRIPT BELOW HERE ##################################################################################### 

bareCaseInlet <- ave(NormProtAreaAgg$CI.B.1, NormProtAreaAgg$CI.B.2) # create vector for each site/treatment that averages the replicates (of which there are 2)
bareFidalgoBay <- ave(NormProtAreaAgg$FB.B.1, NormProtAreaAgg$FB.B.2)
bareWillapaBay <- ave(NormProtAreaAgg$WB.B.1, NormProtAreaAgg$WB.B.2)
bareSkokomishRiver <- ave(NormProtAreaAgg$SK.B.1, NormProtAreaAgg$SK.B.2)
barePortGamble <- ave(NormProtAreaAgg$PG.B.1, NormProtAreaAgg$PG.B.2)
eelgrassCaseInlet <- ave(NormProtAreaAgg$CI.E.1, NormProtAreaAgg$CI.E.2)
eelgrassFidalgoBay <- ave(NormProtAreaAgg$FB.E.1, NormProtAreaAgg$FB.E.2)
eelgrassSkokomishRiver <- ave(NormProtAreaAgg$SK.E.1, NormProtAreaAgg$SK.E.2)
eelgrassPortGamble <- ave(NormProtAreaAgg$PG.E.1, NormProtAreaAgg$PG.E.2)

NormProtAreaAggAveraged <- data.frame(NormProtAreaAgg$Protein.Name, bareCaseInlet, bareFidalgoBay, barePortGamble, bareSkokomishRiver, bareWillapaBay, eelgrassCaseInlet, eelgrassFidalgoBay, eelgrassPortGamble, eelgrassSkokomishRiver) # combine site/treatment vectors into new dataframe
head(NormProtAreaAggAveraged)

### EDITING GEOID NAME TO MATCH ANNOTATED FILE ### 
NormProtAreaAggAveragedGeoID <- NormProtAreaAggAveraged # copy database, create new one to be used to join to annotated file
GeoID.a <- (gsub("cds.", "", NormProtAreaAggAveragedGeoID$NormProtAreaAgg.Protein.Name, fixed=TRUE)) # remove cds. from protein ID
GeoID.b <- (gsub("\\|m.*", "", GeoID.a)) # remove |m.#### from protein ID
noquote(GeoID.b) # remove quotes from resulting protein ID
NormProtAreaAggAveragedGeoID[1] <- GeoID.b # replace the newly created vector of protein ID's in the full database
head(NormProtAreaAggAveragedGeoID) # confirm changes
nrow(NormProtAreaAggAveragedGeoID) # confirm # rows still same
write.table(NormProtAreaAggAveragedGeoID, "2017-06-30_All-Proteins.tab", quote=F, row.names = F, sep="\t") # Save to file

### UPLOAD ANNOTATED GEODUCK PROTEOME FROM URL ###
GeoduckAnnotations <- read.table(
  "https://raw.githubusercontent.com/sr320/paper-pano-go/master/data-results/Geo-v3-join-uniprot-all0916-condensed.tab",
  sep="\t", header=TRUE, fill=TRUE, stringsAsFactors = FALSE, quote="") #fill=TRUE for empty spaces 
nrow(GeoduckAnnotations) # check that the # rows matches the source (I know this from uploading to Galaxy; not sure how else to do it)
ncol(GeoduckAnnotations) # check that the # columns matches the source
head(GeoduckAnnotations) # inspect; although it's easier to "view" in RStudio
write.table(GeoduckAnnotations, "GeoduckAnnotatinos.tab", quote=F, row.names = F, sep="\t") # Save to file


# MERGE ANNOTATIONS WITH ALL MY PROTEINS ### 
AnnotatedProteins <- merge(x = NormProtAreaAggAveragedGeoID, y = GeoduckAnnotations, by.x="NormProtAreaAgg.Protein.Name", by.y="GeoID")
head(AnnotatedProteins) #confirm merge
nrow(AnnotatedProteins) #count number of proteins that were matched with the list of annotations
write.table(AnnotatedProteins, "2017-06-30_All-Annotated-Proteins.tab", quote=F, row.names = F, sep="\t") #write out .tab file

#### FULL HEATMAP ####

#Install package
install.packages("pheatmap")
library(pheatmap)

#Data should be log(x) or log(x+1) transformed for this analysis. I invert my data (since I had normalized by TIC so values were <1), then log transform.

#Invert data then log transformation for heat map
area.protID.t.inv <- (1/area.protID.t)
is.na(area.protID.t.inv) <- sapply(area.protID.t.inv, is.infinite) # change "Inf" to "NA"
area.protID.t.inv.log <- data.trans(area.protID.t.inv, method = 'log', plot = FALSE)
ncol(area.protID.t.inv.log)
nrow(area.protID.t.inv.log)

#Create heatmap
pheatmap(area.protID.t.inv.log, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'median', show_rownames = T, show_colnames = F)

#Export preliminary heatmap as a .png
png(filename = "2017-07-04_Heatmap-by-median.png")
pheatmap(area.protID.t.inv.log, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'median', show_rownames = T, show_colnames = F)
dev.off()

### Find COEFFICIENT OF VARIANCES & MEANS ###

ProtVariance <- c(var(area2.tra2[-1], area2.tra2[1], na.rm=TRUE))
ProtMeans <- c(colMeans(area2.tra2[-1], na.rm=TRUE))
ProtVar.Means <- data.frame(,y)

ProtVar.Means <- data.frame(NormProtAreaAggAveraged[1], ProtMeans, ProtVariance)
nrow(NormProtAreaAggAveraged[1])
length(ProtMeans)

length(ProtVariance)

head(ProtMeans)
ncol(ProtMeans)
plot(x=ProtMeans, y=ProtVariance, main= "Protein Variance by Mean", xlab= "Mean", ylab= "CV", type="p")
ncol(ProtMeans)
nrow(ProtMeans)
ncol(ProtVariance)
nrow(ProtVariance)


### Experimenting 

SRM.data.test <- SRM.data[1:27,]
cat(SRM.data.test[1,3], SRM.data.test[1,4]) 
# should write a loop to concatenate contents of cells in column 3 with contents in column 4. Fill 
nrow(SRM.data.test)
SRM.data.test[,121] <- 
  
  nrow(SRM.data.test)

length(unique(SRM.data[,3]))
SRM.data <- SRM.data[-1,]
head(SRM.data) # confirm column names changed correctly
ncol(SRM.data)

ncol(SRM.data.test)
SRM.data.test.numeric <- SRM.data.test

### BONEYARD 

### SUM TRANSITION AREA BY PROTEIN ###
# The "Area" values (which have already been normalized by TIC) are peak area for each transition. Sum them to determine total area for each protein

NormProtAreaAgg <- aggregate(cbind(FB.E.1, CI.E.1, PG.B.1, SK.E.1, FB.B.1, WB.B.1, SK.B.1, CI.B.1, PG.E.1, WB.B.2, PG.E.2, FB.E.2, FB.B.2, CI.B.2, SK.E.2, PG.B.2, SK.B.2, CI.E.2) ~ Protein.Name, FUN = sum, data = NormProtArea, na.rm = TRUE, na.action = NULL)
head(NormProtAreaAgg) #confirming that proteins are now summed

### Useful functions from R class notes

# str(SRM.data) # examine a data frame specs
# names() # extract names from data frame
# head(data.frame, n=x) where x=# of rows
# ls() # see what is stored in working directory 


set.seed(2)
community_matrix=matrix(
  sample(1:100,300,replace=T),nrow=10,
  dimnames=list(paste("community",1:10,sep=""),paste("sp",1:30,sep="")))

example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                     k=2) # The number of reduced dimensions
example_NMDS=metaMDS(community_matrix,k=2,trymax=100)
stressplot(example_NMDS)
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)
