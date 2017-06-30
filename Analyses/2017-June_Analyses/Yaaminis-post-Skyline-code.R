### Import Dataset ### 
### Note: this dataset has already been massaged in Excel: sum area by protein, normalized by TIC, all n/a's removed (replaced with zeros), and columns renamed
setwd("~/Documents/Roberts Lab/Geoduck-DNR/Analyses/2017-June_Analyses")
NormProtArea <- read.csv("2017-06-02_SKYLINE-Total-Protein-Area-NORM.csv", header=TRUE, na.strings = "0") # import local file
head(NormProtArea)

### SUM TRANSITION AREA BY PROTEIN ### 
# The "Area" values (which have already been normalized by TIC) are peak area for each transition. Sum them to determine total area for each protein

NormProtAreaAgg <- aggregate(cbind(FB.E.1, CI.E.1, PG.B.1, SK.E.1, FB.B.1, WB.B.1, SK.B.1, CI.B.1, PG.E.1, WB.B.2, PG.E.2, FB.E.2, FB.B.2, CI.B.2, SK.E.2, PG.B.2, SK.B.2, CI.E.2) ~ Protein.Name, FUN = sum, data = NormProtArea, na.rm = TRUE, na.action = NULL)
head(NormProtAreaAgg) #confirming that proteins are now summed

#### AVERAGE REPLICATES #### 
# FYI No eelgrass sample for Willapa Bay

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

#### CREATE NMDS PLOT ####

#Load the source file for the biostats package

source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
library(vegan)

#Make sure first column of protein names is recognized as row names instead of values - this takes the first column containin protein names, and assigns it to the row "header" names 
area.protID2 <- NormProtAreaAggAveraged[-1]
rownames(area.protID2) <- NormProtAreaAggAveraged[,1]
head(area.protID2)

#Transpose the file so that rows and columns are switched and normalized by log(x+1)
nrow(NormProtAreaAggAveraged)
ncol(NormProtAreaAggAveraged)

area2.t <- t(area.protID2[,1:9])
area2.tra1 <- (area2.t+1)
area2.tra2 <- data.trans(area2.tra1, method = 'log', plot = FALSE)
ncol(area2.tra2)
nrow(area2.tra2)

#Make MDS dissimilarity matrix
proc.nmds <- metaMDS(area2.t, distance = 'bray', k = 2, trymax = 100, autotransform = FALSE)

#Make figure
fig.nmds <- ordiplot(proc.nmds, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=1)
#bare=circle
#eelgrass=triangle
#case=red
#fidalgo=blue
#willapa=black
#skokomish=green
#gamble=magenta

points(fig.nmds, 'sites', col=c('red', 'blue', 'black', 'green', 'magenta','red', 'blue', 'black', 'green', 'magenta'), pch=c(rep(16,5), rep(17,5)))
legend(0.185,0.1, pch=c(rep(16,5), 1, 2), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'))

#### FULL HEATMAP ####

#Install package
install.packages("pheatmap")
library(pheatmap)

#Data should be log(x) or log(x+1) transformed for this analysis, so I'll use my area2.tra dataset.

#Create heatmap
pheatmap(area2.tra2, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'median', show_rownames = T, show_colnames = F)

#Export preliminary heatmap as a .png
png(filename = "Heatmap-by-median.png")
pheatmap(area2.tra2, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'median', show_rownames = T, show_colnames = F)
dev.off()

### Find COEFFICIENT OF VARIANCES & MEDIANS ###

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

#### MERGE WITH GO TERMS ####

#Upload table with proteins and accession codes.
proteinAccessionCodes <- read.table(file = "background-proteome-accession.txt", header = FALSE, col.names = c("averageAreaAdjusted.proteins", "accession", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"))
proteinAccessionCodes <- within(proteinAccessionCodes, rm("V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")) #Removing the extra columns
names(proteinAccessionCodes) <- c("averageProteinAreas.protein", "accession")
head(proteinAccessionCodes) #View uploaded data and confirm changes
head(averageProteinAreasMerged)

#Merge Skyline output with Uniprot information
skylineProteinAccession <- merge(x = averageProteinAreasMerged, y = proteinAccessionCodes, by = "averageProteinAreas.protein")
head(skylineProteinAccession) #confirm merge

#### WRITE OUT MERGED TABLE ####

#Write out alltreatments_DEG_Uniprot as a tab file. Remove row and column names using "row.names" and "col.names" arguments
write.table(skylineProteinAccession, "2017-06-13-Skyline-ProteinAccession-nohead.txt", col.names = F, row.names = F)
