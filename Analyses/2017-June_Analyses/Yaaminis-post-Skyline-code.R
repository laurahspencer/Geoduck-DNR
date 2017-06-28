### Import Dataset ### 
### Note: this dataset has already been massaged in Excel: sum area by protein, normalized by TIC, all n/a's removed (replaced with blanks), and columns renamed
NormProtArea <- read.csv("2017-06-02_SKYLINE-Total-Protein-Area-NORM.csv") # import local file; working directory must be set prior to executing this line
head(NormProtArea) 

#### AVERAGE REPLICATES #### 
# FYI No eelgrass sample for Willapa Bay

bareCaseInlet <- ave(NormProtArea$CI.B.1, NormProtArea$CI.B.2)
head(bareCaseInlet)
bareFidalgoBay <- ave(NormProtArea$FB.B.1, NormProtArea$FB.B.2)
bareWillapaBay <- ave(NormProtArea$WB.B.1, NormProtArea$WB.B.2)
bareSkokomishRiver <- ave(NormProtArea$SK.B.1, NormProtArea$SK.B.2)
barePortGamble <- ave(NormProtArea$PG.B.1, NormProtArea$PG.B.2)

eelgrassCaseInlet <- ave(NormProtArea$CI.E.1, NormProtArea$CI.E.2)
eelgrassFidalgoBay <- ave(NormProtArea$FB.E.1, NormProtArea$FB.E.2)
eelgrassSkokomishRiver <- ave(NormProtArea$SK.E.1, NormProtArea$SK.E.2)
eelgrassPortGamble <- ave(NormProtArea$PG.E.1, NormProtArea$PG.E.2)

NormProtAreaAveraged <- data.frame(NormProtArea$Protein.Name, bareCaseInlet, bareFidalgoBay, barePortGamble, bareSkokomishRiver, bareWillapaBay, eelgrassCaseInlet, eelgrassFidalgoBay, eelgrassPortGamble, eelgrassSkokomishRiver, eelgrassWillapaBay)

#### CREATE NMDS PLOT ####

head(NormProtAreaAveraged) #Confirm changes

#Load the source file for the biostats package

source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
install.packages("vegan") #Install vegan package
library(vegan)

#Make sure first column of protein names is recognized as row names instead of values
area.protID2 <- NormProtAreaAveraged[-1]
rownames(area.protID2) <- NormProtAreaAveraged[,1]

#Transpose the file so that rows and columns are switched and normalized by log(x+1)
area2.t <- t(area.protID2[,1:10])
area2.tra <- (area2.t+1)
area2.tra <- data.trans(area2.tra, method = 'log', plot = FALSE)

#Make MDS dissimilarity matrix
proc.nmds <- metaMDS(area2.tra, distance = 'bray', k = 2, trymax = 100, autotransform = FALSE)

#Make figure
fig.nmds <- ordiplot(proc.nmds, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=0.5)
#bare=circle
#eelgrass=triangle
#case=red
#fidalgo=blue
#willapa=black
#skokomish=green
#gamble=magenta

points(fig.nmds, 'sites', col=c('red', 'blue', 'black', 'green', 'magenta','red', 'blue', 'black', 'green', 'magenta'), pch=c(rep(16,5), rep(17,5)))
legend(-0.045,0.025, pch=c(rep(16,5), 1, 2), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'))

#### FULL HEATMAP ####

#Install package
install.packages("pheatmap")
library(pheatmap)

#Data should be log(x) or log(x+1) transformed for this analysis, so I'll use my area2.tra dataset.

#Create heatmap
pheatmap(area2.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)

#Export preliminary heatmap as a .png
png(filename = "fullHeatmap.png")
pheatmap(area2.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)
dev.off()

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
