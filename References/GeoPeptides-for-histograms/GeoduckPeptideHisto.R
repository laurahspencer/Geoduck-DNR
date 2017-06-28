#import peptide .csv file into R; this file is a combined dataset with 10 geoduck sample data. Samples were randomly selected as the following: 11, 24, 27, 29, 31, 33, 65, 67, 84, 87
PeptideData <- read.csv("HistogramData.csv")

#plot m/z values from .csv data file as histogram with 10 breaks,
hist(PeptideData$MZratio, breaks=10, col=88, main="Geoduck Peptide m/z", xlab="m/z value", ylab="Frequency", plot=TRUE)
