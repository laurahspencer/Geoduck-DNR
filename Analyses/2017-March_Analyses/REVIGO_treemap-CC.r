

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005576","extracellular region",4.572,0.907,0.000,"extracellular region"),
c("GO:0045177","apical part of cell",0.035,0.868,0.000,"apical part of cell"),
c("GO:0016600","flotillin complex",0.001,0.710,0.016,"flotillin complex"),
c("GO:0030137","COPI-coated vesicle",0.015,0.404,0.669,"flotillin complex"),
c("GO:0010008","endosome membrane",0.021,0.399,0.686,"flotillin complex"),
c("GO:0005784","Sec61 translocon complex",0.001,0.436,0.482,"flotillin complex"),
c("GO:0005794","Golgi apparatus",0.265,0.418,0.574,"flotillin complex"),
c("GO:0005829","cytosol",0.807,0.644,0.025,"cytosol"),
c("GO:0030018","Z disc",0.011,0.551,0.284,"cytosol"),
c("GO:0005930","axoneme",0.014,0.562,0.252,"cytosol"),
c("GO:0000805","X chromosome",0.001,0.643,0.056,"X chromosome"),
c("GO:0033698","Rpd3L complex",0.000,0.590,0.354,"X chromosome"),
c("GO:0015629","actin cytoskeleton",0.151,0.511,0.309,"X chromosome"),
c("GO:0030479","actin cortical patch",0.003,0.487,0.581,"X chromosome"),
c("GO:0005737","cytoplasm",38.159,0.775,0.206,"X chromosome"),
c("GO:0005634","nucleus",2.809,0.519,0.202,"X chromosome"),
c("GO:0008290","F-actin capping protein complex",0.007,0.481,0.651,"X chromosome"),
c("GO:0016021","integral component of membrane",35.230,0.891,0.095,"integral component of membrane"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "uniqueness",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
