

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
revigo.data <- rbind(c("GO:0000086","G2/M transition of mitotic cell cycle",0.007,0.880,0.000,"G2/M transition of mitotic cell cycle"),
c("GO:0000077","DNA damage checkpoint",0.011,0.808,0.559,"G2/M transition of mitotic cell cycle"),
c("GO:0007155","cell adhesion",0.564,0.930,0.000,"cell adhesion"),
c("GO:0050829","defense response to Gram-negative bacterium",0.007,0.829,0.000,"defense response to Gram-negative bacterium"),
c("GO:0034599","cellular response to oxidative stress",0.191,0.798,0.385,"defense response to Gram-negative bacterium"),
c("GO:0010815","bradykinin catabolic process",0.000,0.797,0.013,"bradykinin catabolism"),
c("GO:0005980","glycogen catabolic process",0.025,0.742,0.570,"bradykinin catabolism"),
c("GO:0006401","RNA catabolic process",0.274,0.757,0.251,"bradykinin catabolism"),
c("GO:0006096","glycolytic process",0.522,0.624,0.645,"bradykinin catabolism"),
c("GO:0001700","embryonic development via the syncytial blastoderm",0.001,0.807,0.026,"embryonic development via the syncytial blastoderm"),
c("GO:0043046","DNA methylation involved in gamete generation",0.001,0.729,0.410,"embryonic development via the syncytial blastoderm"),
c("GO:0001525","angiogenesis",0.026,0.790,0.530,"embryonic development via the syncytial blastoderm"),
c("GO:0005975","carbohydrate metabolic process",5.984,0.885,0.033,"carbohydrate metabolism"),
c("GO:0000448","cleavage in ITS2 between 5.8S rRNA and LSU-rRNA of tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)",0.000,0.797,0.054,"cleavage in ITS2 between 5.8S rRNA and LSU-rRNA of tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)"),
c("GO:0051693","actin filament capping",0.003,0.768,0.268,"cleavage in ITS2 between 5.8S rRNA and LSU-rRNA of tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)"),
c("GO:0030031","cell projection assembly",0.184,0.759,0.508,"cleavage in ITS2 between 5.8S rRNA and LSU-rRNA of tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)"),
c("GO:0030030","cell projection organization",0.378,0.765,0.373,"cleavage in ITS2 between 5.8S rRNA and LSU-rRNA of tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)"),
c("GO:0006915","apoptotic process",0.247,0.831,0.065,"apoptotic process"),
c("GO:0006103","2-oxoglutarate metabolic process",0.010,0.778,0.066,"2-oxoglutarate metabolism"),
c("GO:0006048","UDP-N-acetylglucosamine biosynthetic process",0.026,0.732,0.193,"2-oxoglutarate metabolism"),
c("GO:0045456","ecdysteroid biosynthetic process",0.000,0.681,0.642,"2-oxoglutarate metabolism"),
c("GO:0009396","folic acid-containing compound biosynthetic process",0.294,0.643,0.524,"2-oxoglutarate metabolism"),
c("GO:0006564","L-serine biosynthetic process",0.082,0.699,0.300,"2-oxoglutarate metabolism"),
c("GO:0006695","cholesterol biosynthetic process",0.002,0.758,0.122,"2-oxoglutarate metabolism"),
c("GO:0044209","AMP salvage",0.066,0.708,0.447,"2-oxoglutarate metabolism"),
c("GO:0034641","cellular nitrogen compound metabolic process",33.428,0.836,0.085,"cellular nitrogen compound metabolism"),
c("GO:0051301","cell division",1.357,0.818,0.098,"cell division"),
c("GO:0007154","cell communication",4.358,0.809,0.135,"cell division"));

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
