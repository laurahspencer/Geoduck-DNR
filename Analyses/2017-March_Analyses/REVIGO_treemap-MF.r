

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
revigo.data <- rbind(c("GO:0003872","6-phosphofructokinase activity",0.048,0.919,0.000,"6-phosphofructokinase activity"),
c("GO:0008289","lipid binding",0.237,0.883,0.000,"lipid binding"),
c("GO:0030234","enzyme regulator activity",0.437,0.930,0.000,"enzyme regulator activity"),
c("GO:0004051","arachidonate 5-lipoxygenase activity",0.000,0.877,0.010,"arachidonate 5-lipoxygenase activity"),
c("GO:0003863","3-methyl-2-oxobutanoate dehydrogenase (2-methylpropanoyl-transferring) activity",0.004,0.811,0.430,"arachidonate 5-lipoxygenase activity"),
c("GO:0003995","acyl-CoA dehydrogenase activity",0.368,0.854,0.168,"arachidonate 5-lipoxygenase activity"),
c("GO:0004029","aldehyde dehydrogenase (NAD) activity",0.042,0.797,0.512,"arachidonate 5-lipoxygenase activity"),
c("GO:0004043","L-aminoadipate-semialdehyde dehydrogenase activity",0.001,0.814,0.190,"arachidonate 5-lipoxygenase activity"),
c("GO:0018478","malonate-semialdehyde dehydrogenase (acetylating) activity",0.002,0.812,0.425,"arachidonate 5-lipoxygenase activity"),
c("GO:0003990","acetylcholinesterase activity",0.001,0.884,0.012,"acetylcholinesterase activity"),
c("GO:0004190","aspartic-type endopeptidase activity",0.513,0.856,0.624,"acetylcholinesterase activity"),
c("GO:0004177","aminopeptidase activity",0.371,0.857,0.140,"acetylcholinesterase activity"),
c("GO:0008409","5'-3' exonuclease activity",0.056,0.876,0.302,"acetylcholinesterase activity"),
c("GO:0016491","oxidoreductase activity",14.657,0.923,0.024,"oxidoreductase activity"),
c("GO:0044548","S100 protein binding",0.001,0.836,0.029,"S100 protein binding"),
c("GO:0048365","Rac GTPase binding",0.002,0.833,0.316,"S100 protein binding"),
c("GO:0042169","SH2 domain binding",0.002,0.833,0.328,"S100 protein binding"),
c("GO:0003779","actin binding",0.080,0.815,0.400,"S100 protein binding"),
c("GO:0005509","calcium ion binding",0.365,0.848,0.045,"calcium ion binding"),
c("GO:0005507","copper ion binding",0.371,0.847,0.281,"calcium ion binding"),
c("GO:0051537","2 iron, 2 sulfur cluster binding",0.441,0.836,0.047,"2 iron, 2 sulfur cluster binding"),
c("GO:0003723","RNA binding",5.860,0.823,0.063,"RNA binding"),
c("GO:0003677","DNA binding",13.924,0.818,0.497,"RNA binding"),
c("GO:0019003","GDP binding",0.003,0.818,0.091,"GDP binding"),
c("GO:0005524","ATP binding",13.885,0.757,0.369,"GDP binding"),
c("GO:0050661","NADP binding",0.744,0.791,0.165,"GDP binding"),
c("GO:0005525","GTP binding",1.692,0.762,0.546,"GDP binding"),
c("GO:0051287","NAD binding",1.108,0.787,0.627,"GDP binding"));

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
