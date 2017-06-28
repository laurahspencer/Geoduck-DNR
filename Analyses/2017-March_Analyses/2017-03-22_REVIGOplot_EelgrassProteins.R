

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","Proteins_Represented","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0006096","glycolytic process", 0.545, 4.780,-2.169, 4.844,-6.4498,0.591,0.000),
c("GO:0006891","intra-Golgi vesicle-mediated transport", 0.064, 0.070, 5.934, 3.916,-2.8621,0.882,0.999),
c("GO:0051595","response to methylglyoxal", 0.001, 0.264, 5.325, 1.964,-1.5440,0.948,0.000),
c("GO:0098609","cell-cell adhesion", 0.251,-3.626, 5.167, 4.507,-1.4134,0.948,0.000),
c("GO:0016575","histone deacetylation", 0.048,-2.423,-6.151, 3.787,-5.2652,0.651,0.000),
c("GO:0045769","negative regulation of asymmetric cell division", 0.000,-5.530, 2.343, 1.447,-1.5440,0.808,0.102),
c("GO:0010971","positive regulation of G2/M transition of mitotic cell cycle", 0.009,-1.532, 1.850, 3.046,-1.9829,0.768,0.128),
c("GO:0003341","cilium movement", 0.023, 0.472,-0.229, 3.467,-5.2652,0.726,0.000),
c("GO:0010501","RNA secondary structure unwinding", 0.025, 3.553,-8.363, 3.507,-2.1537,0.831,0.000),
c("GO:0051823","regulation of synapse structural plasticity", 0.001,-5.468,-3.321, 2.167,-1.5440,0.724,0.256),
c("GO:0046166","glyceraldehyde-3-phosphate biosynthetic process", 0.001, 6.700, 1.770, 2.283,-3.2132,0.764,0.321),
c("GO:0090611","ubiquitin-independent protein catabolic process via the multivesicular body sorting pathway", 0.001, 0.914,-8.670, 2.130,-1.5440,0.832,0.352),
c("GO:0015966","diadenosine tetraphosphate biosynthetic process", 0.001, 6.839,-0.495, 1.869,-1.5440,0.752,0.352),
c("GO:0006424","glutamyl-tRNA aminoacylation", 0.065, 3.750,-4.959, 3.922,-2.9162,0.617,0.353),
c("GO:0006108","malate metabolic process", 0.088, 5.460,-3.146, 4.051,-2.1587,0.722,0.000),
c("GO:0010817","reg. of hormone levels", 0.161,-6.798, 0.601, 4.314,-1.3710,0.908,0.000),
c("GO:0007010","cytoskeleton organization", 0.786,-5.030,-4.904, 5.004,-1.9791,0.746,0.390),
c("GO:0006002","fructose 6-phosphate metabolic process", 0.060, 6.295, 2.928, 3.885,-2.5266,0.839,0.414),
c("GO:0034063","stress granule assembly", 0.005,-4.852,-5.355, 2.786,-1.8393,0.752,0.000),
c("GO:0000054","ribosomal subunit export from nucleus", 0.037,-4.471,-2.861, 3.673,-1.5440,0.784,0.455),
c("GO:0097056","selenocysteinyl-tRNA(Sec) biosynthetic process", 0.034, 4.948,-6.792, 3.643,-2.1587,0.754,0.470),
c("GO:0030317","flagellated sperm motility", 0.017, 2.354, 2.849, 3.333,-3.7570,0.750,0.000),
c("GO:1901675","negative regulation of histone H3-K27 acetylation", 0.000,-3.185,-5.609, 1.519,-1.5440,0.707, 0.000),
c("GO:0036159","inner dynein arm assembly", 0.005,-3.318,-4.372, 2.803,-3.1226,0.620,0.545),
c("GO:0033962","cytoplasmic mRNA processing body assembly", 0.011,-4.700,-5.567, 3.145,-1.6142,0.745,0.587),
c("GO:0006433","prolyl-tRNA aminoacylation", 0.052, 4.128,-5.153, 3.825,-1.3710,0.621,0.599),
c("GO:0006434","seryl-tRNA aminoacylation", 0.053, 3.816,-5.258, 3.836,-2.3846,0.621,0.600),
c("GO:0006436","tryptophanyl-tRNA aminoacylation", 0.054, 4.114,-4.895, 3.843,-1.3710,0.621,0.601),
c("GO:0007409","axonogenesis", 0.118,-2.811,-3.808, 4.179,-1.6990,0.659,0.609),
c("GO:0021873","forebrain neuroblast division", 0.002,-4.826, 2.147, 2.316,-1.5440,0.796,0.629),
c("GO:0006094","gluconeogenesis", 0.262, 6.752,-2.543, 4.527,-2.1770,0.749,0.648),
c("GO:0007018","microtubule-based movement", 0.287, 1.096,-0.589, 4.567,-1.8687,0.694,0.000),
c("GO:0006099","tricarboxylic acid cycle", 0.469, 4.421,-2.531, 4.780,-4.6799,0.683,0.731),
c("GO:0051653","spindle localization", 0.020, 0.114, 2.379, 3.406,-1.3710,0.726,0.742),
c("GO:0006418","tRNA aminoacylation for protein translation", 1.099, 3.321,-4.292, 5.149,-2.3232,0.558,0.749),
c("GO:0007051","spindle organization", 0.083,-2.103,-2.816, 4.029,-1.6142,0.584,0.771),
c("GO:0006890","retrograde vesicle-mediated transport, Golgi to ER", 0.047, 0.956, 5.933, 3.777,-2.3401,0.883,0.772),
c("GO:0006098","pentose-phosphate shunt", 0.287, 5.226,-1.025, 4.566,-1.9078,0.659,0.800));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$Proteins_Represented <- as.numeric( as.character(one.data$Proteins_Represented) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = Proteins_Represented), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("darkred", "tomato2", "seagreen3", "darkslateblue"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = Proteins_Represented), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.001, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), label.padding = unit(0.5, "lines"), check_overlap = TRUE, colour = I(alpha("black", 0.85)), size = 5.5, fontface = "bold");
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

ggsave("2017-03-22_REVIGOplot_EelgrassProteins.pdf");
