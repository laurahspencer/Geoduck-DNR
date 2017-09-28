######## Run ANOSIM for each protein separately. Requires reshape2() library

# Prepare data for anosim by protein, isolating the melted & annotated area data
Arachidonate <- data.melted.plus[grepl(c("Arachidonate"), data.melted.plus$Protein.Name),]
Catalase <- data.melted.plus[grepl(c("Catalase"), data.melted.plus$Protein.Name),]
Cytochrome <- data.melted.plus[grepl(c("Cytochrome"), data.melted.plus$Protein.Name),]
Glycogen <- data.melted.plus[grepl(c("Glycogen"), data.melted.plus$Protein.Name),]
HSP70 <- data.melted.plus[grepl(c("Heat shock 70"), data.melted.plus$Protein.Name),]
HSP90 <- data.melted.plus[grepl(c("HSP 90"), data.melted.plus$Protein.Name),]
Peroxiredoxin <- data.melted.plus[grepl(c("Peroxiredoxin"), data.melted.plus$Protein.Name),]
PDI <- data.melted.plus[grepl(c("Protein disulfide-isomerase"), data.melted.plus$Protein.Name),]
Puromycin <- data.melted.plus[grepl(c("Puromycin"), data.melted.plus$Protein.Name),]
Rab.11B <- data.melted.plus[grepl(c("Rab-11B"), data.melted.plus$Protein.Name),]
NAK <- data.melted.plus[grepl(c("Sodium/potassium"), data.melted.plus$Protein.Name),]
Superoxide <- data.melted.plus[grepl(c("Superoxide"), data.melted.plus$Protein.Name),]
Trifunctional <- data.melted.plus[grepl(c("Trifunctional"), data.melted.plus$Protein.Name),]

# Convert data to long form by casting it, making first column row names
Arachidonate.4anosim <- data.frame(merge(dcast(Arachidonate[,c(1,2,3,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE", all.y = TRUE), row.names = 1)
Catalase.4anosim <- data.frame(merge(dcast(Catalase[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
Cytochrome.4anosim <- data.frame(merge(dcast(Cytochrome[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
Glycogen.4anosim <- data.frame(merge(dcast(Glycogen[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
HSP70.4anosim <- data.frame(merge(dcast(HSP70[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
HSP90.4anosim <- data.frame(merge(dcast(HSP90[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
Peroxiredoxin.4anosim <- data.frame(merge(dcast(Peroxiredoxin[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
PDI.4anosim <- data.frame(merge(dcast(PDI[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
Puromycin.4anosim <- data.frame(merge(dcast(Puromycin[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
Rab.11B.4anosim <- data.frame(merge(dcast(Rab.11B[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
NAK.4anosim <- data.frame(merge(dcast(NAK[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
Superoxide.4anosim <- data.frame(merge(dcast(Superoxide[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)
Trifunctional.4anosim <- data.frame(merge(dcast(Trifunctional[,c(1,2,6)], SAMPLE~Pep.Trans), samples4anosim, by.x="SAMPLE", by.y="SAMPLE"), row.names = 1)

# ANOSIM between sites, each protein

# P=0.42129 
Arachidonate.vegdist <- vegdist(Arachidonate.4anosim[,1:9], 'bray', na.rm=TRUE)
Arachidonate.ANOSIM <-anosim(Arachidonate.vegdist, Arachidonate.4anosim$SITE, permutations = 2000)
summary(Arachidonate.ANOSIM)
# png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ANOSIM%03d.png")
plot(Arachidonate.ANOSIM)
dev.off()

# P=0.74863 
Catalase.vegdist <- vegdist(Catalase.4anosim[,1:9], 'bray', na.rm=TRUE)
Catalase.ANOSIM <-anosim(Catalase.vegdist, grouping=Catalase.4anosim$SITE, permutations = 2000)
summary(Catalase.ANOSIM)
# png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ANOSIM%03d.png")
plot(Catalase.ANOSIM)
# dev.off()

# P=0.81259
Cytochrome.vegdist <- vegdist(Cytochrome.4anosim[,1:9], 'bray', na.rm=TRUE)
Cytochrome.ANOSIM <-anosim(Cytochrome.vegdist, grouping=Cytochrome.4anosim$SITE, permutations = 2000)
summary(Cytochrome.ANOSIM)
# png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ANOSIM%03d.png")
plot(Cytochrome.ANOSIM)
# dev.off()

# P=0.70465 
Glycogen.vegdist <- vegdist(Glycogen.4anosim[,1:9], 'bray', na.rm=TRUE)
Glycogen.ANOSIM <-anosim(Glycogen.vegdist, grouping=Glycogen.4anosim$SITE, permutations = 2000)
summary(Glycogen.ANOSIM)
# png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ANOSIM%03d.png")
plot(Glycogen.ANOSIM)
# dev.off()

# P=0.00049975 
HSP70.vegdist <- vegdist(HSP70.4anosim[,1:6], 'bray', na.rm=TRUE)
HSP70.ANOSIM <-anosim(HSP70.vegdist, grouping=HSP70.4anosim$SITE, permutations = 2000)
summary(HSP70.ANOSIM)
png("Analyses/2017-September_SRM-results/2017-09-16_ANOSIM-HSP70%03d.png")
plot(HSP70.ANOSIM)
dev.off()

# P=0.023488 
HSP90.vegdist <- vegdist(HSP90.4anosim[,1:9], 'bray', na.rm=TRUE)
HSP90.ANOSIM <-anosim(HSP90.vegdist, grouping=HSP90.4anosim$SITE, permutations = 2000)
summary(HSP90.ANOSIM)
png("Analyses/2017-September_SRM-results/2017-09-16_ANOSIM-HSP90%03d.png")
plot(HSP90.ANOSIM)
dev.off()

#P=0.4053 
Peroxiredoxin.vegdist <- vegdist(Peroxiredoxin.4anosim[,1:3], 'bray', na.rm=TRUE)
Peroxiredoxin.ANOSIM <-anosim(Peroxiredoxin.vegdist, grouping=Peroxiredoxin.4anosim$SITE, permutations = 2000)
summary(Peroxiredoxin.ANOSIM)
# png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ANOSIM%03d.png")
plot(Peroxiredoxin.ANOSIM)
# dev.off()

#P=0.028986 
PDI.vegdist <- vegdist(PDI.4anosim[,1:6], 'bray', na.rm=TRUE)
PDI.ANOSIM <-anosim(PDI.vegdist, grouping=PDI.4anosim$SITE, permutations = 2000)
summary(PDI.ANOSIM)
png("Analyses/2017-September_SRM-results/2017-09-16_ANOSIM-PDI0%03d.png")
plot(PDI.ANOSIM)
dev.off()

#P=0.56122 
Puromycin.vegdist <- vegdist(Puromycin.4anosim[,1:9], 'bray', na.rm=TRUE)
Puromycin.ANOSIM <-anosim(Puromycin.vegdist, grouping=Puromycin.4anosim$SITE, permutations = 2000)
summary(Puromycin.ANOSIM)
# png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ANOSIM%03d.png")
plot(PDI.ANOSIM)
# dev.off()

#P= 0.81959 
Rab.11B.vegdist <- vegdist(Rab.11B.4anosim[,1:3], 'bray', na.rm=TRUE)
Rab.11B.ANOSIM <-anosim(Rab.11B.vegdist, grouping=Rab.11B.4anosim$SITE, permutations = 2000)
summary(Rab.11B.ANOSIM)
# png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ANOSIM%03d.png")
plot(Rab.11B.ANOSIM)
# dev.off()

#P=0.64118 
NAK.vegdist <- vegdist(NAK.4anosim[,1:6], 'bray', na.rm=TRUE)
NAK.ANOSIM <-anosim(NAK.vegdist, grouping=NAK.4anosim$SITE, permutations = 2000)
summary(NAK.ANOSIM)
# png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ANOSIM%03d.png")
plot(NAK.ANOSIM)
# dev.off()

#P=0.88556 
Superoxide.vegdist <- vegdist(Superoxide.4anosim[,1:5], 'bray', na.rm=TRUE)
Superoxide.ANOSIM <-anosim(Superoxide.vegdist, grouping=Superoxide.4anosim$SITE, permutations = 2000)
summary(Superoxide.ANOSIM)
# png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ANOSIM%03d.png")
plot(Superoxide.ANOSIM)
# dev.off()

#P=0.089455 
View(Trifunctional.4anosim)
Trifunctional.vegdist <- vegdist(Trifunctional.4anosim[,1:9], 'bray', na.rm=TRUE)
Trifunctional.ANOSIM <-anosim(Trifunctional.vegdist, grouping=Trifunctional.4anosim$SITE, permutations = 2000)
summary(Trifunctional.ANOSIM)
png("Analyses/2017-September_SRM-results/2017-09-16_ANOSIM-TRIFUNCTIONAL0%03d.png")
plot(Trifunctional.ANOSIM)
dev.off()

########## Violin plots of proteins found to be significantly different between sites (ANOSIM)
library(ggplot2)
levels(data.melted.plus$SITE)
data.melted.plus$SITE <- factor(data.melted.plus$SITE, levels=c("WB", "CI", "PG", "FB"))

# HSP 90
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-boxplot-HSP90.png", width = 400, height = 500)
ggplot(subset(data.melted.plus, Pep.Trans %in% "GVVDSEDLPLNISR y7"), aes(x=SITE, y=Area, fill=SITE)) + 
  geom_boxplot(color="black", position = position_dodge()) + 
  ggtitle("Heat shock 90 \nabundance by site") + 
  theme(plot.title = element_text(size=22), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill=FALSE) +
  scale_fill_discrete(labels=c("Willapa Bay", "Case Inlet", "Fidalgo Bay", "Port Gamble")) + 
  ylab("Protein Abundance (Peak Intensity)") +
  coord_flip()
dev.off()


# HSP 70
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-boxplot-HSP70.png", width = 400, height = 500)
ggplot(subset(data.melted.plus, Pep.Trans %in% "IINEPTAAALAYGLDK y12"), aes(x=SITE, y=Area, fill=SITE)) + 
  geom_boxplot(color="black", position = position_dodge()) + 
  ggtitle("Heat shock 70 \nabundance by site") + 
  theme(plot.title = element_text(size=22), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill=FALSE) +
  scale_fill_discrete(labels=c("Willapa Bay", "Case Inlet", "Fidalgo Bay", "Port Gamble")) + 
  ylab("Protein Abundance (Peak Intensity)") +
  coord_flip()
dev.off()


# PDI
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-boxplot-PDI.png", width = 400, height = 500)
ggplot(subset(data.melted.plus, Pep.Trans %in% "DNVVVIGFFK y5"), aes(x=SITE, y=Area, fill=SITE)) + 
  geom_boxplot(color="black", position = position_dodge()) + 
  ggtitle("Protein Disulfide Isomerase \nabundance by site") + 
  theme(plot.title = element_text(size=22), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill=FALSE) +
  scale_fill_discrete(labels=c("Willapa Bay", "Case Inlet", "Fidalgo Bay", "Port Gamble")) + 
  ylab("Protein Abundance (Peak Intensity)") +
  coord_flip()
dev.off()


# Trifunctional Enzyme Subunit B
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-violinplot-TESB.png", width = 400, height = 500)
ggplot(subset(data.melted.plus, Protein.Name %in% SRM.stats[360,2]), aes(x=SITE, y=Area, fill=SITE)) + 
  geom_violin(color="black", position = position_dodge()) + 
  ggtitle("Trifunctional enzyme subunit beta \nabundance by site") + 
  theme(plot.title = element_text(size=22), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill=FALSE) +
  scale_fill_discrete(labels=c("Willapa Bay", "Case Inlet", "Fidalgo Bay", "Port Gamble")) + 
  ylab("Protein Abundance (Peak Intensity)") +
  coord_flip()
dev.off()
