# Script #3 in data processing for NOT NORMALIZED DATA

######## CALCULATE & PLOT MEAN & STANDARD ERROR FOR SAMPLES BY SITE FOR EACH PROTEIN ########
# Use data4anosim.noNA dataset 

#melt data to prepare for ggplot
library(reshape2)
View(data4anosim.noNA)
data.melted <- melt(data4anosim.noNA, id=c("SAMPLE", "SITE", "TREATMENT", "BOTH"), variable.name = "Transition", value.name = "Area")
View(data.melted)

####### check out total abundance & other stats by site
TotAbund.SITE <- do.call(data.frame, aggregate(Area ~ SITE, data.melted, function(x) c(sum=sum(x), sd=sd(x), range=range(x), min=min(x), max=max(x))))
View(TotAbund.SITE)

#Bar plot showing total abundance by site, with error bars
library(ggplot2)
TotSum.bar <- ggplot(TotAbund.SITE, aes(x=SITE, y=Area.sum, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Total Abundance by Site") +
  geom_errorbar(aes(ymin=Area.sum-Area.sd, ymax=Area.sum+Area.sd), width=.2, position=position_dodge(.9)) 
ggplot(data.melted, aes(SITE, Area)) + geom_boxplot(aes(colour=SITE)) + ggtitle("Abundance Distribution by Site")
ggplot(data.melted, aes(SITE, Area)) + geom_violin(aes(colour=SITE)) + ggtitle("Abundance Distribution by Site")
ggplot(data.melted, aes(SITE, Area)) + geom_boxplot(aes(colour=SAMPLE)) + ggtitle("Abundance Distribution in each sample, grouped by site")
ggplot(data.melted, aes(SITE, Area)) + geom_violin(aes(colour=SAMPLE)) + ggtitle("Abundance Distribution in each sample, grouped by site")
ggplot(data.melted, aes(SITE, Area)) + geom_boxplot(aes(colour=TREATMENT)) + ggtitle("Abundance Distribution in each treatment, grouped by site")
ggplot(data.melted, aes(SITE, Area)) + geom_violin(aes(colour=TREATMENT)) + ggtitle("Abundance Distribution in each treatment, grouped by site")


#Log+1 transform summary data for plotting purposes
data.melted.log <- data.melted
data.melted.log$Area <- log(data.melted.log$Area+1)
View(data.melted.log)
TotAbund.SITE.log <- do.call(data.frame, aggregate(Area ~ SITE, data.melted.log, function(x) c(sum=sum(x), sd=sd(x), range=range(x), min=min(x), max=max(x))))
TotSum.log.bar <- ggplot(TotAbund.SITE.log, aes(x=SITE, y=Area.sum, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Total Abundance by Site, log+1 transformed") +
  geom_errorbar(aes(ymin=Area.sum-Area.sd, ymax=Area.sum+Area.sd), width=.2, position=position_dodge(.9)) 
area.log.box.site <- ggplot(data.melted.log, aes(SITE, Area)) + geom_boxplot(aes(colour=SITE)) + ggtitle("Abundance Distribution in each site, log+1 transformed")
area.log.viol.site <- ggplot(data.melted.log, aes(SITE, Area)) + geom_violin(aes(colour=SITE)) + ggtitle("Abundance Distribution in each site, log+1 transformed")
area.log.box.smpl <- ggplot(data.melted.log, aes(SITE, Area)) + geom_boxplot(aes(fill=TREATMENT, colour=SAMPLE)) + ggtitle("Abundance Distribution in each sample log+1 transformed, by site") + guides(colour=FALSE)
area.log.viol.smpl <- ggplot(data.melted.log, aes(SITE, Area)) + geom_violin(aes(colour=SAMPLE)) + ggtitle("Abundance Distribution in each sample log+1 transformed, by site")
area.log.box.trmt <- ggplot(data.melted.log, aes(SITE, Area)) + geom_boxplot(aes(colour=TREATMENT)) + ggtitle("Abundance Distribution in each treatment log+1 transformed, by site")
area.log.viol.trmt <- ggplot(data.melted.log, aes(SITE, Area)) + geom_violin(aes(colour=TREATMENT)) + ggtitle("Abundance Distribution in each treatment log+1 transformed, by site")

png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-All-Sample-plot%03d.png")
TotSum.bar
area.log.box.site
area.log.box.smpl
area.log.box.trmt
dev.off()
#######

# Merge protein names back to abundance data
SRM.proteins <- data.frame(SRM.data.numeric[2:142,1:4]) #protein name to each transition
SRM.proteins[,1] <- sub(" cds.*", "", SRM.proteins[,1])
data.melted.plus <- merge(x=data.melted, y=SRM.proteins, by.x = "Transition", by.y = "row.names", all.x=TRUE, all.y=FALSE)
colnames(data.melted.plus)[1] <- "Pep.Trans"
# write.csv(data.melted.plus, file="Analyses/2017-September_SRM-results/2017-09-04_SRM-data-notNORM-melted-annotated.csv")

# Write program with summary statistics http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_summary <- function(data, varname, groupnames){
  require(plyr)
  require(plotrix)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      st.err = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
SRM.stats <-  data_summary(data.melted.plus, varname="Area", groupnames=c("SITE", "Protein.Name", "Peptide.Sequence", "Pep.Trans")) #run program with my data
View(SRM.stats)
write.csv(SRM.stats, file="Analyses/2017-September_SRM-results/2017-09-04_SRM-data-notNORM-summary-stats.csv")

#### PLOT SITE MEANS FOR EACH PROTEIN, BROKEN INTO PROTEIN, PEPTIDE & TRANSITION
library(ggplot2)

# Arachidonate 5-lipoxygenase
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ArachidonatePro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Arachidonate 5-lipoxygenase"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Arachidonate 5-lipoxygenase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ArachidonatePep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Arachidonate 5-lipoxygenase"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Arachidonate 5-lipoxygenase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-ArachidonateTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Arachidonate 5-lipoxygenase"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Arachidonate 5-lipoxygenase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Catalase
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-CatalasePro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Catalase"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Catalase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-CatalasePep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Catalase"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Catalase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-CatalaseTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Catalase"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Catalase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Cytochrome P450
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-CytochromePro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Cytochrome P450"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Cytochrome P450") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-CytochromePep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Cytochrome P450"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Cytochrome P450") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-CytochromeTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Cytochrome P450"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Cytochrome P450") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Glycogen phosphorylase, muscle form
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-GlycogenPro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Glycogen phosphorylase, muscle form"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Glycogen phosphorylase, muscle form") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-GlycogenPep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Glycogen phosphorylase, muscle form"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Glycogen phosphorylase, muscle form") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-GlycogenTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Glycogen phosphorylase, muscle form"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Glycogen phosphorylase, muscle form") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Heat shock 70 kDa protein
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-HSP70Pro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock 70 kDa protein"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 70") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-HSP70Pep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock 70 kDa protein"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 70") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-HSP70Tran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock 70 kDa protein"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 70") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Heat shock protein HSP 90-alpha
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-HSP90Pro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock protein HSP 90-alpha"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 90") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-HSP90Pep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock protein HSP 90-alpha"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 90") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-HSP90Tran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock protein HSP 90-alpha"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 90") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Peroxiredoxin-1

png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-PeroxiredoxinPro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Peroxiredoxin-1"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Peroxiredoxin-1") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-PeroxiredoxinPep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Peroxiredoxin-1"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Peroxiredoxin-1") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-PeroxiredoxinTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Peroxiredoxin-1"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Peroxiredoxin-1") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Protein disulfide-isomerase (PDI)

png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-PDIPro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Protein disulfide-isomerase (PDI)"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Protein disulfide-isomerase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-PDIPep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Protein disulfide-isomerase (PDI)"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Protein disulfide-isomerase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-PDITran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Protein disulfide-isomerase (PDI)"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Protein disulfide-isomerase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Puromycin-sensitive aminopeptidase
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-PuromycinPro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Puromycin-sensitive aminopeptidase"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Puromycin-sensitive aminopeptidase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-PuromycinPep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Puromycin-sensitive aminopeptidase"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Puromycin-sensitive aminopeptidase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-PuromycinTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Puromycin-sensitive aminopeptidase"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Puromycin-sensitive aminopeptidase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Ras-related protein Rab-11B
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-RasrelatedRabPro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Ras-related protein Rab-11B"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Ras-related protein Rab-11B") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-RasrelatedRabPep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Ras-related protein Rab-11B"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Ras-related protein Rab-11B") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-RasrelatedRabTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Ras-related protein Rab-11B"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Ras-related protein Rab-11B") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Sodium/potassium-transporting ATPase subunit alpha-4 
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-NAKtransprotingPro.png")
ggplot(subset(SRM.stats, Protein.Name %in% SRM.stats[435,2]), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Na/K-transporting ATPase subunit alpha-4") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-NAKtransprotingPep.png")
ggplot(subset(SRM.stats, Protein.Name %in% SRM.stats[435,2]), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Na/K-transporting ATPase subunit alpha-4") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-NAKtransprotingTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% SRM.stats[435,2]), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Na/K-transporting ATPase subunit alpha-4") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Superoxide dismutase

png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-SuperoxidePro.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Superoxide dismutase"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Superoxide dismutase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-SuperoxidePep.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Superoxide dismutase"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Superoxide dismutase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-SuperoxideTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% "Superoxide dismutase"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Superoxide dismutase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

# Trifunctional enzyme subunit beta, mitochondrial
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-TrifunctEnzymePro.png")
ggplot(subset(SRM.stats, Protein.Name %in% SRM.stats[460,2]), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Trifunctional enzyme subunit beta, mitochondrial") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-TrifunctEnzymePep.png")
ggplot(subset(SRM.stats, Protein.Name %in% SRM.stats[460,2]), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Trifunctional enzyme subunit beta, mitochondrial") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()
png("Analyses/2017-September_SRM-results/2017-09-04_NotNORM-plot-TrifunctEnzymeTran.png")
ggplot(subset(SRM.stats, Protein.Name %in% SRM.stats[460,2]), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Trifunctional enzyme subunit beta, mitochondrial") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))
dev.off()

