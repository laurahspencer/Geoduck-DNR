# Script #4 in data processing for PRTC-Normalized SRM data

######## CALCULATE & PLOT MEAN & STANDARD ERROR FOR SAMPLES BY SITE FOR EACH PROTEIN ########
# Use data4anosim.noNA dataset 

View(data4anosim.noNA)

#melt un-meaned data
library(reshape2)
data.melted <- melt(data4anosim.noNA, id=c("SAMPLE", "SITE", "TREATMENT", "BOTH"), variable.name = "Transition", value.name = "Area")
View(data.melted)

# Merge protein names back to abundance data
SRM.proteins <- data.frame(SRM.data.numeric[2:142,1:4]) #protein name to each transition
SRM.proteins[,1] <- sub(" cds.*", "", SRM.proteins[,1])
data.melted.plus <- merge(x=data.melted, y=SRM.proteins, by.x = "Transition", by.y = "row.names", all.x=TRUE, all.y=FALSE)
colnames(data.melted.plus)[1] <- "Pep.Trans"
# write.csv(data.melted.plus, file="Analyses/2017-September_SRM-results/2017-09-04_SRM-data-PRTCNorm-melted-annotated.csv")

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
# write.csv(SRM.stats, file="Analyses/2017-September_SRM-results/2017-09-04_SRM-data-PRTCNorm-summary-stats.csv")

#### PLOT SITE MEANS FOR EACH PROTEIN, BROKEN INTO PROTEIN, PEPTIDE & TRANSITION
library(ggplot2)

# Arachidonate 5-lipoxygenase
ggplot(subset(SRM.stats, Protein.Name %in% "Arachidonate 5-lipoxygenase"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Arachidonate 5-lipoxygenase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Arachidonate 5-lipoxygenase"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Arachidonate 5-lipoxygenase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Arachidonate 5-lipoxygenase"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Arachidonate 5-lipoxygenase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

# Catalase
ggplot(subset(SRM.stats, Protein.Name %in% "Catalase"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Catalase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Catalase"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Catalase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Catalase"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Catalase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

# Cytochrome P450
ggplot(subset(SRM.stats, Protein.Name %in% "Cytochrome P450"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Cytochrome P450") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Cytochrome P450"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Cytochrome P450") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Cytochrome P450"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Cytochrome P450") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

# Glycogen phosphorylase, muscle form

ggplot(subset(SRM.stats, Protein.Name %in% "Glycogen phosphorylase, muscle form"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Glycogen phosphorylase, muscle form") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Glycogen phosphorylase, muscle form"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Glycogen phosphorylase, muscle form") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Glycogen phosphorylase, muscle form"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Glycogen phosphorylase, muscle form") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

# Heat shock 70 kDa protein

ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock 70 kDa protein"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 70") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock 70 kDa protein"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 70") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock 70 kDa protein"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 70") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))


# Heat shock protein HSP 90-alpha
ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock protein HSP 90-alpha"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 90") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock protein HSP 90-alpha"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 90") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Heat shock protein HSP 90-alpha"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Heat shock 90") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))


# Peroxiredoxin-1

ggplot(subset(SRM.stats, Protein.Name %in% "Peroxiredoxin-1"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Peroxiredoxin-1") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Peroxiredoxin-1"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Peroxiredoxin-1") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Peroxiredoxin-1"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Peroxiredoxin-1") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))


# Protein disulfide-isomerase (PDI)

ggplot(subset(SRM.stats, Protein.Name %in% "Protein disulfide-isomerase (PDI)"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Protein disulfide-isomerase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Protein disulfide-isomerase (PDI)"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Protein disulfide-isomerase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Protein disulfide-isomerase (PDI)"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Protein disulfide-isomerase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))


# Puromycin-sensitive aminopeptidase

ggplot(subset(SRM.stats, Protein.Name %in% "Puromycin-sensitive aminopeptidase"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Puromycin-sensitive aminopeptidase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Puromycin-sensitive aminopeptidase"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Puromycin-sensitive aminopeptidase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Puromycin-sensitive aminopeptidase"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Puromycin-sensitive aminopeptidase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))


# Ras-related protein Rab-11B

ggplot(subset(SRM.stats, Protein.Name %in% "Ras-related protein Rab-11B"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Ras-related protein Rab-11B") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Ras-related protein Rab-11B"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Ras-related protein Rab-11B") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Ras-related protein Rab-11B"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Ras-related protein Rab-11B") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))


# Sodium/potassium-transporting ATPase subunit alpha-4 #NOT WORKING

ggplot(subset(SRM.stats, Protein.Name %in% "Sodium/potassium-transporting ATPase subunit alpha-4"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Na/K-transporting ATPase subunit alpha-4") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Sodium/potassium-transporting ATPase subunit alpha-4"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Na/K-transporting ATPase subunit alpha-4") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Sodium/potassium-transporting ATPase subunit alpha-4"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Na/K-transporting ATPase subunit alpha-4") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

# Superoxide dismutase

ggplot(subset(SRM.stats, Protein.Name %in% "Superoxide dismutase"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Superoxide dismutase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Superoxide dismutase"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Superoxide dismutase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Superoxide dismutase"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Superoxide dismutase") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

# Trifunctional enzyme subunit beta, mitochondrial

ggplot(subset(SRM.stats, Protein.Name %in% "Trifunctional enzyme subunit beta, mitochondrial"), aes(x=Protein.Name, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Trifunctional enzyme subunit beta, mitochondrial") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Trifunctional enzyme subunit beta, mitochondrial"), aes(x=Peptide.Sequence, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Trifunctional enzyme subunit beta, mitochondrial") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

ggplot(subset(SRM.stats, Protein.Name %in% "Trifunctional enzyme subunit beta, mitochondrial"), aes(x=Pep.Trans, y=Area, fill=SITE)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + ggtitle("Trifunctional enzyme subunit beta, mitochondrial") +
  geom_errorbar(aes(ymin=Area-st.err, ymax=Area+st.err), width=.2, position=position_dodge(.9))

