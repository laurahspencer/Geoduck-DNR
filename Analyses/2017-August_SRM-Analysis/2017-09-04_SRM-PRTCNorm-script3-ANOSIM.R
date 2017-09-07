# Script #3 of PRTC-Normalized SRM data analysis

################################################################################################################
# RUN ANALYSIS OF SIMILARITY (ANOSIM) ON DIFFERENT ITERATIONS OF DATA

###### PREPARE DATA FOR ANOSIM ##########

CI.b <- data.frame(SAMPLE=CI.B.samples, SITE=rep("CI", times=length(CI.B.samples)), TREATMENT=rep("Bare", times=length(CI.B.samples)), BOTH=rep("CI-Bare", times=length(CI.B.samples)))
CI.e <- data.frame(SAMPLE=CI.E.samples, SITE=rep("CI", times=length(CI.E.samples)), TREATMENT=rep("Eelgrass", times=length(CI.E.samples)), BOTH=rep("CI-Eel", times=length(CI.E.samples)))
PG.b <- data.frame(SAMPLE=PG.B.samples, SITE=rep("PG", times=length(PG.B.samples)), TREATMENT=rep("Bare", times=length(PG.B.samples)), BOTH=rep("PG-Bare", times=length(PG.B.samples)))
PG.e <- data.frame(SAMPLE=PG.E.samples, SITE=rep("PG", times=length(PG.E.samples)), TREATMENT=rep("Eelgrass", times=length(PG.E.samples)), BOTH=rep("PG-Eel", times=length(PG.E.samples)))
WB.b <- data.frame(SAMPLE=WB.B.samples, SITE=rep("WB", times=length(WB.B.samples)), TREATMENT=rep("Bare", times=length(WB.B.samples)), BOTH=rep("WB-Bare", times=length(WB.B.samples)))
WB.e <- data.frame(SAMPLE=WB.E.samples, SITE=rep("WB", times=length(WB.E.samples)), TREATMENT=rep("Eelgrass", times=length(WB.E.samples)), BOTH=rep("WB-Eel", times=length(WB.E.samples)))
FB.b <- data.frame(SAMPLE=FB.B.samples, SITE=rep("FB", times=length(FB.B.samples)), TREATMENT=rep("Bare", times=length(FB.B.samples)), BOTH=rep("FB-Bare", times=length(FB.B.samples)))
FB.e <- data.frame(SAMPLE=FB.E.samples, SITE=rep("FB", times=length(FB.E.samples)), TREATMENT=rep("Eelgrass", times=length(FB.E.samples)), BOTH=rep("FB-Eel", times=length(FB.E.samples)))

samples4anosim <- rbind.data.frame(CI.b, CI.e, PG.b, PG.e, WB.b, WB.e, FB.e, FB.b, stringsAsFactors = TRUE)
samples4anosim$SAMPLE <- as.character(samples4anosim$SAMPLE)

###############
# ANOSIM of data (not log transformed, no zeros in data (instead, ignore NAs))
data4anosim <- cbind.data.frame(SRM.data.mean.t[order(rownames(SRM.data.mean.t)),], samples4anosim[order(samples4anosim$SAMPLE),])
data4anosim$SITE <- as.factor(data4anosim$SITE)
data4anosim$TREATMENT <- as.factor(data4anosim$TREATMENT)
data4anosim$BOTH <- as.factor(data4anosim$BOTH)

# ANOSIM between sites
sdms.vegdist <- vegdist(data4anosim[,-116:-119], 'bray', na.rm=TRUE)
site.anos<-anosim(sdms.vegdist, grouping=data4anosim$SITE, permutations = 2000)
summary(site.anos)
plot(site.anos)

# ANOSIM between treatments
treatment.anos<-anosim(sdms.vegdist, grouping=data4anosim$TREATMENT, permutations = 2000)
summary(treatment.anos)
plot(treatment.anos)

# ANOSIM between both site/treatments
siteANDtreatment.anos<-anosim(sdms.vegdist, grouping=data4anosim$BOTH, permutations = 2000)
summary(siteANDtreatment.anos)
plot(siteANDtreatment.anos)

############
# ANOSIM of data (not log transformed, including zeros where no peak found)
data4anosim.noNA <- cbind.data.frame(SRM.data.mean.t.noNA[order(rownames(SRM.data.mean.t.noNA)),], samples4anosim[order(samples4anosim$SAMPLE),])

# ANOSIM between sites, no NA
sdms.noNA.vegdist <- vegdist(data4anosim.noNA[,-116:-119], 'bray', na.rm=TRUE)
site.noNA.anos<-anosim(sdms.vegdist, grouping=data4anosim.noNA$SITE, permutations = 2000)
summary(site.noNA.anos)
plot(site.noNA.anos)

# ANOSIM between treatments, no NA
treatment.noNA.anos<-anosim(sdms.noNA.vegdist, grouping=data4anosim.noNA$TREATMENT, permutations = 2000)
summary(treatment.noNA.anos)
plot(treatment.noNA.anos)

# ANOSIM between both site/treatments, no NA
siteANDtreatment.noNA.anos <- anosim(sdms.noNA.vegdist, grouping=data4anosim.noNA$BOTH, permutations = 2000)
siteANDtreatment.noNA.anos
summary(siteANDtreatment.noNA.anos)
plot(siteANDtreatment.noNA.anos)

#############
# ANOSIM of data after log+1 transformation  
data4anosim.log <- cbind.data.frame(SRM.data.mean.t.log[order(rownames(SRM.data.mean.t.log)),], samples4anosim[order(samples4anosim$SAMPLE),])

# ANOSIM between sites, log+1 transf.
sdms.log.vegdist <- vegdist(data4anosim.log[,-116:-119], 'bray', na.rm=TRUE)
site.log.anos<-anosim(sdms.log.vegdist, grouping=data4anosim.log$SITE, permutations = 2000)
summary(site.log.anos)
plot(site.log.anos)

# ANOSIM between treatments, log+1 transf.
treatment.log.anos<-anosim(sdms.log.vegdist, grouping=data4anosim.log$TREATMENT, permutations = 2000)
summary(treatment.log.anos)
plot(treatment.log.anos)

# ANOSIM between site/treatments, log+1 transf.
siteANDtreatment.log.anos <- anosim(sdms.log.vegdist, grouping=data4anosim.log$BOTH, permutations = 2000)
summary(siteANDtreatment.log.anos)
plot(siteANDtreatment.log.anos)
