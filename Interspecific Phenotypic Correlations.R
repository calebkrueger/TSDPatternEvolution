# Interspecific Phenotypic Correlations

# Load the necessary R packages for all analyses
library("embryogrowth")
library("phytools")
library("RRPP")
library("ggplot2")
library("coda")

############################################################
# Data filtering

# Read in sex-ratio data from the ROSIE database
data<-read.csv("https://raw.githubusercontent.com/calebkrueger/ROSIE/main/Database.csv")

# Filter to include data from constant incubation experiments without chemical manipulations
# Exclude species with GSD, multiple instances of the same data, and influential/large data points
data<-subset(data, SDM=="TSD" & Incubation_Setup=="Constant" & Data_Elsewhere==0 & Total_Sexed<1000 & is.na(data$Chemical_Treatment))

# Exclude data with large reported fluctuations in incubation temperatures
data<-subset(data, is.na(data$Fluctuation) | Fluctuation<=1)

# Filter out rows without temperature data or numbers of males and/or females
data<-data[complete.cases(data$Mean_Temperature),]
data<-data[complete.cases(data$Males),]
data<-data[complete.cases(data$Females),]

# Exclude species with one data point or no mixed-sex data
data<-data[!data$Species %in% c("Cuora_trifasciata",
                                "Cuora_zhoui",
                                "Geochelone_elegans",
                                "Graptemys_barbouri",
                                "Graptemys_nigrinoda",
                                "Graptemys_pulchra",
                                "Manouria_emys",
                                "Manouria_impressa",
                                "Mauremys_nigricans",
                                "Peltocephalus_dumerilianus",
                                "Terrapene_ornata",
                                "Testudo_kleinmanni"),]

# Exclude multi-species data
data<-data[!data$Species %in% "Graptemys_pseudogeographica/ouachitensis",]

# Read in file with sex-determining mechanism classifications from ROSIE
SDM<-read.csv("https://raw.githubusercontent.com/calebkrueger/ROSIE/main/SDM.csv")

# Subset out data on species with pattern II TSD
forII<-subset(SDM, Type=="II" | Type=="Ia, II")
II<-data[data$Species %in% forII$Species,]

# Subset out data on species with pattern Ia TSD
forIa<-subset(SDM, Type=="Ia")
Ia<-data[data$Species %in% forIa$Species,]

# Separate data by species into lists for pattern II and pattern Ia species
species.II<-list()
for(i in unique(II$Species)) {
  species.II[[i]]<-II[II$Species==i,]
}
species.Ia<-list()
for(i in unique(Ia$Species)) {
  species.Ia[[i]]<-Ia[Ia$Species==i,]
}

# Remove some irrelevant data from lower temps so models can fit
species.II$Pelomedusa_subrufa<-subset(species.II$Pelomedusa_subrufa, species.II$Pelomedusa_subrufa$Mean_Temperature>25)
species.II$Pelusios_castaneus<-subset(species.II$Pelusios_castaneus, species.II$Pelusios_castaneus$Mean_Temperature>25)

# Remove data on K. s. hippocrepis with unclear TSD pattern
species.II$Kinosternon_subrubrum<-subset(species.II$Kinosternon_subrubrum, species.II$Kinosternon_subrubrum$Subspecies=="hippocrepis")

############################################################
# Sex-ratio reaction norm calculation

# Create empty vectors for reaction norm data
tpiv.upper<-NULL
tpiv.lower<-NULL
tpiv.upper.97.5<-NULL
tpiv.upper.2.5<-NULL
tpiv.lower.97.5<-NULL
tpiv.lower.2.5<-NULL

# Fit reaction norms for pattern II TSD species and extract pivotal temperature values + credible intervals
for(i in species.II){
  log<-tsd(i, males=i$Males, females=i$Females, temperatures=i$Mean_Temperature, males.freq=F, equation="logistic",
           parameters.initial=c(P_low=min(i$Mean_Temperature)+2, S_low=0.9, P_high=max(i$Mean_Temperature)-2, S_high=-0.9))
  log.p<-tsd_MHmcmc_p(log, accept=T)
  log.p$Min[1]<-log.p$Prior1[1]-5
  log.p$Max[1]<-log.p$Prior1[1]+5
  log.p$Min[3]<-log.p$Prior1[3]-5
  log.p$Max[3]<-log.p$Prior1[3]+5
  log.result<-tsd_MHmcmc(result=log, parametersMCMC=log.p, n.iter=100000, adaptive=T)
  # Commented out code for visualizing the reaction norm fit and lower pivotal temperature prior + posterior
  #plot(log.result, parameters="P_low", xlim=c(log.p$Min[1], log.p$Max[1]), main=i[1,4])
  #plot(log, resultmcmc=log.result, main=i[1,4], xlim=c(18,35), show.PTRT=F, show.observations=T)
  CI<-P_TRT(resultmcmc=log.result, equation="logistic", replicate.CI=100000)
  tpiv.upper<-c(tpiv.upper, as.numeric(CI$P_TRT_quantiles[2,8]))
  tpiv.lower<-c(tpiv.lower, as.numeric(CI$P_TRT_quantiles[2,4]))
  tpiv.upper.97.5<-c(tpiv.upper.97.5, as.numeric(CI$P_TRT_quantiles[3,8]))
  tpiv.upper.2.5<-c(tpiv.upper.2.5, as.numeric(CI$P_TRT_quantiles[1,8]))
  tpiv.lower.97.5<-c(tpiv.lower.97.5, as.numeric(CI$P_TRT_quantiles[3,4]))
  tpiv.lower.2.5<-c(tpiv.lower.2.5, as.numeric(CI$P_TRT_quantiles[1,4]))
}

# Fit reaction norms for pattern Ia TSD species and extract pivotal temperature values + credible intervals
for(i in species.Ia) {
  log<-tsd(i, males=i$Males, females=i$Females, temperatures=i$Mean_Temperature, equation="logistic")
  log.p<-tsd_MHmcmc_p(log, accept=T)
  log.p$Min[1]<-log.p$Prior1[1]-5
  log.p$Max[1]<-log.p$Prior1[1]+5
  log.result<-tsd_MHmcmc(result=log, parametersMCMC=log.p, n.iter=100000, adaptive=T)
  # Commented out code for visualizing the reaction norm fit and pivotal temperature prior + posterior
  #plot(log.result, parameters="P", xlim=c(log.p$Min[1], log.p$Max[1]), main=i[1,4])
  #plot(log, resultmcmc = log.result, main = i[1,4], xlim=c(18,35), show.PTRT=F, show.observations=T, lwd=3, col="gray55")
  CI<-P_TRT(resultmcmc=log.result, equation="logistic", replicate.CI=100000)
  tpiv.upper<-c(tpiv.upper, as.numeric(CI$P_TRT_quantiles[2,4]))
  tpiv.lower<-c(tpiv.lower, as.numeric(CI$P_TRT_quantiles[2,4]))
  tpiv.upper.97.5<-c(tpiv.upper.97.5, as.numeric(CI$P_TRT_quantiles[3,4]))
  tpiv.upper.2.5<-c(tpiv.upper.2.5, as.numeric(CI$P_TRT_quantiles[1,4]))
  tpiv.lower.97.5<-c(tpiv.lower.97.5, as.numeric(CI$P_TRT_quantiles[3,4]))
  tpiv.lower.2.5<-c(tpiv.lower.2.5, as.numeric(CI$P_TRT_quantiles[1,4]))
}

# Put data in a dataframe
tpiv<-data.frame(matrix(NA, nrow = length(species.II)+length(species.Ia), ncol = 10))
colnames(tpiv)<-c("Species", "Pattern", "Upper_Tpiv", "Upper_Tpiv_2.5", "Upper_Tpiv_97.5", "Upper_Tpiv_CI",
                  "Lower_Tpiv", "Lower_Tpiv_2.5", "Lower_Tpiv_97.5", "Lower_Tpiv_CI")
names<-NULL
for(i in species.II) {
  names<-c(names, as.character(i[1,4]))
}
for(i in species.Ia) {
  names<-c(names, as.character(i[1,4]))
}
tpiv$Species<-names
pattern<-c(rep("II",length(species.II)),
       rep("Ia",length(species.Ia)))
tpiv$Pattern<-as.factor(pattern)
tpiv$Upper_Tpiv<-tpiv.upper
tpiv$Upper_Tpiv_2.5<-tpiv.upper.2.5
tpiv$Upper_Tpiv_97.5<-tpiv.upper.97.5
tpiv$Lower_Tpiv<-tpiv.lower
tpiv$Lower_Tpiv_2.5<-tpiv.lower.2.5
tpiv$Lower_Tpiv_97.5<-tpiv.lower.97.5
tpiv$Upper_Tpiv_CI<-as.numeric(tpiv$Upper_Tpiv_97.5-tpiv$Upper_Tpiv_2.5)
tpiv$Lower_Tpiv_CI<-as.numeric(tpiv$Lower_Tpiv_97.5-tpiv$Lower_Tpiv_2.5)
tpiv$Difference<-as.numeric((tpiv$Upper_Tpiv)-(tpiv$Lower_Tpiv))  #Calculate width of reaction norm

# Remove P. expansa due to uncertainty in TSD pattern
tpiv<-tpiv[!tpiv$Species=="Podocnemis_expansa",]

# Create subset of pattern II species
tpiv.II<-subset(tpiv, tpiv$Pattern=="II")

############################################################
# Phylogenetic regressions

# Read in maximum clade credibility chronogram from Thomson et al. (2021) and match with SDM data
# Create separate trees with just pattern II species and all TSD species
thom<-read.tree("Thomson.nwk", tree.names=T)
thom<-drop.tip(thom, setdiff(thom$tip.label, tpiv$Species))
thom.II<-drop.tip(thom, setdiff(thom$tip.label, tpiv.II$Species))
thom<-ladderize(thom)
thom.II<-ladderize(thom.II)

# Match order of Tpiv data to order of tips in Thomson phylogeny
tpiv<-tpiv[tpiv$Species %in% thom$tip.label,]
tpiv<-tpiv[match(thom$tip.label, tpiv$Species),]
tpiv.II<-tpiv.II[tpiv.II$Species %in% thom.II$tip.label,]
tpiv.II<-tpiv.II[match(thom.II$tip.label, tpiv.II$Species),]

# Convert to rrpp dataframe
tpiv<-rrpp.data.frame(tpiv)
tpiv.II<-rrpp.data.frame(tpiv.II)

# Calculate variance-covariance matrix for analysis from the Thomson phylogeny
cov<-vcv(thom)
cov.II<-vcv(thom.II)

# Fit phylogenetic linear regressions with RRPP
model.upper.thom<-lm.rrpp(Difference~Upper_Tpiv, data=tpiv.II, Cov=cov.II, iter=9999)
anova(model.upper.thom)
model.lower.thom<-lm.rrpp(Difference~Lower_Tpiv, data=tpiv.II, Cov=cov.II, iter=9999)
anova(model.lower.thom)

# Repeat analysis with 100 trees from posterior of Thomson et al.
tree<-read.tree("Thomson_Posterior.nwk", tree.names=T)
tree.II<-list()
for(i in 1:length(tree)){
  tree.II[[i]]<-drop.tip(tree[[i]], setdiff(tree[[i]]$tip.label, tpiv.II$Species))
  tree.II[[i]]<-ladderize(tree.II[[i]])
}

up<-list()
low<-list()
for(i in 1:length(tree.II)){
  cov.II.i<-vcv(tree.II[[i]])
  up[[i]]<-lm.rrpp(Difference~Upper_Tpiv, data=tpiv.II, Cov=cov.II.i, iter=9999)
  low[[i]]<-lm.rrpp(Difference~Lower_Tpiv, data=tpiv.II, Cov=cov.II.i, iter=9999)
}

# Extract R2 and p values from posterior tree analysis
up.post<-data.frame(matrix(NA, nrow=101, ncol=2))
low.post<-data.frame(matrix(NA, nrow=101, ncol=2))
colnames(up.post)<-colnames(low.post)<-c("R2","p")
for(i in 1:length(up)){
  up.post[i,1]<-anova(up[[i]])$table[1,4]
  up.post[i,2]<-anova(up[[i]])$table[1,7]
  low.post[i,1]<-anova(low[[i]])$table[1,4]
  low.post[i,2]<-anova(low[[i]])$table[1,7]
}
up.post[101,1]<-anova(model.upper.thom)$table[1,4]
up.post[101,2]<-anova(model.upper.thom)$table[1,7]
low.post[101,1]<-anova(model.lower.thom)$table[1,4]
low.post[101,2]<-anova(model.lower.thom)$table[1,7]

# Check range of correlations and p-values across the 101 dated phylogenies
median(up.post$R2)
HPDinterval(as.mcmc(up.post$R2))
median(up.post$p)
HPDinterval(as.mcmc(up.post$p))
median(low.post$R2)
HPDinterval(as.mcmc(low.post$R2))
median(low.post$p)
HPDinterval(as.mcmc(low.post$p))

# Plot regressions with 95% CI from Thomson phylogeny and regression lines from 100 posterior trees
png(filename = "Thomson Correlations.png", res = 300, width = 2310, height = 2100)

par(mfrow=c(2,1))
par(mar=c(5,6,1,2))
new.up<-seq((min(tpiv.II$Upper_Tpiv)-1), (max(tpiv.II$Upper_Tpiv)+1), by=0.05)
ci.up<-predict(model.upper.thom, newdata=data.frame(Upper_Tpiv=new.up), confidence=0.95)
plot(tpiv.II$Upper_Tpiv, tpiv.II$Difference, xlab=NA, yaxt="n", ylab=NA, pch=16, cex=1.8, cex.axis=1.5, cex.lab=2)
polygon(c(new.up,rev(new.up)),c(ci.up$lcl,rev(ci.up$ucl)), col=alpha("gray80",0.5), border=NA)
for(i in 1:length(up)){
  abline(a=up[[i]]$LM$gls.coefficients[[1]], b=up[[i]]$LM$gls.coefficients[[2]], lwd=3, col=alpha("gray40", alpha=0.1))
}
par(new=TRUE)
plot(tpiv.II$Upper_Tpiv, tpiv.II$Difference, xlab=NA, ylab=NA, pch=16, cex=1.8, cex.axis=1.5, cex.lab=2)
abline(a=model.upper.thom$LM$gls.coefficients[1], b=model.upper.thom$LM$gls.coefficients[2], lwd=4)
title(xlab=expression(paste(Upper~T[piv]," (",degree,"C)")), line=3.5, cex.lab=1.8)
title(ylab=expression(paste("Breadth (",degree,"C)")), line=2.4, cex.lab=1.8)

par(mar=c(5,6,1,2))
new.low<-seq((min(tpiv.II$Lower_Tpiv)-1), (max(tpiv.II$Lower_Tpiv)+1), by=0.05)
ci.low<-predict(model.lower.thom, newdata=data.frame(Lower_Tpiv=new.low), confidence=0.95)
plot(tpiv.II$Lower_Tpiv, tpiv.II$Difference, xlab=NA, ylab=NA, pch=16, cex=1.8, cex.axis=1.5, cex.lab=2)
polygon(c(new.low,rev(new.low)),c(ci.low$lcl,rev(ci.low$ucl)), col=alpha("gray80",0.5), border=NA)
for(i in 1:length(low)){
  abline(a=low[[i]]$LM$gls.coefficients[[1]], b=low[[i]]$LM$gls.coefficients[[2]], lwd=3, col=alpha("gray40", alpha=0.1))
}
par(new=TRUE)
plot(tpiv.II$Lower_Tpiv, tpiv.II$Difference, xlab=NA, ylab=NA, pch=16, cex=1.8, cex.axis=1.5, cex.lab=2)
abline(a=model.lower.thom$LM$gls.coefficients[1], b=model.lower.thom$LM$gls.coefficients[2], lwd=4)
title(xlab=expression(paste(Lower~T[piv]," (",degree,"C)")), line=3.5, cex.lab=1.8)
title(ylab=expression(paste("Breadth (",degree,"C)")), line=2.4, cex.lab=1.8)

dev.off()

# Compare Tpiv distributions between pattern Ia and pattern II
model.thom<-lm.rrpp(Upper_Tpiv~Pattern, data=tpiv, Cov=cov, iter=9999)
anova(model.thom)

# Repeat ANOVA with 100 posterior trees
tree.new<-list()
for(i in 1:length(tree)){
  tree.new[[i]]<-drop.tip(tree[[i]], setdiff(tree[[i]]$tip.label, tpiv$Species))
  tree.new[[i]]<-ladderize(tree.new[[i]])
}

out<-list()
for(i in 1:length(tree.new)){
  cov.i<-vcv(tree.new[[i]])
  out[[i]]<-lm.rrpp(Upper_Tpiv~Pattern, data=tpiv, Cov=cov.i, iter=9999)
}

out.post<-data.frame(matrix(NA, nrow=101, ncol=2))
colnames(out.post)<-c("R2","p")
for(i in 1:length(out)){
  out.post[i,1]<-anova(out[[i]])$table[1,4]
  out.post[i,2]<-anova(out[[i]])$table[1,7]
}
out.post[101,1]<-anova(model.thom)$table[1,4]
out.post[101,2]<-anova(model.thom)$table[1,7]
median(out.post$R2)
HPDinterval(as.mcmc(out.post$R2))
median(out.post$p)
HPDinterval(as.mcmc(out.post$p))

# Plot as density maps rather than boxplots + violin plots to show difference in counts between patterns
tpiv<-data.frame(matrix(NA, nrow = length(species.II)+length(species.Ia), ncol = 10))
colnames(tpiv)<-c("Species", "Pattern", "Upper_Tpiv", "Upper_Tpiv_2.5", "Upper_Tpiv_97.5", "Upper_Tpiv_CI",
                  "Lower_Tpiv", "Lower_Tpiv_2.5", "Lower_Tpiv_97.5", "Lower_Tpiv_CI")
names<-NULL
for(i in species.II) {
  names<-c(names, as.character(i[1,4]))
}
for(i in species.Ia) {
  names<-c(names, as.character(i[1,4]))
}
tpiv$Species<-names
pattern<-c(rep("II",length(species.II)),
           rep("Ia",length(species.Ia)))
tpiv$Pattern<-as.factor(pattern)
tpiv$Upper_Tpiv<-tpiv.upper
tpiv$Upper_Tpiv_2.5<-tpiv.upper.2.5
tpiv$Upper_Tpiv_97.5<-tpiv.upper.97.5
tpiv$Lower_Tpiv<-tpiv.lower
tpiv$Lower_Tpiv_2.5<-tpiv.lower.2.5
tpiv$Lower_Tpiv_97.5<-tpiv.lower.97.5
tpiv$Upper_Tpiv_CI<-as.numeric(tpiv$Upper_Tpiv_97.5-tpiv$Upper_Tpiv_2.5)
tpiv$Lower_Tpiv_CI<-as.numeric(tpiv$Lower_Tpiv_97.5-tpiv$Lower_Tpiv_2.5)
tpiv$Difference<-as.numeric((tpiv$Upper_Tpiv)-(tpiv$Lower_Tpiv))
tpiv<-tpiv[!tpiv$Species=="Podocnemis_expansa",] #Remove due to uncertain pattern

png(filename = "Thomson Tpiv Comparison Densities.png", res = 300, width = 2310, height = 2100)

ggplot(tpiv, aes(x=Upper_Tpiv, y=after_stat(count), fill=Pattern)) +
  geom_density(color=NA) +
  geom_boxplot(aes(x=Upper_Tpiv, fill=Pattern, y=-3), inherit.aes=F, show.legend=F,
               width=4, lwd=1.2, outlier.size=3, position=position_dodge(width=5)) +
  theme_classic() +
  scale_fill_manual(values=alpha(c("dodgerblue4","steelblue1"),0.8), name="Pattern") +
  xlab(expression(paste(T[piv]," (",degree,"C)"))) +
  ylab("# Species") +
  xlim(25, 33) +
  theme(axis.title.y=element_text(size=24, color="black"),
        axis.text=element_text(size=20, color="black"),
        axis.title.x=element_text(size=24, color="black"),
        legend.position=c(0.88,0.88),
        legend.key.size=unit(8, "mm"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=22, color="black"),
        plot.margin = unit(c(0.5,10,0.5,0.5),"mm"))

dev.off()

# Repeat comparison without outliers and C. rossignonii
tpiv<-tpiv[!tpiv$Species=="Kinosternon_stejnegeri",]
tpiv<-tpiv[!tpiv$Species=="Pelomedusa_subrufa",]
tpiv<-tpiv[!tpiv$Species=="Pelusios_castaneus",]
tpiv<-tpiv[!tpiv$Species=="Chelydra_rossignonii",]
tpiv<-tpiv[tpiv$Species %in% thom$tip.label,]
thom<-drop.tip(thom, setdiff(thom$tip.label, tpiv$Species))
tpiv<-tpiv[match(thom$tip.label, tpiv$Species),]
tpiv<-rrpp.data.frame(tpiv)
cov<-vcv(thom)

model.thom<-lm.rrpp(Upper_Tpiv~Pattern, data=tpiv, Cov=cov, iter=9999)
anova(model.thom)

tree.new<-list()
for(i in 1:length(tree)){
  tree.new[[i]]<-drop.tip(tree[[i]], setdiff(tree[[i]]$tip.label, tpiv$Species))
  tree.new[[i]]<-ladderize(tree.new[[i]])
}

out<-list()
for(i in 1:length(tree.new)){
  cov.i<-vcv(tree.new[[i]])
  out[[i]]<-lm.rrpp(Upper_Tpiv~Pattern, data=tpiv, Cov=cov.i, iter=9999)
}

out.post<-data.frame(matrix(NA, nrow=101, ncol=2))
colnames(out.post)<-c("R2","p")
for(i in 1:length(out)){
  out.post[i,1]<-anova(out[[i]])$table[1,4]
  out.post[i,2]<-anova(out[[i]])$table[1,7]
}
out.post[101,1]<-anova(model.thom)$table[1,4]
out.post[101,2]<-anova(model.thom)$table[1,7]
median(out.post$R2)
HPDinterval(as.mcmc(out.post$R2))
median(out.post$p)
HPDinterval(as.mcmc(out.post$p))

############################################################
# Calculate correlation between (upper) Tpiv and liabilities simulated under threshold model

# Recreate Tpiv dataframe
tpiv<-data.frame(matrix(NA, nrow = length(species.II)+length(species.Ia), ncol = 10))
colnames(tpiv)<-c("Species", "Pattern", "Upper_Tpiv", "Upper_Tpiv_2.5", "Upper_Tpiv_97.5", "Upper_Tpiv_CI",
                  "Lower_Tpiv", "Lower_Tpiv_2.5", "Lower_Tpiv_97.5", "Lower_Tpiv_CI")
names<-NULL
for(i in species.II) {
  names<-c(names, as.character(i[1,4]))
}
for(i in species.Ia) {
  names<-c(names, as.character(i[1,4]))
}
tpiv$Species<-names
pattern<-c(rep("II",length(species.II)),
           rep("Ia",length(species.Ia)))
tpiv$Pattern<-as.factor(pattern)
tpiv$Upper_Tpiv<-tpiv.upper
tpiv$Upper_Tpiv_2.5<-tpiv.upper.2.5
tpiv$Upper_Tpiv_97.5<-tpiv.upper.97.5
tpiv$Lower_Tpiv<-tpiv.lower
tpiv$Lower_Tpiv_2.5<-tpiv.lower.2.5
tpiv$Lower_Tpiv_97.5<-tpiv.lower.97.5
tpiv$Upper_Tpiv_CI<-as.numeric(tpiv$Upper_Tpiv_97.5-tpiv$Upper_Tpiv_2.5)
tpiv$Lower_Tpiv_CI<-as.numeric(tpiv$Lower_Tpiv_97.5-tpiv$Lower_Tpiv_2.5)
tpiv$Difference<-as.numeric((tpiv$Upper_Tpiv)-(tpiv$Lower_Tpiv))
tpiv<-tpiv[!tpiv$Species=="Podocnemis_expansa",] #Remove due to uncertain pattern

# Read in phylogeny
thom<-read.tree("Thomson.nwk", tree.names=T)
thom<-drop.tip(thom, setdiff(thom$tip.label, tpiv$Species))
thom<-ladderize(thom)

dat<-data.frame(cbind(tpiv$Pattern,tpiv$Upper_Tpiv))
rownames(dat)<-tpiv$Species
colnames(dat)<-c("Pattern","Tpiv")
dat<-dat[match(thom$tip.label,rownames(dat)),]
dat$Pattern<-as.factor(dat$Pattern)
dat$Tpiv<-as.numeric(dat$Tpiv)

cov.out.full<-threshBayes(thom, dat, c("discrete","continuous"), ngen=10.1e6, control=list(sample=100,quiet=T), plot=F)
plot(cov.out.full)
cov.out.full
HPDinterval(as.mcmc(cov.out.full$par[(1e5/100)+1:nrow(cov.out.full$par),"r"]))

# Repeat analysis without outliers and Chelydra rossignonii
thom<-drop.tip(thom, c("Kinosternon_stejnegeri","Pelomedusa_subrufa","Pelusios_castaneus"))
thom<-drop.tip(thom, "Chelydra_rossignonii")
dat<-data.frame(cbind(tpiv$Pattern,tpiv$Upper_Tpiv))
rownames(dat)<-tpiv$Species
colnames(dat)<-c("Pattern","Tpiv")
dat<-dat[match(thom$tip.label,rownames(dat)),]
dat$Pattern<-as.factor(dat$Pattern)
dat$Tpiv<-as.numeric(dat$Tpiv)

cov.out<-threshBayes(thom, dat, c("discrete","continuous"), ngen=10.1e6, control=list(sample=100,quiet=T), plot=F)
plot(cov.out)
cov.out
HPDinterval(as.mcmc(cov.out$par[(1e5/100)+1:nrow(cov.out$par),"r"]))