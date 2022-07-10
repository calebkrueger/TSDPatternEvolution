# The Macroevolution of TSD Patterns

# Load the necessary packages
library("corHMM")
library("phytools")
library("geiger")
library("coda")
library("parallel")
library("plotrix")

############################################################
# Fit hidden Markov models (HMMs)

# Read in file with sex-determining mechanism classifications from ROSIE
data<-read.csv("https://raw.githubusercontent.com/calebkrueger/ROSIE/main/SDM.csv")
data<-subset(data, SDM_Confidence==1)
data[data$Species=="Kinosternon_subrubrum",8:9]<-NA
data<-subset(data, Type_Confidence==1)

#Change format for corHMM functions
traits<-data.frame(data$Species,as.factor(data$Type))
colnames(traits)<-c("Species","SDM")

# Read in phylogeny from Thomson et al. (2021) and match data order with SDM data
thom<-read.tree("Thomson.nwk",tree.names=T)
thom<-drop.tip(thom, setdiff(thom$tip.label, traits$Species))
thom<-ladderize(thom)
traits<-traits[traits$Species %in% thom$tip.label,]
traits<-traits[match(thom$tip.label,traits$Species),]

# Create transition rate matrices to fit
leg.and.rate<-getStateMat4Dat(traits)

base<-dropStateMatPars(leg.and.rate$rate.mat, c(2,3,5,6,9,12))
model1<-equateStateMatPars(base, c(1,2,3,4,5,6))
model2<-equateStateMatPars(base, list(c(1,2),c(3,4,5,6)))
model3<-equateStateMatPars(base, list(c(1,2),c(3,4),c(5,6)))
model4<-equateStateMatPars(base, list(c(3,4),c(5,6)))
model5<-equateStateMatPars(base, c(3,4,5,6))

# Fit models with several potential transition rate matrices
fitER<-corHMM(thom, traits, 1, model="ER", nstarts=10, root.p="NULL")
fitSYM<-corHMM(thom, traits, 1, model="SYM", nstarts=10, root.p="NULL")
fitARD<-corHMM(thom, traits, 1, model="ARD", nstarts=10, root.p="NULL")
fitmodel1<-corHMM(thom, traits, 1, model1, nstarts=10, root.p="NULL")
fitmodel2<-corHMM(thom, traits, 1, model2, nstarts=10, root.p="NULL")
fitmodel3<-corHMM(thom, traits, 1, model3, nstarts=10, root.p="NULL")
fitmodel4<-corHMM(thom, traits, 1, model4, nstarts=10, root.p="NULL")
fitmodel5<-corHMM(thom, traits, 1, model5, nstarts=10, root.p="NULL")

# Compare model fit via AIC and AIC weights
AIC<-setNames(c(fitER$AIC, fitSYM$AIC, fitARD$AIC, fitmodel1$AIC, fitmodel2$AIC, fitmodel3$AIC, fitmodel4$AIC, fitmodel5$AIC),
              c("ER","SYM","ARD","model1","model2","model3","model4","model5"))
AIC
aic.w(AIC)

# Visualize ancestral state estimates
colnames(fitmodel5$states)<-c("Ia","II","XY","ZW")
colnames(fitmodel5$tip.states)<-c("Ia","II","XY","ZW")

cols<-setNames(c("dodgerblue4","steelblue1","gray65","gray80")[1:length(unique(colnames(fitmodel5$tip.states)))],sort(unique(colnames(fitmodel5$tip.states))))
h<-max(nodeHeights(thom))
offset.factor<-1.03
plotTree(geiger::rescale(thom,model="depth",depth=offset.factor*h),
         color="transparent",ftype="i",type="fan",fsize=0.8,cex=0.4,lwd=4,mar=c(0,0,0,0), 
         maxY=110.75)
par(fg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
draw.circle(0,0,max(node.depth.edgelength(thom)),col=NA,border="gray85")
draw.circle(0,0,max(node.depth.edgelength(thom))-50,col=NA,border="gray85")
draw.circle(0,0,max(node.depth.edgelength(thom))-100,col=NA,border="gray85")
draw.circle(0,0,max(node.depth.edgelength(thom))-150,col=NA,border="gray85")
draw.circle(0,0,max(node.depth.edgelength(thom))-200,col=NA,border="gray85")
plotTree(thom, type="fan", colors=cols, fsize=0.8,cex=0.4,ftype="i",lwd=2,mar=c(0,0,0,0), add=T, 
         xlim=obj$x.lim,ylim=obj$y.lim, maxY=110.75)
tiplabels(pie=fitmodel5$tip.states, piecol=cols, cex=0.43)
nodelabels(pie=fitmodel5$states, piecol=cols, cex=0.43)
nodelabels(pie=fitmodel5$states[1,], node=110, piecol=cols, cex=1.2)
par(fg="black")
obj<-axis(1,pos=3,at=seq(max(nodeHeights(thom)),max(nodeHeights(thom))-100,by=-50),cex.axis=1,labels=F,lwd=1.8,tcl=-0.2)
text(obj,rep(-7,length(obj)),max(nodeHeights(thom))-obj,cex=0.8)
add.simmap.legend(colors=cols[1:4], shape="circle", prompt=F, x=max(nodeHeights(thom))+110, y=22, fsize=1.2)

# Use best-fit model (model5) to create HMM
Qs<-list(model5, model5)
cat.mat<-getRateCatMat(2)
cat.mat<-equateStateMatPars(cat.mat, 1:2)
mat<-getFullMat(Qs, cat.mat)

# For Windows users, this code is necessary to run random restarts in parallel
# Mac users can simply set the 'n.cores' argument in the corHMM function
setCores<-round(detectCores() * 0.85)
cl<-parallel::makeCluster(getOption("cl.cores", setCores))
cl.pkg<-parallel::clusterEvalQ(cl, library(corHMM))
parallel::clusterExport(cl, "thom")
parallel::clusterExport(cl, "traits")
parallel::clusterExport(cl, "mat")
NRUNS <- 10
fit <- parallel::parLapply(cl = cl, 1:NRUNS, function(i) {
  corHMM(thom, traits, 2, mat, nstarts=10, root.p="NULL")
})
stopCluster(cl)

# Extract run with lowest AIC
fithmm<-fit[[1]]
for(i in 1:length(fit)){
  if(fit[[i]][[2]]<fithmm[[2]]){
    fithmm<-fit[[i]]
  }else{
  }
}

# Visualize ancestral state estimates
# NOTE: root node indicates state not observed in any extant taxa
plot(thom, type="fan")
nodelabels(pie=fithmm$states, cex=0.4)

# Fit additional HMMs that avoid assigning unobserved state to root node
mat.noII2Ia<-dropStateMatPars(mat, 1)
mat.noIa2II<-dropStateMatPars(mat, 2)
mat.noTSD<-dropStateMatPars(mat, c(1,2))
mat.noTSD.eqGSD<-equateStateMatPars(mat.noTSD, c(2,5))

setCores<-round(detectCores() * 0.85)
cl<-parallel::makeCluster(getOption("cl.cores", setCores))
cl.pkg<-parallel::clusterEvalQ(cl, library(corHMM))
parallel::clusterExport(cl, "thom")
parallel::clusterExport(cl, "traits")
parallel::clusterExport(cl, "mat.noIa2II")
NRUNS <- 10
fit <- parallel::parLapply(cl = cl, 1:NRUNS, function(i) {
  corHMM(thom, traits, 2, mat.noIa2II, nstarts=10, root.p="NULL")
})
stopCluster(cl)
fitnoIa2II<-fit[[1]]
for(i in 1:length(fit)){
  if(fit[[i]][[2]]<fitnoIa2II[[2]]){
    fitnoIa2II<-fit[[i]]
  }else{
  }
}

setCores<-round(detectCores() * 0.85)
cl<-parallel::makeCluster(getOption("cl.cores", setCores))
cl.pkg<-parallel::clusterEvalQ(cl, library(corHMM))
parallel::clusterExport(cl, "thom")
parallel::clusterExport(cl, "traits")
parallel::clusterExport(cl, "mat.noII2Ia")
NRUNS <- 10
fit <- parallel::parLapply(cl = cl, 1:NRUNS, function(i) {
  corHMM(thom, traits, 2, mat.noII2Ia, nstarts=10, root.p="NULL")
})
stopCluster(cl)
fitnoII2Ia<-fit[[1]]
for(i in 1:length(fit)){
  if(fit[[i]][[2]]<fitnoII2Ia[[2]]){
    fitnoII2Ia<-fit[[i]]
  }else{
  }
}

setCores<-round(detectCores() * 0.85)
cl<-parallel::makeCluster(getOption("cl.cores", setCores))
cl.pkg<-parallel::clusterEvalQ(cl, library(corHMM))
parallel::clusterExport(cl, "thom")
parallel::clusterExport(cl, "traits")
parallel::clusterExport(cl, "mat.noTSD")
NRUNS <- 10
fit <- parallel::parLapply(cl = cl, 1:NRUNS, function(i) {
  corHMM(thom, traits, 2, mat.noTSD, nstarts=10, root.p="NULL")
})
stopCluster(cl)
fitnoTSD<-fit[[1]]
for(i in 1:length(fit)){
  if(fit[[i]][[2]]<fitnoTSD[[2]]){
    fitnoTSD<-fit[[i]]
  }else{
  }
}

setCores<-round(detectCores() * 0.85)
cl<-parallel::makeCluster(getOption("cl.cores", setCores))
cl.pkg<-parallel::clusterEvalQ(cl, library(corHMM))
parallel::clusterExport(cl, "thom")
parallel::clusterExport(cl, "traits")
parallel::clusterExport(cl, "mat.noTSD.eqGSD")
NRUNS <- 10
fit <- parallel::parLapply(cl = cl, 1:NRUNS, function(i) {
  corHMM(thom, traits, 2, mat.noTSD.eqGSD, nstarts=10, root.p="NULL")
})
stopCluster(cl)
fitnoTSD.eqGSD<-fit[[1]]
for(i in 1:length(fit)){
  if(fit[[i]][[2]]<fitnoTSD.eqGSD[[2]]){
    fitnoTSD.eqGSD<-fit[[i]]
  }else{
  }
}

# Compare model fit via AIC and AIC weights
# Including model5 from previous comparisons to determine if HMM fits better than standard Mk model
AIC<-setNames(c(fitnoII2Ia$AIC, fitnoIa2II$AIC, fitnoTSD$AIC, fitnoTSD.eqGSD$AIC, fitmodel5$AIC),
              c("noII2Ia","noIa2II","noTSD","noTSD.eqGSD","model5"))
AIC
aic.w(AIC)

# Visualize ancestral state estimates from best-fit model (fitnoTSD.eqGSD)
plot(thom, type="fan")
nodelabels(pie=fitnoTSD.eqGSD$states, cex=0.4)

# Iterate fitting best-fit model over 100 trees from posterior
tree<-read.tree("Thomson_Posterior.nwk", tree.names=T)
for(i in 1:length(tree)){
  tree[[i]]<-drop.tip(tree[[i]], setdiff(tree[[i]]$tip.label, traits$Species))
  tree[[i]]<-ladderize(tree[[i]])
}
out<-list()
for(i in 1:length(tree)){
  traits.i<-traits[match(tree[[i]]$tip.label,traits$Species),]
  setCores<-round(detectCores() * 0.85)
  cl<-parallel::makeCluster(getOption("cl.cores", setCores))
  cl.pkg<-parallel::clusterEvalQ(cl, library(corHMM))
  parallel::clusterExport(cl, "i")
  parallel::clusterExport(cl, "tree")
  parallel::clusterExport(cl, "traits.i")
  parallel::clusterExport(cl, "mat.noTSD.eqGSD")
  NRUNS <- 10
  fit <- parallel::parLapply(cl = cl, 1:NRUNS, function(i) {
    corHMM(tree[[i]], traits.i, 2, mat.noTSD.eqGSD, nstarts=10, root.p="NULL")
  })
  stopCluster(cl)
  fit.temp<-fit[[1]]
  for(j in 1:length(fit)){
    if(fit[[j]][[2]]<fit.temp[[2]]){
      fit.temp<-fit[[j]]
    }else{
    }
  }
  out[[i]]<-fit.temp
}

# Add output from using Thomson phylogeny (maximum clade credibility chronogram) to list
out[[101]]<-fitnoTSD.eqGSD

# Extract ancestral states for nodes down to most common ancestor of each family
nodes<-c(110,201,213,216,214,202,111,112,113,122,180,197,181,182,184,185,123,152,124,125,137)
node.states<-replicate(length(out), data.frame(matrix(NA, nrow=length(nodes), ncol=4)), simplify=F)
for(i in 1:length(node.states)){
  rownames(node.states[[i]])<-nodes
  colnames(node.states[[i]])<-c("Ia","II","XY","ZW")
}

for(i in 1:length(out)){
  for(j in unique(nodes)){
    node.states[[i]][rownames(node.states[[i]])==j,1]<-out[[i]]$states[j-109,1]+out[[i]]$states[j-109,5]
    node.states[[i]][rownames(node.states[[i]])==j,2]<-out[[i]]$states[j-109,2]+out[[i]]$states[j-109,6]
    node.states[[i]][rownames(node.states[[i]])==j,3]<-out[[i]]$states[j-109,3]+out[[i]]$states[j-109,7]
    node.states[[i]][rownames(node.states[[i]])==j,4]<-out[[i]]$states[j-109,4]+out[[i]]$states[j-109,8]
  }
}

anc.states<-data.frame(matrix(NA, nrow=length(nodes), ncol=12))
rownames(anc.states)<-nodes
colnames(anc.states)<-c("Ia", "Ia Low", "Ia High", "II", "II Low", "II High", "XY", "XY Low", "XY High", "ZW", "ZW Low", "ZW High")

Ia<-NULL
II<-NULL
XY<-NULL
ZW<-NULL
for(i in unique(nodes)){
  for(j in 1:length(node.states)){
    Ia<-c(Ia, node.states[[j]][rownames(node.states[[j]])==i,1])
    II<-c(II, node.states[[j]][rownames(node.states[[j]])==i,2])
    XY<-c(XY, node.states[[j]][rownames(node.states[[j]])==i,3])
    ZW<-c(ZW, node.states[[j]][rownames(node.states[[j]])==i,4])
  }
  anc.states[rownames(anc.states)==i,"Ia"]<-median(Ia)
  anc.states[rownames(anc.states)==i,"Ia Low"]<-HPDinterval(as.mcmc(Ia))[1]
  anc.states[rownames(anc.states)==i,"Ia High"]<-HPDinterval(as.mcmc(Ia))[2]
  anc.states[rownames(anc.states)==i,"II"]<-median(II)
  anc.states[rownames(anc.states)==i,"II Low"]<-HPDinterval(as.mcmc(II))[1]
  anc.states[rownames(anc.states)==i,"II High"]<-HPDinterval(as.mcmc(II))[2]
  anc.states[rownames(anc.states)==i,"XY"]<-median(XY)
  anc.states[rownames(anc.states)==i,"XY Low"]<-HPDinterval(as.mcmc(XY))[1]
  anc.states[rownames(anc.states)==i,"XY High"]<-HPDinterval(as.mcmc(XY))[2]
  anc.states[rownames(anc.states)==i,"ZW"]<-median(ZW)
  anc.states[rownames(anc.states)==i,"ZW Low"]<-HPDinterval(as.mcmc(ZW))[1]
  anc.states[rownames(anc.states)==i,"ZW High"]<-HPDinterval(as.mcmc(ZW))[2]
  Ia<-NULL
  II<-NULL
  XY<-NULL
  ZW<-NULL
}

# View median probability and 95% interval for state at ancestral node
anc.states[1,]

tip.states<-out[[1]]$tip.states[,1:4]
colnames(tip.states)<-c("Ia","II","XY","ZW")

# Plot median ancestral state values on Thomson phylogeny
thom<-read.tree("Thomson.nwk", tree.names=T)
thom<-drop.tip(thom, setdiff(thom$tip.label, out[[1]]$phy$tip.label))
thom<-ladderize(thom)
tip.states<-tip.states[match(thom$tip.label, rownames(tip.states)),]

png(filename = "ASE_Thomson_HMM.png", res = 300, width = 180, height = 170, units="mm")
par(mar = c(0,0,0,0), mfrow = c(1,1))

cols<-setNames(c("dodgerblue4","steelblue1","gray65","gray80")[1:length(unique(colnames(tip.states)))],sort(unique(colnames(tip.states))))
h<-max(nodeHeights(thom))
offset.factor<-1.03
plotTree(geiger::rescale(thom,model="depth",depth=offset.factor*h),
         color="transparent",ftype="i",type="fan",fsize=0.55,cex=1,lwd=1,mar=c(0,0,0,0), 
         maxY=110.75)
par(fg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
draw.circle(0,0,max(node.depth.edgelength(thom)),col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-50,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-100,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-150,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-200,col=NA,border="gray90")
plotTree(thom, type="fan", colors=cols, fsize=0.55,cex=1,ftype="i",lwd=1,mar=c(0,0,0,0), add=T, 
         xlim=obj$x.lim,ylim=obj$y.lim, maxY=110.75)
tiplabels(pie=tip.states, piecol=cols, cex=0.4)
for(i in unique(nodes)){
  if(i==110){
    nodelabels(pie=anc.states[rownames(anc.states)==i,c("Ia","II","XY","ZW")], node=i, piecol=cols, cex=1.2)
  }else{
    nodelabels(pie=anc.states[rownames(anc.states)==i,c("Ia","II","XY","ZW")], node=i, piecol=cols, cex=0.4) 
  }
}
par(fg="black")
obj<-axis(1,pos=3,at=seq(max(nodeHeights(thom)),max(nodeHeights(thom))-100,by=-50),cex.axis=1,labels=F,lwd=1.2,tcl=-0.2)
text(obj,rep(-7,length(obj)),max(nodeHeights(thom))-obj,cex=0.6)
add.simmap.legend(colors=cols[1:4], shape="circle", prompt=F, x=max(nodeHeights(thom))+115, y=18, fsize=0.75)

dev.off()

############################################################
# Fit threshold model for Testudines only

# Read in file with sex-determining mechanism classifications from ROSIE
data<-read.csv("https://raw.githubusercontent.com/calebkrueger/ROSIE/main/SDM.csv")
data<-subset(data, SDM_Confidence==1)
data[data$Species=="Kinosternon_subrubrum",8:9]<-NA
sdms<-as.factor(t(data$Type)); names(sdms)<-data$Species
PP<-to.matrix(sdms, seq=c("Ia","II","XY","ZW"))
for(i in unique(rownames(PP))){
  if(data[data$Species==i,9]==0 | is.na(data[data$Species==i,9])){
    if(data[data$Species==i,6]=="GSD"){
      PP[i,3:4]<-c(0.5,0.5)
    }else{
      PP[i,1:2]<-c(0.5,0.5)
    }
  }
}
PP.thresh<-subset(PP[,1:2], PP[,1]==1 | PP[,2]==1 | PP[,1]==0.5)

# Read in the file containing a random sample of trees from the posterior of Thomson et al.
tree<-read.tree("Thomson_Posterior.nwk", tree.names=T)
for(i in 1:length(tree)){
  tree[[i]]<-drop.tip(tree[[i]], setdiff(tree[[i]]$tip.label, rownames(PP.thresh)))
  tree[[i]]<-ladderize(tree[[i]])
}

# Read in the Thomson phylogeny
thom<-read.tree("Thomson.nwk", tree.names=T)
thom<-drop.tip(thom, setdiff(thom$tip.label, rownames(PP.thresh)))
thom<-ladderize(thom)

# Perform ancestral state reconstruction under the threshold model with each tree
out<-list()
for(i in 1:length(tree)){
  PP.thresh.i<-PP.thresh[match(tree[[i]]$tip.label,rownames(PP.thresh)),]
  out[[i]]<-ancThresh(tree[[i]],PP.thresh.i,ngen=1e6,model="BM",control=list(print=FALSE))
}
PP.thresh<-PP.thresh[match(thom$tip.label, rownames(PP.thresh)),]
out[[101]]<-ancThresh(thom,PP.thresh,ngen=1e6,model="BM",control=list(print=FALSE))

# Extract ancestral states for select nodes
nodes<-c(99,188,194,189,100,101,167,184,168,169,171,172,102,140,103,125,104)
node.states<-replicate(length(out), data.frame(matrix(NA, nrow=length(nodes), ncol=length(colnames(PP.thresh)))), simplify=F)
for(i in 1:length(node.states)){
  row.names(node.states[[i]])<-nodes
  colnames(node.states[[i]])<-c("Ia","II")
}

for(i in 1:length(out)){
  for(j in unique(nodes)){
    node.states[[i]][row.names(node.states[[i]])==j,1]<-out[[i]]$ace[row.names(out[[i]]$ace)==j,1]
    node.states[[i]][row.names(node.states[[i]])==j,2]<-out[[i]]$ace[row.names(out[[i]]$ace)==j,2]
  }
}

anc.states<-data.frame(matrix(NA, nrow=length(nodes), ncol=3*length(colnames(PP.thresh))))
rownames(anc.states)<-nodes
colnames(anc.states)<-c("Ia", "Ia Low", "Ia High", "II", "II Low", "II High")

Ia<-NULL
II<-NULL
for(i in unique(nodes)){
  for(j in 1:length(node.states)){
    Ia<-c(Ia, node.states[[j]][rownames(node.states[[j]])==i,1])
    II<-c(II, node.states[[j]][rownames(node.states[[j]])==i,2])
  }
  anc.states[rownames(anc.states)==i,"Ia"]<-median(Ia)
  anc.states[rownames(anc.states)==i,"Ia Low"]<-HPDinterval(as.mcmc(Ia))[1]
  anc.states[rownames(anc.states)==i,"Ia High"]<-HPDinterval(as.mcmc(Ia))[2]
  anc.states[rownames(anc.states)==i,"II"]<-median(II)
  anc.states[rownames(anc.states)==i,"II Low"]<-HPDinterval(as.mcmc(II))[1]
  anc.states[rownames(anc.states)==i,"II High"]<-HPDinterval(as.mcmc(II))[2]
  Ia<-NULL
  II<-NULL
}

# View median probability and 95% interval for state at ancestral node
anc.states[1,]

# Extract reconstructed states for tips, including those with TSD but unknown pattern
tips<-list()
for(i in 1:length(out)){
  tips[[i]]<-plot(out[[i]])[1:length(out[[i]]$tree$tip.label),]
}

tip.states<-data.frame(matrix(NA, nrow=length(out[[1]]$tree$tip.label), ncol=length(colnames(PP.thresh))))
rownames(tip.states)<-out[[1]]$tree$tip.label
colnames(tip.states)<-c("Ia","II")

tips.Ia<-NULL
tips.II<-NULL
for(i in unique(out[[1]]$tree$tip.label)){
  for(j in 1:length(tips)){
    tips.Ia<-c(tips.Ia, tips[[j]][rownames(tips[[j]])==i,1])
    tips.II<-c(tips.II, tips[[j]][rownames(tips[[j]])==i,2])
  }
  tip.states[rownames(tip.states)==i,"Ia"]<-median(tips.Ia)
  tip.states[rownames(tip.states)==i,"II"]<-median(tips.II)
  tips.Ia<-NULL
  tips.II<-NULL
}

# Plot median ancestral state values on Thomson phylogeny
thom<-read.tree("Thomson.nwk", tree.names=T)
thom<-drop.tip(thom, setdiff(thom$tip.label, out[[1]]$tree$tip.label))
thom<-ladderize(thom)
tip.states<-tip.states[match(thom$tip.label, rownames(tip.states)),]

png(filename = "ASE_Thomson_Thresh.png", res = 450, width = 180, height = 170, units="mm")
par(mar = c(0,0,0,0), mfrow = c(1,1))

cols<-setNames(c("dodgerblue4","steelblue1")[1:length(unique(colnames(PP.thresh)))],sort(unique(colnames(PP.thresh))))
h<-max(nodeHeights(thom))
offset.factor<-1.03
plotTree(geiger::rescale(thom,model="depth",depth=offset.factor*h),
         color="transparent",ftype="i",type="fan",fsize=0.55,cex=1,lwd=1, 
         maxY=99.5)
par(fg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
draw.circle(0,0,max(node.depth.edgelength(thom)),col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-50,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-100,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-150,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-200,col=NA,border="gray90")
plotTree(thom, type="fan", colors=cols, fsize=0.55,cex=1,ftype="i",lwd=1, add=T, 
         xlim=obj$x.lim,ylim=obj$y.lim, maxY=99.5)
tiplabels(pie=tip.states, piecol=cols, cex=0.43)
for(i in unique(nodes)){
  if(i==99){
    nodelabels(pie=anc.states[rownames(anc.states)==i,c("Ia","II")], node=i, piecol=cols, cex=1.2)
  }else{
    nodelabels(pie=anc.states[rownames(anc.states)==i,c("Ia","II")], node=i, piecol=cols, cex=0.43) 
  }
}
par(fg="black")
obj<-axis(1,pos=3,at=seq(max(nodeHeights(thom)),max(nodeHeights(thom))-100,by=-50),cex.axis=1,labels=F,lwd=1.2,tcl=-0.2)
text(obj,rep(-7,length(obj)),max(nodeHeights(thom))-obj,cex=0.6)
add.simmap.legend(colors=cols[1:2], shape="circle", prompt=F, x=max(nodeHeights(thom))+115, y=6, fsize=0.75)

dev.off()

############################################################
# Fit threshold model for Testudines and Crocodilia

# Read in file with sex-determining mechanism classifications from ROSIE
data<-read.csv("https://raw.githubusercontent.com/calebkrueger/ROSIE/main/SDM.csv")
data<-subset(data, SDM_Confidence==1)
data[data$Species=="Kinosternon_subrubrum",8:9]<-NA

# Read in SDM data for crocodilians
croc.sdm<-read.csv("Crocodilia_SDM.csv")
data<-data[,1:9]
data<-rbind(data, croc.sdm)
sdms<-as.factor(t(data$Type)); names(sdms)<-data$Species
PP<-to.matrix(sdms, seq=c("Ia","Ib","II","XY","ZW"))
for(i in unique(rownames(PP))){
  if(data[data$Species==i,9]==0 | is.na(data[data$Species==i,9])){
    if(data[data$Species==i,6]=="GSD"){
      PP[i,4:5]<-c(0.5,0.5)
    }else{
      if(data[data$Species==i,1]=="Testudines"){
        PP[i,c(1,3)]<-c(0.5,0.5)
      }else{
        PP[i,2:3]<-c(0.5,0.5)
      }
    }
  }
}
PP.thresh<-subset(PP[,1:3], PP[,1]==1 | PP[,2]==1 | PP[,3]==1 | PP[,1]==0.5 | PP[,2]==0.5)

# Read in the file containing a random sample of trees from the posterior of Thomson et al.
tree<-read.tree("Thomson_Posterior.nwk", tree.names=T)

# Add crocodilian data from TimeTree.org to each tree
# Time of divergence (261.37mya) from TimeTree
crocs<-read.tree("crocodilia_species.nwk", tree.names=T)
for(i in 1:length(tree)){
  tree[[i]]$root.edge<-261.37
  crocs$root.edge<-261.37-max(node.depth.edgelength(crocs))
  tree[[i]]<-bind.tree(tree[[i]], crocs, where="root", position=261.37-max(node.depth.edgelength(tree[[i]])))
}

thom<-read.tree("Thomson.nwk", tree.names=T)
thom<-drop.tip(thom, tip="Alligator_mississippiensis")
thom<-drop.tip(thom, tip="Gallus_gallus")
thom$root.edge<-261.37
crocs$root.edge<-261.37-max(node.depth.edgelength(crocs))
thom<-bind.tree(thom, crocs, where="root", position=261.37-max(node.depth.edgelength(thom)))

# Remove species without data
for(i in 1:length(tree)){
  tree[[i]]<-drop.tip(tree[[i]], setdiff(tree[[i]]$tip.label, rownames(PP.thresh)))
  tree[[i]]<-ladderize(tree[[i]])
}

thom<-drop.tip(thom, setdiff(thom$tip.label, rownames(PP.thresh)))
thom<-ladderize(thom)

# Perform ancestral state reconstruction under the threshold model with each tree
out<-list()
for(i in 1:length(tree)){
  PP.thresh.i<-PP.thresh[match(tree[[i]]$tip.label,rownames(PP.thresh)),]
  out[[i]]<-ancThresh(tree[[i]],PP.thresh.i,sequence=c("Ia","II","Ib"),ngen=1e6,model="BM",control=list(print=FALSE))
}
PP.thresh<-PP.thresh[match(thom$tip.label, rownames(PP.thresh)),]
out[[101]]<-ancThresh(thom,PP.thresh,sequence=c("Ia","II","Ib"),ngen=1e6,model="BM",control=list(print=FALSE))

# Extract ancestral states for select nodes within Testudines
nodes<-c(118,207,213,208,119,120,186,203,187,188,190,191,121,159,122,144,123)
node.states<-replicate(length(out), data.frame(matrix(NA, nrow=length(nodes), ncol=length(colnames(PP.thresh)))), simplify=F)
for(i in 1:length(node.states)){
  row.names(node.states[[i]])<-nodes
  colnames(node.states[[i]])<-c("Ia","Ib","II")
}

for(i in 1:length(out)){
  for(j in unique(nodes)){
    node.states[[i]][row.names(node.states[[i]])==j,1]<-out[[i]]$ace[row.names(out[[i]]$ace)==j,1]
    node.states[[i]][row.names(node.states[[i]])==j,2]<-out[[i]]$ace[row.names(out[[i]]$ace)==j,2]
    node.states[[i]][row.names(node.states[[i]])==j,3]<-out[[i]]$ace[row.names(out[[i]]$ace)==j,3]
  }
}

anc.states<-data.frame(matrix(NA, nrow=length(nodes), ncol=3*length(colnames(PP.thresh))))
rownames(anc.states)<-nodes
colnames(anc.states)<-c("Ia", "Ia Low", "Ia High", "II", "II Low", "II High", "Ib", "Ib Low", "Ib High")

Ia<-NULL
II<-NULL
Ib<-NULL
for(i in unique(nodes)){
  for(j in 1:length(node.states)){
    Ia<-c(Ia, node.states[[j]][rownames(node.states[[j]])==i,1])
    II<-c(II, node.states[[j]][rownames(node.states[[j]])==i,2])
    Ib<-c(Ib, node.states[[j]][rownames(node.states[[j]])==i,3])
  }
  anc.states[rownames(anc.states)==i,"Ia"]<-median(Ia)
  anc.states[rownames(anc.states)==i,"Ia Low"]<-HPDinterval(as.mcmc(Ia))[1]
  anc.states[rownames(anc.states)==i,"Ia High"]<-HPDinterval(as.mcmc(Ia))[2]
  anc.states[rownames(anc.states)==i,"II"]<-median(II)
  anc.states[rownames(anc.states)==i,"II Low"]<-HPDinterval(as.mcmc(II))[1]
  anc.states[rownames(anc.states)==i,"II High"]<-HPDinterval(as.mcmc(II))[2]
  anc.states[rownames(anc.states)==i,"Ib"]<-median(Ib)
  anc.states[rownames(anc.states)==i,"Ib Low"]<-HPDinterval(as.mcmc(Ib))[1]
  anc.states[rownames(anc.states)==i,"Ib High"]<-HPDinterval(as.mcmc(Ib))[2]
  Ia<-NULL
  II<-NULL
  Ib<-NULL
}

# View median probability and 95% interval for state at ancestral node
anc.states[1,]

# Extract reconstructed states for tips, including those with TSD but unknown pattern
tips<-list()
for(i in 1:length(out)){
  tips[[i]]<-plot(out[[i]])[1:length(out[[i]]$tree$tip.label),]
}

tip.states<-data.frame(matrix(NA, nrow=length(out[[1]]$tree$tip.label), ncol=length(colnames(PP.thresh))))
rownames(tip.states)<-out[[1]]$tree$tip.label
colnames(tip.states)<-c("Ia","II","Ib")

tips.Ia<-NULL
tips.II<-NULL
tips.Ib<-NULL
for(i in unique(out[[1]]$tree$tip.label)){
  for(j in 1:length(tips)){
    tips.Ia<-c(tips.Ia, tips[[j]][rownames(tips[[j]])==i,1])
    tips.II<-c(tips.II, tips[[j]][rownames(tips[[j]])==i,2])
    tips.Ib<-c(tips.Ib, tips[[j]][rownames(tips[[j]])==i,3])
  }
  tip.states[rownames(tip.states)==i,"Ia"]<-median(tips.Ia)
  tip.states[rownames(tip.states)==i,"II"]<-median(tips.II)
  tip.states[rownames(tip.states)==i,"Ib"]<-median(tips.Ib)
  tips.Ia<-NULL
  tips.II<-NULL
  tips.Ib<-NULL
}

# Plot median ancestral state values on Thomson phylogeny
tip.states<-tip.states[match(thom$tip.label, rownames(tip.states)),]

png(filename = "ASE_Thomson_Crocs_Thresh.png", res = 300, width = 180, height = 170, units="mm")
par(mar = c(0,0,0,0), mfrow = c(1,1))

cols<-setNames(c("dodgerblue4","steelblue1","lightsteelblue1")[1:length(unique(colnames(PP.thresh)))],sort(unique(colnames(PP.thresh))))
h<-max(nodeHeights(thom))
offset.factor<-1.87
plotTree(geiger::rescale(thom,model="depth",depth=offset.factor*h),
         color="transparent",ftype="i",type="fan",fsize=0.55,cex=1,lwd=1,mar=c(0,0,0,0),
         maxY=118)
par(fg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
draw.circle(0,0,max(node.depth.edgelength(thom)),col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-50,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-100,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-150,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-200,col=NA,border="gray90")
draw.circle(0,0,max(node.depth.edgelength(thom))-250,col=NA,border="gray90")
plotTree(thom, type="fan", colors=cols, fsize=0.55,cex=1,ftype="i",lwd=1,mar=c(0,0,0,0), add=T, 
         xlim=obj$x.lim,ylim=obj$y.lim, maxY=118)
tiplabels(pie=tip.states, piecol=cols, cex=0.37)
for(i in unique(nodes)){
  if(i==118){
    nodelabels(pie=anc.states[rownames(anc.states)==i,c("Ia","II","Ib")], node=i, piecol=cols, cex=1.1)
  }else{
    nodelabels(pie=anc.states[rownames(anc.states)==i,c("Ia","II","Ib")], node=i, piecol=cols, cex=0.37) 
  }
}
par(fg="black")
obj<-axis(1,pos=3,at=seq(max(nodeHeights(thom)),max(nodeHeights(thom))-150,by=-50),cex.axis=1,labels=F,lwd=1.2,tcl=-0.2)
text(obj,rep(-10,length(obj)),max(nodeHeights(thom))-obj,cex=0.6)
add.simmap.legend(leg=c("Ia","II","Ib"), colors=cols, shape="circle", prompt=F, x=max(nodeHeights(thom))+140, y=15, fsize=0.75)

dev.off()
