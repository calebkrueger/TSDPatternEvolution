# Quantitative Genetics of Pattern II TSD

# Load the necessary packages
library("bootstrap")

# Read in sex-ratio data, including logit-transformed data
data<-read.csv("Chelydra_data.csv")

############################################################
# Calculate family mean genetic correlations

# Calculate family sex ratios (proportion male), counting intersexes as 1/2 male
data$sr21.5<-(data$Males21.5+(data$Intersex21.5/2))/data$n21.5
data$sr27.5<-(data$Males27.5+(data$Intersex27.5/2))/data$n27.5

# Calculate genetic correlation of untransformed sex ratios following Janzen 1992 for comparison
(cov(data$sr21.5,data$sr27.5))/((var(data$sr21.5)*var(data$sr27.5))^(1/2))

# Bootstrap to get 95% CI limits
theta<-function(x,data) {(cov(data[x,11],data[x,12]))/((var(data[x,11])*var(data[x,12]))^(1/2))}
r.jboot<-bootstrap(1:30,100000,theta,data)
quantile(r.jboot$thetastar, c(0.025,0.975))

# Calculate corrected family sex ratios for 0s and 1s
data$cor.sr21.5<-NA
data$cor.sr27.5<-NA
for(i in 1:nrow(data)){
  if(data[i,]$sr21.5==0){
    data[i,]$cor.sr21.5<-1/(2*data[i,]$n21.5)
  }else{
    if(data[i,]$sr21.5==1){
      data[i,]$cor.sr21.5<-1-(1/(2*data[i,]$n21.5))
    }else{
      data[i,]$cor.sr21.5<-data[i,]$sr21.5
    }
  }
}
for(i in 1:nrow(data)){
  if(data[i,]$sr27.5==0){
    data[i,]$cor.sr27.5<-1/(2*data[i,]$n27.5)
  }else{
    if(data[i,]$sr27.5==1){
      data[i,]$cor.sr27.5<-1-(1/(2*data[i,]$n27.5))
    }else{
      data[i,]$cor.sr27.5<-data[i,]$sr27.5
    }
  }
}

# Logit transform family sex ratios
data$log.sr21.5<-log(data$cor.sr21.5/(1-data$cor.sr21.5))
data$log.sr27.5<-log(data$cor.sr27.5/(1-data$cor.sr27.5))

# Calculate weight factors
data$wf21.5<-1/(data$n21.5*data$cor.sr21.5*(1-data$cor.sr21.5))
data$wf27.5<-1/(data$n27.5*data$cor.sr27.5*(1-data$cor.sr27.5))

# Calculate population sex ratios
data$pop.sr21.5<-(sum(data$Males21.5)+(sum(data$Intersex21.5)/2))/sum(data$n21.5)
data$pop.sr27.5<-(sum(data$Males27.5)+(sum(data$Intersex27.5)/2))/sum(data$n27.5)

# Logit transform population sex ratios
data$log.pop.sr21.5<-log(data$pop.sr21.5/(1-data$pop.sr21.5))
data$log.pop.sr27.5<-log(data$pop.sr27.5/(1-data$pop.sr27.5))

# Calculate family mean correlation
(sum((data$log.sr21.5-data$log.pop.sr21.5)*(data$log.sr27.5-data$log.pop.sr27.5)*(0.5*(data$wf21.5+data$wf27.5))))/(sum((0.5*(data$wf21.5+data$wf27.5))))/
  sqrt((sum(((data$log.sr21.5-data$log.pop.sr21.5)^2)*data$wf21.5))/(sum(data$wf21.5))*(sum(((data$log.sr27.5-data$log.pop.sr27.5)^2)*data$wf27.5))/(sum(data$wf27.5)))

# Bootstrap to get 95% CI limits
theta<-function(x, data) {(sum((data[x,15]-data[x,21])*(data[x,16]-data[x,22])*(0.5*(data[x,17]+data[x,18]))))/(sum((0.5*(data[x,17]+data[x,18]))))/
    sqrt((sum(((data[x,15]-data[x,21])^2)*data[x,17]))/(sum(data[x,17]))*(sum(((data[x,16]-data[x,22])^2)*data[x,18]))/(sum(data[x,18])))}
rboot<-bootstrap(1:30,100000,theta,data)
quantile(rboot$thetastar, c(0.025,0.975))

############################################################
# Calculate heritabilities of sex ratios

# Calculate pB for Appendix I in Bull et al. 1982 for data at 21.5
((sum((data[,11]-data[,19])^2))*(mean(data[,6]))/((nrow(data)-1)*(mean(data[,6])-1))-(((mean(data[,19])*(1-mean(data[,19]))))/(mean(data[,6])-1)))/
(mean(data[,19])*(1-mean(data[,19])))

# Bootstrap to get 95% CI for pB
theta<-function(x, data) {((sum((data[x,11]-data[x,19])^2))*(mean(data[x,6]))/((nrow(data)-1)*(mean(data[x,6])-1))-(((mean(data[x,19])*(1-mean(data[x,19]))))/(mean(data[x,6])-1)))/(mean(data[x,19])*(1-mean(data[x,19])))}
hboot_21.5<-bootstrap(1:30,100000,theta,data)
quantile(hboot_21.5$thetastar, c(0.025, 0.975))

# Calculate pB for Appendix I in Bull et al. 1982 for data at 27.5
((sum((data[,12]-data[,20])^2))*(mean(data[,10]))/((nrow(data)-1)*(mean(data[,10])-1))-(((mean(data[,20])*(1-mean(data[,20]))))/(mean(data[,10])-1)))/
  (mean(data[,20])*(1-mean(data[,20])))

# Bootstrap to get 95% CI for pB
theta<-function(x, data) {((sum((data[x,12]-data[x,20])^2))*(mean(data[x,10]))/((nrow(data)-1)*(mean(data[x,10])-1))-(((mean(data[x,20])*(1-mean(data[x,20]))))/(mean(data[x,10])-1)))/(mean(data[x,20])*(1-mean(data[x,20])))}
hboot_27.5<-bootstrap(1:30,100000,theta,data)
quantile(hboot_27.5$thetastar, c(0.025, 0.975))

# pB values are used to estimate px in Appendix I of Bull et al. 1982
# Heritability is then 2 * px