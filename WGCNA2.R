# network analysis ####

library(zoo)
library(stats)

makeEQTN<-funciton(x,y,a,b){
  require(zoo)
  require(stats)
  power<-c(1:b)
  power<-paste("*x^", power, sep="")
  
  model<-lm(rollmean(y, a)~poly(rollmean(x, a), b))
  reg<-summary(model)
  coefs<-reg$coefficients[1:(b+1)]
  eqtn<-coefs
  eqtn[2:b]<-paste(eqtn[2:b], power, sep="")
  eqtn<-paste(eqtn[1], eqtn[2], eqtn[3], eqtn[4], eqtn[5], eqtn[6], sep="+")
  noquote(eqtn)
}
eqtn[2:6]<-paste(coefs[2:6], power, sep="")
fit<-function(x,y,eqtn){}


eqtn[1]
# choose model
model <- lm(y ~ poly(q,3))

# extract values 


# hierarchical clustering ####
library("ggplot2")
library("ggdendro")
ggdendrogram()
dend <- as.dendrogram(hc)
plot(as.phylo(hc), type = "unrooted", cex = 0.6,
     no.margin = TRUE)
dend %>% set("branches_k_color", k = 2) %>% 
  plot(main = "Default colors")