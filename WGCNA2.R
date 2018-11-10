# network analysis ####
library(zoo)
library(stats)

# a=number averaged in rolling mean
# b=power
# x= value of x
# y= value of y

makeEQTN<-function(x,y,a,b){
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
#eqtn[2:6]<-paste(coefs[2:6], power, sep="")
fit<-function(x,y,eqtn){
  #variance explained/total variance
}


eqtn[1]
# choose model
model <- lm(y ~ poly(q,3))

# extract values 
# test case ####
a<-c(1:50, rnorm(10,25,2))
b<-c(rnorm(25, 30, 20),rnorm(25, 150, 30),rnorm(10, 100, 20))
eq<-makeEQTN(a,b,10,5)

# hierarchical clustering ####
library("ggplot2")
library("ggdendro")
ggdendrogram()
dend <- as.dendrogram(hc)
plot(as.phylo(hc), type = "unrooted", cex = 0.6,
     no.margin = TRUE)
dend %>% set("branches_k_color", k = 2) %>% 
  plot(main = "Default colors")