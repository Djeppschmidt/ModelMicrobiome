# validation/calibration run ####

# import species from /sppModel_independent.R
# species list= AllSpp
library(reshape2)

# integrate ...
integrate<-function(fun, Factors){
  #fun<-Randomspp.H
  require(quantmod)#test this is actually doing what I think it is...
  out<-matrix(data=NA, nrow=nrow(Factors), ncol = length(fun))
  for(i in 1:length(fun)) {
  for(row in 1:nrow(Factors)){
    out[row,i]<-do.call(get(fun[[i]]), list(Factors[row,1],Factors[row,2],Factors[row,3],Factors[row,4],Factors[row,5]))#this line get()
  }}
  
  Sites<-c("Site1","Site2","Site3","Site4","Site5","Site6","Site7","Site8","Site9","Site10","Site11","Site12","Site13","Site14","Site15","Site16","Site17","Site18","Site19","Site20","Site21","Site22","Site23","Site24","Site25","Site26","Site27","Site28","Site29","Site30")
  rownames(out)<-Sites
  colnames(out)<-fun
  out[out<0]<-0
  out<-round(out)
  out}


generateNewEnvrionment<-function(seed){
require(reshape2)
  #set.seed(seed)
f1c1<-c(5,5,5,5,5,5)
f1c2<-c(1,3,10,15,3,15)
f1c3<-c(0.5,0.5,3,3,1,5)
F1.frame<-mapply(rnorm, f1c1,f1c2,f1c3)
F1<-melt(F1.frame)

#F2
f2c1<-c(5,5,5,5,5,5)
f2c2<-c(34,30,50,55,35,60)
f2c3<-c(0.5,0.5,3,3,1,5)
F2.frame<-mapply(rnorm, f2c1,f2c2,f2c3)
F2<-melt(F2.frame)

#F3
f3c1<-c(5,5,5,5,5,5)
f3c2<-c(1,3,10,15,3,15)
f3c3<-c(0.5,0.5,3,3,1,5)
F3.frame<-mapply(rnorm, f3c1,f3c2,f3c3)
F3<-melt(F3.frame)

#F4
f4c1<-c(5,5,5,5,5,5)
f4c2<-c(1,3,10,15,3,15)
f4c3<-c(0.5,0.5,3,3,1,5)
F4.frame<-mapply(rnorm, f4c1,f4c2,f4c3)
F4<-melt(F4.frame)

#F5
f5c1<-c(5,5,5,5,5,5)
f5c2<-c(1,3,10,15,3,15)
f5c3<-c(0.5,0.5,3,3,1,5)
F5.frame<-mapply(rnorm, f5c1,f5c2,f5c3)
F5<-melt(F5.frame)

Factors<-data.frame(F1$value,F2$value,F3$value,F4$value,F5$value)
Sites<-c("Site1","Site2","Site3","Site4","Site5","Site6","Site7","Site8","Site9","Site10","Site11","Site12","Site13","Site14","Site15","Site16","Site17","Site18","Site19","Site20","Site21","Site22","Site23","Site24","Site25","Site26","Site27","Site28","Site29","Site30")
rownames(Factors)<-Sites
colnames(Factors)<-c("F1","F2","F3","F4","F5")
#head(Factors)

# output response table####
return(Factors)
}

# test ####
#test integrate
test<-integrate(Randomspp.H, Factors) #problem in spp function list...
test<-integrate(Randomspp.L, Factors)

test<-integrate(AllSpp, Factors) # still generating not rounded numbers...


# make key table ####
library(phyloseq)
otu<-otu_table(t(test), taxa_are_rows = TRUE)

# high group
sums.tax<-taxa_sums(otu)
high<-sums.tax>1000
med<-sums.tax<1000&sums.tax>100
low<-sums.tax<100
key<-data.frame(high,med,low)
