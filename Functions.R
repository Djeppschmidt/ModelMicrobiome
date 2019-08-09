# x = otu table with taxa as rows
# y = number of replicates desired for each treatment condition
# depth = simulated sequencing depth. Default = 250
# var = variance of sequencing totals. Default = 75
# adj = vector of probablility adjustments to simulate technical bias. Default = 1

sampling<-function(x, depth = 250, var = 75, adj = rep(1, ncol(x))){
  #data.frame(table(sample(rownames(x), round(rnorm(1, depth, var)), replace=T, adj*x/sum(x))))
  t1<-data.frame(unlist(table(sample(rownames(x), round(rnorm(1, depth, var)), replace=T, prob= adj*x[,1]/sum(x[,1])))))
  site.1<-t1$Freq
  names(site.1)<-as.character(t1$Var1)
  t1<-data.frame(unlist(table(sample(rownames(x), round(rnorm(1, depth, var)), replace=T, prob= adj*x[,2]/sum(x[,2])))))
  site.2<-t1$Freq
  names(site.2)<-as.character(t1$Var1)
  t1<-merge(as.data.frame(site.1), as.data.frame(site.2), by=0, all=T)
  rownames(t1)<-t1$Row.names
  t1<-subset(t1, select=c(names(t1)[names(t1)!="Row.names"]))
for (i in 3:ncol(x)){
  t2<-data.frame(unlist(table(sample(rownames(x), round(rnorm(1, depth, var)), replace=T, prob= adj*x[,(i)]/sum(x[,(i)])))))
  t.2<-t2$Freq
  names(t.2)<-t2$Var1
  t1<-merge(as.data.frame(t1), as.data.frame(t.2), by=0, all=T)
  rownames(t1)<-t1$Row.names
  t1<-subset(t1, select=c(names(t1)[names(t1)!="Row.names"]))
  names(t1)[names(t1)=="t.2"]<-paste("site", (i), sep=".")
    }
  t1[is.na(t1)]<-0
  t1
}

model.rarefy<-function(x, y, ...){
  out<-sampling(x, ...)
  colnames(out)<-paste(colnames(out), "1", sep=".")
  if (y>1) {
  for (i in 2:y){
    out2<-sampling(x, ...)
    colnames(out2)<-paste(colnames(out2), (i), sep=".")
    out<-merge(as.data.frame(out), as.data.frame(out2), by=0, all=T)
    rownames(out)<-out$Row.names
    out<-subset(out, select=c(names(out)[names(out)!="Row.names"]))
  }
    out[is.na(out)]<-0
  out
  } else {out}
}


# wrapper for making spp list: ####
make.spList<-function(n, replace=F){
  AllSpp<-c(paste0("spp", c(1:1724), sep=""))
  AllSpp<-lapply(AllSpp, get)
  AllSpp<-unlist(AllSpp)
  
  C<-base::sample(AllSpp, n, replace)
  C
}

make.comm<-function(Comm1, Factors){
otu<-matrix(data=NA, nrow=nrow(Factors), ncol = length(Comm1))
Sites<-c(paste0("Site", 1:30))
for(i in 1:length(Comm1)) {
  for(row in 1:nrow(Factors)){
   otu[row,i]<-do.call(Comm1[[i]], list(Factors[row,1],Factors[row,2],Factors[row,3],Factors[row,4],Factors[row,5]))
      }
}}

Build a script to run the following platforms:

DESeq2 #
model.deseq2<-function(){}
Limma Trend #
model.limma<-function(){}

EdgeR #essentially same as deseq2
# start with sample phyloseq object
model.edgeR<-function(PS){
Library(edgeR)
dgList <- DGEList(counts=Counts, genes=rownames(Counts))
dgList <- calcNormFactors(dgList, method="TMM") # tmm is not appropriate for microbial metabarcoding studies!!
designMat <- model.matrix(~Factor)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat) # other options are trended (as in limma-voom) and whole dataset...
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=4)
deGenes <- decideTestsDGE(lrt, p=0.001)
}
BBSeq # specifically designed for whole genome sequencing; inappropriate normalization methods

DSS #
model.DSS<-function(){
  
}
BaySeq #
model.BaySeq<-function(ps){
 require(baySeq) 
  require(phyloseq)
  tb<-as.data.matrix(as.matrix(otu_table(ps)))
  replicates<-sample_data(ps)$Factor
  groups<-list("NDE"=c(rep(1,ncol(tb))), "DE"=c(as.factor(sample_data(ps)$Factor)))
  CD<-new("countData", data=tb, replicates=replicates, groups=groups)
  libsizes(CD) <- getLibsizes(CD)
   CD <- getPriors.NB(CD, samplesize = 500, estimation = "QL", cl = NULL)# using a negative binomial distribution
   CD <- getLikelihoods(CD, cl = NULL, bootStraps = 3, verbose = FALSE)
   
}

ShrinkBayes #
model.ShrinkBays<-function(){
  require()
}
PoissonSeq #
model.PoissonSeq<-function(){
  require(PoissonSeq)
  a<-PS.Main()
  b<-
}


#define species functions
# 100 high response, linear functions Factor 1 ####
# 50 positive, 50 negative

#linear high responses factor A
spp1<-function(a,b,c,d,e) {1000+a*1000+(0*(a+b+c+d+e))}
spp2<-function(a,b,c,d,e) {5000+a*1000+(0*(a+b+c+d+e))}
spp3<-function(a,b,c,d,e) {10000+a*1000+(0*(a+b+c+d+e))}
spp4<-function(a,b,c,d,e) {15000+a*1000+(0*(a+b+c+d+e))}
spp5<-function(a,b,c,d,e) {20000+a*1000+(0*(a+b+c+d+e))}
spp6<-function(a,b,c,d,e) {25000+a*1000+(0*(a+b+c+d+e))}
spp7<-function(a,b,c,d,e) {30000+a*1000+(0*(a+b+c+d+e))}
spp8<-function(a,b,c,d,e) {35000+a*1000+(0*(a+b+c+d+e))}
spp9<-function(a,b,c,d,e) {40000+a*1000+(0*(a+b+c+d+e))}
spp10<-function(a,b,c,d,e) {50000+a*1000+(0*(a+b+c+d+e))}

spp11<-function(a,b,c,d,e) {1000-a*1000+(0*(a+b+c+d+e))}
spp12<-function(a,b,c,d,e) {5000-a*1000+(0*(a+b+c+d+e))}
spp13<-function(a,b,c,d,e) {10000-a*1000+(0*(a+b+c+d+e))}
spp14<-function(a,b,c,d,e) {15000-a*1000+(0*(a+b+c+d+e))}
spp15<-function(a,b,c,d,e) {20000-a*1000+(0*(a+b+c+d+e))}
spp16<-function(a,b,c,d,e) {25000-a*1000+(0*(a+b+c+d+e))}
spp17<-function(a,b,c,d,e) {30000-a*1000+(0*(a+b+c+d+e))}
spp18<-function(a,b,c,d,e) {35000-a*1000+(0*(a+b+c+d+e))}
spp19<-function(a,b,c,d,e) {40000-a*1000+(0*(a+b+c+d+e))}
spp20<-function(a,b,c,d,e) {50000-a*1000+(0*(a+b+c+d+e))}

spp21<-function(a,b,c,d,e) {10000+a*0.01+(0*(a+b+c+d+e))}
spp22<-function(a,b,c,d,e) {10000+a*0.1+(0*(a+b+c+d+e))}
spp23<-function(a,b,c,d,e) {10000+a+(0*(a+b+c+d+e))}
spp24<-function(a,b,c,d,e) {10000+a*10+(0*(a+b+c+d+e))}
spp25<-function(a,b,c,d,e) {10000+a*50+(0*(a+b+c+d+e))}
spp26<-function(a,b,c,d,e) {10000+a*100+(0*(a+b+c+d+e))}
spp27<-function(a,b,c,d,e) {10000+a*200+(0*(a+b+c+d+e))}
spp28<-function(a,b,c,d,e) {10000+a*500+(0*(a+b+c+d+e))}
spp29<-function(a,b,c,d,e) {10000+a*2000+(0*(a+b+c+d+e))}
spp30<-function(a,b,c,d,e) {10000+a*10000+(0*(a+b+c+d+e))}

spp31<-function(a,b,c,d,e) {10000-a*0.01+(0*(a+b+c+d+e))}
spp32<-function(a,b,c,d,e) {10000-a*0.1+(0*(a+b+c+d+e))}
spp33<-function(a,b,c,d,e) {10000-a+(0*(a+b+c+d+e))}
spp34<-function(a,b,c,d,e) {10000-a*10+(0*(a+b+c+d+e))}
spp35<-function(a,b,c,d,e) {10000-a*50+(0*(a+b+c+d+e))}
spp36<-function(a,b,c,d,e) {10000-a*100+(0*(a+b+c+d+e))}
spp37<-function(a,b,c,d,e) {10000-a*200+(0*(a+b+c+d+e))}
spp38<-function(a,b,c,d,e) {10000-a*500+(0*(a+b+c+d+e))}
spp39<-function(a,b,c,d,e) {10000-a*2000+(0*(a+b+c+d+e))}
spp40<-function(a,b,c,d,e) {10000-a*10000+(0*(a+b+c+d+e))}

spp41<-function(a,b,c,d,e) {1000+a*0.01+(0*(a+b+c+d+e))}
spp42<-function(a,b,c,d,e) {2000+a*0.1+(0*(a+b+c+d+e))}
spp43<-function(a,b,c,d,e) {5000+a+(0*(a+b+c+d+e))}
spp44<-function(a,b,c,d,e) {10000+a*10+(0*(a+b+c+d+e))}
spp45<-function(a,b,c,d,e) {15000+a*50+(0*(a+b+c+d+e))}
spp46<-function(a,b,c,d,e) {20000+a*100+(0*(a+b+c+d+e))}
spp47<-function(a,b,c,d,e) {25000+a*200+(0*(a+b+c+d+e))}
spp48<-function(a,b,c,d,e) {30000+a*500+(0*(a+b+c+d+e))}
spp49<-function(a,b,c,d,e) {40000+a*2000+(0*(a+b+c+d+e))}
spp50<-function(a,b,c,d,e) {50000+a*10000+(0*(a+b+c+d+e))}

spp51<-function(a,b,c,d,e) {50000-a*0.01+(0*(a+b+c+d+e))}
spp52<-function(a,b,c,d,e) {40000-a*0.1+(0*(a+b+c+d+e))}
spp53<-function(a,b,c,d,e) {30000-a+(0*(a+b+c+d+e))}
spp54<-function(a,b,c,d,e) {25000-a*10+(0*(a+b+c+d+e))}
spp55<-function(a,b,c,d,e) {20000-a*50+(0*(a+b+c+d+e))}
spp56<-function(a,b,c,d,e) {15000-a*100+(0*(a+b+c+d+e))}
spp57<-function(a,b,c,d,e) {10000-a*200+(0*(a+b+c+d+e))}
spp58<-function(a,b,c,d,e) {5000-a*500+(0*(a+b+c+d+e))}
spp59<-function(a,b,c,d,e) {2000-a*2000+(0*(a+b+c+d+e))}
spp60<-function(a,b,c,d,e) {1000-a*10000+(0*(a+b+c+d+e))}

spp61<-function(a,b,c,d,e) {1000-a*0.01+(0*(a+b+c+d+e))}
spp62<-function(a,b,c,d,e) {2000-a*0.1+(0*(a+b+c+d+e))}
spp63<-function(a,b,c,d,e) {5000-a+(0*(a+b+c+d+e))}
spp64<-function(a,b,c,d,e) {10000-a*10+(0*(a+b+c+d+e))}
spp65<-function(a,b,c,d,e) {15000-a*50+(0*(a+b+c+d+e))}
spp66<-function(a,b,c,d,e) {20000-a*100+(0*(a+b+c+d+e))}
spp67<-function(a,b,c,d,e) {25000-a*200+(0*(a+b+c+d+e))}
spp68<-function(a,b,c,d,e) {30000-a*500+(0*(a+b+c+d+e))}
spp69<-function(a,b,c,d,e) {40000-a*2000+(0*(a+b+c+d+e))}
spp70<-function(a,b,c,d,e) {50000-a*10000+(0*(a+b+c+d+e))}

spp71<-function(a,b,c,d,e) {50000+a*0.01+(0*(a+b+c+d+e))}
spp72<-function(a,b,c,d,e) {40000+a*0.1+(0*(a+b+c+d+e))}
spp73<-function(a,b,c,d,e) {30000+a+(0*(a+b+c+d+e))}
spp74<-function(a,b,c,d,e) {25000+a*10+(0*(a+b+c+d+e))}
spp75<-function(a,b,c,d,e) {20000+a*50+(0*(a+b+c+d+e))}
spp76<-function(a,b,c,d,e) {15000+a*100+(0*(a+b+c+d+e))}
spp77<-function(a,b,c,d,e) {10000+a*200+(0*(a+b+c+d+e))}
spp78<-function(a,b,c,d,e) {5000+a*500+(0*(a+b+c+d+e))}
spp79<-function(a,b,c,d,e) {2000+a*2000+(0*(a+b+c+d+e))}
spp80<-function(a,b,c,d,e) {1000+a*10000+(0*(a+b+c+d+e))}


##########################################
##                                      ##
##  linear low responses factor A       ##
##                                      ##
##########################################

spp81<-function(a,b,c,d,e) {10+a+(0*(a+b+c+d+e))}
spp82<-function(a,b,c,d,e) {10+a*2+(0*(a+b+c+d+e))}
spp83<-function(a,b,c,d,e) {10+a*3+(0*(a+b+c+d+e))}
spp84<-function(a,b,c,d,e) {10+a*4+(0*(a+b+c+d+e))}
spp85<-function(a,b,c,d,e) {10+a*5+(0*(a+b+c+d+e))}
spp86<-function(a,b,c,d,e) {10+a*6+(0*(a+b+c+d+e))}
spp87<-function(a,b,c,d,e) {10+a*7+(0*(a+b+c+d+e))}
spp88<-function(a,b,c,d,e) {10+a*8+(0*(a+b+c+d+e))}
spp89<-function(a,b,c,d,e) {10+a*9+(0*(a+b+c+d+e))}
spp90<-function(a,b,c,d,e) {10+a*10+(0*(a+b+c+d+e))}

spp91<-function(a,b,c,d,e) {10+a*11+(0*(a+b+c+d+e))}
spp92<-function(a,b,c,d,e) {10+a*20+(0*(a+b+c+d+e))}
spp93<-function(a,b,c,d,e) {10+a*30+(0*(a+b+c+d+e))}
spp94<-function(a,b,c,d,e) {10+a*40+(0*(a+b+c+d+e))}
spp95<-function(a,b,c,d,e) {10+a*50+(0*(a+b+c+d+e))}
spp96<-function(a,b,c,d,e) {10+a*100+(0*(a+b+c+d+e))}
spp97<-function(a,b,c,d,e) {10+a*110+(0*(a+b+c+d+e))}
spp98<-function(a,b,c,d,e) {10+a*120+(0*(a+b+c+d+e))}
spp99<-function(a,b,c,d,e) {10+a*140+(0*(a+b+c+d+e))}
spp100<-function(a,b,c,d,e) {10+a*150+(0*(a+b+c+d+e))}

spp101<-function(a,b,c,d,e) {10-a+(0*(a+b+c+d+e))}
spp102<-function(a,b,c,d,e) {10-a*2+(0*(a+b+c+d+e))}
spp103<-function(a,b,c,d,e) {10-a*3+(0*(a+b+c+d+e))}
spp104<-function(a,b,c,d,e) {10-a*4+(0*(a+b+c+d+e))}
spp105<-function(a,b,c,d,e) {10-a*5+(0*(a+b+c+d+e))}
spp106<-function(a,b,c,d,e) {10-a*6+(0*(a+b+c+d+e))}
spp107<-function(a,b,c,d,e) {10-a*7+(0*(a+b+c+d+e))}
spp108<-function(a,b,c,d,e) {10-a*8+(0*(a+b+c+d+e))}
spp109<-function(a,b,c,d,e) {10-a*9+(0*(a+b+c+d+e))}
spp110<-function(a,b,c,d,e) {10-a*10+(0*(a+b+c+d+e))}

spp111<-function(a,b,c,d,e) {10-a*11+(0*(a+b+c+d+e))}
spp112<-function(a,b,c,d,e) {10-a*20+(0*(a+b+c+d+e))}
spp113<-function(a,b,c,d,e) {10-a*30+(0*(a+b+c+d+e))}
spp114<-function(a,b,c,d,e) {10-a*40+(0*(a+b+c+d+e))}
spp115<-function(a,b,c,d,e) {10-a*50+(0*(a+b+c+d+e))}
spp116<-function(a,b,c,d,e) {10-a*100+(0*(a+b+c+d+e))}
spp117<-function(a,b,c,d,e) {10-a*110+(0*(a+b+c+d+e))}
spp118<-function(a,b,c,d,e) {10-a*120+(0*(a+b+c+d+e))}
spp119<-function(a,b,c,d,e) {10-a*140+(0*(a+b+c+d+e))}
spp120<-function(a,b,c,d,e) {10-a*150+(0*(a+b+c+d+e))}

spp121<-function(a,b,c,d,e) {50+a+(0*(a+b+c+d+e))}
spp122<-function(a,b,c,d,e) {50+a*2+(0*(a+b+c+d+e))}
spp123<-function(a,b,c,d,e) {50+a*3+(0*(a+b+c+d+e))}
spp124<-function(a,b,c,d,e) {50+a*4+(0*(a+b+c+d+e))}
spp125<-function(a,b,c,d,e) {50+a*5+(0*(a+b+c+d+e))}
spp126<-function(a,b,c,d,e) {50+a*6+(0*(a+b+c+d+e))}
spp127<-function(a,b,c,d,e) {50+a*7+(0*(a+b+c+d+e))}
spp128<-function(a,b,c,d,e) {50+a*8+(0*(a+b+c+d+e))}
spp129<-function(a,b,c,d,e) {50+a*9+(0*(a+b+c+d+e))}
spp130<-function(a,b,c,d,e) {50+a*10+(0*(a+b+c+d+e))}

spp131<-function(a,b,c,d,e) {50+a*11+(0*(a+b+c+d+e))}
spp132<-function(a,b,c,d,e) {50+a*20+(0*(a+b+c+d+e))}
spp133<-function(a,b,c,d,e) {50+a*30+(0*(a+b+c+d+e))}
spp134<-function(a,b,c,d,e) {50+a*40+(0*(a+b+c+d+e))}
spp135<-function(a,b,c,d,e) {50+a*50+(0*(a+b+c+d+e))}
spp136<-function(a,b,c,d,e) {50+a*100+(0*(a+b+c+d+e))}
spp137<-function(a,b,c,d,e) {50+a*110+(0*(a+b+c+d+e))}
spp138<-function(a,b,c,d,e) {50+a*120+(0*(a+b+c+d+e))}
spp139<-function(a,b,c,d,e) {50+a*140+(0*(a+b+c+d+e))}
spp140<-function(a,b,c,d,e) {50+a*150+(0*(a+b+c+d+e))}

spp141<-function(a,b,c,d,e) {50-a+(0*(a+b+c+d+e))}
spp142<-function(a,b,c,d,e) {50-a*2+(0*(a+b+c+d+e))}
spp143<-function(a,b,c,d,e) {50-a*3+(0*(a+b+c+d+e))}
spp144<-function(a,b,c,d,e) {50-a*4+(0*(a+b+c+d+e))}
spp145<-function(a,b,c,d,e) {50-a*5+(0*(a+b+c+d+e))}
spp146<-function(a,b,c,d,e) {50-a*6+(0*(a+b+c+d+e))}
spp147<-function(a,b,c,d,e) {50-a*7+(0*(a+b+c+d+e))}
spp148<-function(a,b,c,d,e) {50-a*8+(0*(a+b+c+d+e))}
spp149<-function(a,b,c,d,e) {50-a*9+(0*(a+b+c+d+e))}
spp150<-function(a,b,c,d,e) {50-a*10+(0*(a+b+c+d+e))}

spp151<-function(a,b,c,d,e) {50-a*11+(0*(a+b+c+d+e))}
spp152<-function(a,b,c,d,e) {50-a*20+(0*(a+b+c+d+e))}
spp153<-function(a,b,c,d,e) {50-a*30+(0*(a+b+c+d+e))}
spp154<-function(a,b,c,d,e) {50-a*40+(0*(a+b+c+d+e))}
spp155<-function(a,b,c,d,e) {50-a*50+(0*(a+b+c+d+e))}
spp156<-function(a,b,c,d,e) {50-a*100+(0*(a+b+c+d+e))}
spp157<-function(a,b,c,d,e) {50-a*110+(0*(a+b+c+d+e))}
spp158<-function(a,b,c,d,e) {50-a*120+(0*(a+b+c+d+e))}
spp159<-function(a,b,c,d,e) {50-a*140+(0*(a+b+c+d+e))}
spp160<-function(a,b,c,d,e) {50-a*150+(0*(a+b+c+d+e))}

##########################################
##                                      ##
##  linear medium responses factor A    ##
##                                      ##
##########################################

spp161<-function(a,b,c,d,e) {100+a+(0*(a+b+c+d+e))}
spp162<-function(a,b,c,d,e) {100+a*2+(0*(a+b+c+d+e))}
spp163<-function(a,b,c,d,e) {100+a*3+(0*(a+b+c+d+e))}
spp164<-function(a,b,c,d,e) {100+a*4+(0*(a+b+c+d+e))}
spp165<-function(a,b,c,d,e) {100+a*5+(0*(a+b+c+d+e))}
spp166<-function(a,b,c,d,e) {100+a*6+(0*(a+b+c+d+e))}
spp167<-function(a,b,c,d,e) {100+a*7+(0*(a+b+c+d+e))}
spp168<-function(a,b,c,d,e) {100+a*8+(0*(a+b+c+d+e))}
spp169<-function(a,b,c,d,e) {100+a*9+(0*(a+b+c+d+e))}
spp170<-function(a,b,c,d,e) {100+a*10+(0*(a+b+c+d+e))}

spp171<-function(a,b,c,d,e) {100+a*11+(0*(a+b+c+d+e))}
spp172<-function(a,b,c,d,e) {100+a*20+(0*(a+b+c+d+e))}
spp173<-function(a,b,c,d,e) {100+a*30+(0*(a+b+c+d+e))}
spp174<-function(a,b,c,d,e) {100+a*40+(0*(a+b+c+d+e))}
spp175<-function(a,b,c,d,e) {100+a*50+(0*(a+b+c+d+e))}
spp176<-function(a,b,c,d,e) {100+a*100+(0*(a+b+c+d+e))}
spp177<-function(a,b,c,d,e) {100+a*110+(0*(a+b+c+d+e))}
spp178<-function(a,b,c,d,e) {100+a*120+(0*(a+b+c+d+e))}
spp179<-function(a,b,c,d,e) {100+a*140+(0*(a+b+c+d+e))}
spp180<-function(a,b,c,d,e) {100+a*150+(0*(a+b+c+d+e))}

spp181<-function(a,b,c,d,e) {100-a+(0*(a+b+c+d+e))}
spp182<-function(a,b,c,d,e) {100-a*2+(0*(a+b+c+d+e))}
spp183<-function(a,b,c,d,e) {100-a*3+(0*(a+b+c+d+e))}
spp184<-function(a,b,c,d,e) {100-a*4+(0*(a+b+c+d+e))}
spp185<-function(a,b,c,d,e) {100-a*5+(0*(a+b+c+d+e))}
spp186<-function(a,b,c,d,e) {100-a*6+(0*(a+b+c+d+e))}
spp187<-function(a,b,c,d,e) {100-a*7+(0*(a+b+c+d+e))}
spp188<-function(a,b,c,d,e) {100-a*8+(0*(a+b+c+d+e))}
spp189<-function(a,b,c,d,e) {100-a*9+(0*(a+b+c+d+e))}
spp190<-function(a,b,c,d,e) {100-a*10+(0*(a+b+c+d+e))}

spp191<-function(a,b,c,d,e) {100-a*11+(0*(a+b+c+d+e))}
spp192<-function(a,b,c,d,e) {100-a*20+(0*(a+b+c+d+e))}
spp193<-function(a,b,c,d,e) {100-a*30+(0*(a+b+c+d+e))}
spp194<-function(a,b,c,d,e) {100-a*40+(0*(a+b+c+d+e))}
spp195<-function(a,b,c,d,e) {100-a*50+(0*(a+b+c+d+e))}
spp196<-function(a,b,c,d,e) {100-a*100+(0*(a+b+c+d+e))}
spp197<-function(a,b,c,d,e) {100-a*110+(0*(a+b+c+d+e))}
spp198<-function(a,b,c,d,e) {100-a*120+(0*(a+b+c+d+e))}
spp199<-function(a,b,c,d,e) {100-a*140+(0*(a+b+c+d+e))}
spp200<-function(a,b,c,d,e) {100-a*150+(0*(a+b+c+d+e))}

spp201<-function(a,b,c,d,e) {500+a+(0*(a+b+c+d+e))}
spp202<-function(a,b,c,d,e) {500+a*2+(0*(a+b+c+d+e))}
spp203<-function(a,b,c,d,e) {500+a*3+(0*(a+b+c+d+e))}
spp204<-function(a,b,c,d,e) {500+a*4+(0*(a+b+c+d+e))}
spp205<-function(a,b,c,d,e) {500+a*5+(0*(a+b+c+d+e))}
spp206<-function(a,b,c,d,e) {500+a*6+(0*(a+b+c+d+e))}
spp207<-function(a,b,c,d,e) {500+a*7+(0*(a+b+c+d+e))}
spp208<-function(a,b,c,d,e) {500+a*8+(0*(a+b+c+d+e))}
spp209<-function(a,b,c,d,e) {500+a*9+(0*(a+b+c+d+e))}
spp210<-function(a,b,c,d,e) {500+a*10+(0*(a+b+c+d+e))}

spp211<-function(a,b,c,d,e) {500+a*11+(0*(a+b+c+d+e))}
spp212<-function(a,b,c,d,e) {500+a*20+(0*(a+b+c+d+e))}
spp213<-function(a,b,c,d,e) {500+a*30+(0*(a+b+c+d+e))}
spp214<-function(a,b,c,d,e) {500+a*40+(0*(a+b+c+d+e))}
spp215<-function(a,b,c,d,e) {500+a*50+(0*(a+b+c+d+e))}
spp216<-function(a,b,c,d,e) {500+a*100+(0*(a+b+c+d+e))}
spp217<-function(a,b,c,d,e) {500+a*110+(0*(a+b+c+d+e))}
spp218<-function(a,b,c,d,e) {500+a*120+(0*(a+b+c+d+e))}
spp219<-function(a,b,c,d,e) {500+a*140+(0*(a+b+c+d+e))}
spp220<-function(a,b,c,d,e) {500+a*150+(0*(a+b+c+d+e))}

spp221<-function(a,b,c,d,e) {500-a+(0*(a+b+c+d+e))}
spp222<-function(a,b,c,d,e) {500-a*2+(0*(a+b+c+d+e))}
spp223<-function(a,b,c,d,e) {500-a*3+(0*(a+b+c+d+e))}
spp224<-function(a,b,c,d,e) {500-a*4+(0*(a+b+c+d+e))}
spp225<-function(a,b,c,d,e) {500-a*5+(0*(a+b+c+d+e))}
spp226<-function(a,b,c,d,e) {500-a*6+(0*(a+b+c+d+e))}
spp227<-function(a,b,c,d,e) {500-a*7+(0*(a+b+c+d+e))}
spp228<-function(a,b,c,d,e) {500-a*8+(0*(a+b+c+d+e))}
spp229<-function(a,b,c,d,e) {500-a*9+(0*(a+b+c+d+e))}
spp230<-function(a,b,c,d,e) {500-a*10+(0*(a+b+c+d+e))}

spp231<-function(a,b,c,d,e) {500-a*11+(0*(a+b+c+d+e))}
spp232<-function(a,b,c,d,e) {500-a*20+(0*(a+b+c+d+e))}
spp233<-function(a,b,c,d,e) {500-a*30+(0*(a+b+c+d+e))}
spp234<-function(a,b,c,d,e) {500-a*40+(0*(a+b+c+d+e))}
spp235<-function(a,b,c,d,e) {500-a*50+(0*(a+b+c+d+e))}
spp236<-function(a,b,c,d,e) {500-a*100+(0*(a+b+c+d+e))}
spp237<-function(a,b,c,d,e) {500-a*110+(0*(a+b+c+d+e))}
spp238<-function(a,b,c,d,e) {500-a*120+(0*(a+b+c+d+e))}
spp239<-function(a,b,c,d,e) {500-a*140+(0*(a+b+c+d+e))}
spp240<-function(a,b,c,d,e) {500-a*150+(0*(a+b+c+d+e))}

##########################################
##                                      ##
##  linear Hi responses factor B        ##
##                                      ##
##########################################

spp241<-function(a,b,c,d,e) {1000+b*1000+(0*(a+b+c+d+e))}
spp242<-function(a,b,c,d,e) {5000+b*1000+(0*(a+b+c+d+e))}
spp243<-function(a,b,c,d,e) {10000+b*1000+(0*(a+b+c+d+e))}
spp244<-function(a,b,c,d,e) {15000+b*1000+(0*(a+b+c+d+e))}
spp245<-function(a,b,c,d,e) {20000+b*1000+(0*(a+b+c+d+e))}
spp246<-function(a,b,c,d,e) {25000+b*1000+(0*(a+b+c+d+e))}
spp247<-function(a,b,c,d,e) {30000+b*1000+(0*(a+b+c+d+e))}
spp248<-function(a,b,c,d,e) {35000+b*1000+(0*(a+b+c+d+e))}
spp249<-function(a,b,c,d,e) {40000+b*1000+(0*(a+b+c+d+e))}
spp250<-function(a,b,c,d,e) {50000+b*1000+(0*(a+b+c+d+e))}

spp251<-function(a,b,c,d,e) {1000-b*1000+(0*(a+b+c+d+e))}
spp252<-function(a,b,c,d,e) {5000-b*1000+(0*(a+b+c+d+e))}
spp253<-function(a,b,c,d,e) {10000-b*1000+(0*(a+b+c+d+e))}
spp254<-function(a,b,c,d,e) {15000-b*1000+(0*(a+b+c+d+e))}
spp255<-function(a,b,c,d,e) {20000-b*1000+(0*(a+b+c+d+e))}
spp256<-function(a,b,c,d,e) {25000-b*1000+(0*(a+b+c+d+e))}
spp257<-function(a,b,c,d,e) {30000-b*1000+(0*(a+b+c+d+e))}
spp258<-function(a,b,c,d,e) {35000-b*1000+(0*(a+b+c+d+e))}
spp259<-function(a,b,c,d,e) {40000-b*1000+(0*(a+b+c+d+e))}
spp260<-function(a,b,c,d,e) {50000-b*1000+(0*(a+b+c+d+e))}

spp261<-function(a,b,c,d,e) {10000+b*0.01+(0*(a+b+c+d+e))}
spp262<-function(a,b,c,d,e) {10000+b*0.1+(0*(a+b+c+d+e))}
spp263<-function(a,b,c,d,e) {10000+b+(0*(a+b+c+d+e))}
spp264<-function(a,b,c,d,e) {10000+b*10+(0*(a+b+c+d+e))}
spp265<-function(a,b,c,d,e) {10000+b*50+(0*(a+b+c+d+e))}
spp266<-function(a,b,c,d,e) {10000+b*100+(0*(a+b+c+d+e))}
spp267<-function(a,b,c,d,e) {10000+b*200+(0*(a+b+c+d+e))}
spp268<-function(a,b,c,d,e) {10000+b*500+(0*(a+b+c+d+e))}
spp269<-function(a,b,c,d,e) {10000+b*2000+(0*(a+b+c+d+e))}
spp270<-function(a,b,c,d,e) {10000+b*10000+(0*(a+b+c+d+e))}

spp271<-function(a,b,c,d,e) {10000-b*0.01+(0*(a+b+c+d+e))}
spp272<-function(a,b,c,d,e) {10000-b*0.1+(0*(a+b+c+d+e))}
spp273<-function(a,b,c,d,e) {10000-b+(0*(a+b+c+d+e))}
spp274<-function(a,b,c,d,e) {10000-b*10+(0*(a+b+c+d+e))}
spp275<-function(a,b,c,d,e) {10000-b*50+(0*(a+b+c+d+e))}
spp276<-function(a,b,c,d,e) {10000-b*100+(0*(a+b+c+d+e))}
spp277<-function(a,b,c,d,e) {10000-b*200+(0*(a+b+c+d+e))}
spp278<-function(a,b,c,d,e) {10000-b*500+(0*(a+b+c+d+e))}
spp279<-function(a,b,c,d,e) {10000-b*2000+(0*(a+b+c+d+e))}
spp280<-function(a,b,c,d,e) {10000-b*10000+(0*(a+b+c+d+e))}

spp281<-function(a,b,c,d,e) {1000+b*0.01+(0*(a+b+c+d+e))}
spp282<-function(a,b,c,d,e) {2000+b*0.1+(0*(a+b+c+d+e))}
spp283<-function(a,b,c,d,e) {5000+b+(0*(a+b+c+d+e))}
spp284<-function(a,b,c,d,e) {10000+b*10+(0*(a+b+c+d+e))}
spp285<-function(a,b,c,d,e) {15000+b*50+(0*(a+b+c+d+e))}
spp286<-function(a,b,c,d,e) {20000+b*100+(0*(a+b+c+d+e))}
spp287<-function(a,b,c,d,e) {25000+b*200+(0*(a+b+c+d+e))}
spp288<-function(a,b,c,d,e) {30000+b*500+(0*(a+b+c+d+e))}
spp289<-function(a,b,c,d,e) {40000+b*2000+(0*(a+b+c+d+e))}
spp290<-function(a,b,c,d,e) {50000+b*10000+(0*(a+b+c+d+e))}

spp291<-function(a,b,c,d,e) {50000-b*0.01+(0*(a+b+c+d+e))}
spp292<-function(a,b,c,d,e) {40000-b*0.1+(0*(a+b+c+d+e))}
spp293<-function(a,b,c,d,e) {30000-b+(0*(a+b+c+d+e))}
spp294<-function(a,b,c,d,e) {25000-b*10+(0*(a+b+c+d+e))}
spp295<-function(a,b,c,d,e) {20000-b*50+(0*(a+b+c+d+e))}
spp296<-function(a,b,c,d,e) {15000-b*100+(0*(a+b+c+d+e))}
spp297<-function(a,b,c,d,e) {10000-b*200+(0*(a+b+c+d+e))}
spp298<-function(a,b,c,d,e) {5000-b*500+(0*(a+b+c+d+e))}
spp299<-function(a,b,c,d,e) {2000-b*2000+(0*(a+b+c+d+e))}
spp300<-function(a,b,c,d,e) {1000-b*10000+(0*(a+b+c+d+e))}

spp301<-function(a,b,c,d,e) {1000-b*0.01+(0*(a+b+c+d+e))}
spp302<-function(a,b,c,d,e) {2000-b*0.1+(0*(a+b+c+d+e))}
spp303<-function(a,b,c,d,e) {5000-b+(0*(a+b+c+d+e))}
spp304<-function(a,b,c,d,e) {10000-b*10+(0*(a+b+c+d+e))}
spp305<-function(a,b,c,d,e) {15000-b*50+(0*(a+b+c+d+e))}
spp306<-function(a,b,c,d,e) {20000-b*100+(0*(a+b+c+d+e))}
spp307<-function(a,b,c,d,e) {25000-b*200+(0*(a+b+c+d+e))}
spp308<-function(a,b,c,d,e) {30000-b*500+(0*(a+b+c+d+e))}
spp309<-function(a,b,c,d,e) {40000-b*2000+(0*(a+b+c+d+e))}
spp310<-function(a,b,c,d,e) {50000-b*10000+(0*(a+b+c+d+e))}

spp311<-function(a,b,c,d,e) {50000+b*0.01+(0*(a+b+c+d+e))}
spp312<-function(a,b,c,d,e) {40000+b*0.1+(0*(a+b+c+d+e))}
spp313<-function(a,b,c,d,e) {30000+b+(0*(a+b+c+d+e))}
spp314<-function(a,b,c,d,e) {25000+b*10+(0*(a+b+c+d+e))}
spp315<-function(a,b,c,d,e) {20000+b*50+(0*(a+b+c+d+e))}
spp316<-function(a,b,c,d,e) {15000+b*100+(0*(a+b+c+d+e))}
spp317<-function(a,b,c,d,e) {10000+b*200+(0*(a+b+c+d+e))}
spp318<-function(a,b,c,d,e) {5000+b*500+(0*(a+b+c+d+e))}
spp319<-function(a,b,c,d,e) {2000+b*2000+(0*(a+b+c+d+e))}
spp320<-function(a,b,c,d,e) {1000+b*10000+(0*(a+b+c+d+e))}


##########################################
##                                      ##
##  linear low responses factor B       ##
##                                      ##
##########################################

spp321<-function(a,b,c,d,e) {10+b+(0*(a+b+c+d+e))}
spp322<-function(a,b,c,d,e) {10+b*2+(0*(a+b+c+d+e))}
spp323<-function(a,b,c,d,e) {10+b*3+(0*(a+b+c+d+e))}
spp324<-function(a,b,c,d,e) {10+b*4+(0*(a+b+c+d+e))}
spp325<-function(a,b,c,d,e) {10+b*5+(0*(a+b+c+d+e))}
spp326<-function(a,b,c,d,e) {10+b*6+(0*(a+b+c+d+e))}
spp327<-function(a,b,c,d,e) {10+b*7+(0*(a+b+c+d+e))}
spp328<-function(a,b,c,d,e) {10+b*8+(0*(a+b+c+d+e))}
spp329<-function(a,b,c,d,e) {10+b*9+(0*(a+b+c+d+e))}
spp330<-function(a,b,c,d,e) {10+b*10+(0*(a+b+c+d+e))}

spp331<-function(a,b,c,d,e) {10+b*11+(0*(a+b+c+d+e))}
spp332<-function(a,b,c,d,e) {10+b*20+(0*(a+b+c+d+e))}
spp333<-function(a,b,c,d,e) {10+b*30+(0*(a+b+c+d+e))}
spp334<-function(a,b,c,d,e) {10+b*40+(0*(a+b+c+d+e))}
spp335<-function(a,b,c,d,e) {10+b*50+(0*(a+b+c+d+e))}
spp336<-function(a,b,c,d,e) {10+b*100+(0*(a+b+c+d+e))}
spp337<-function(a,b,c,d,e) {10+b*110+(0*(a+b+c+d+e))}
spp338<-function(a,b,c,d,e) {10+b*120+(0*(a+b+c+d+e))}
spp339<-function(a,b,c,d,e) {10+b*140+(0*(a+b+c+d+e))}
spp340<-function(a,b,c,d,e) {10+b*150+(0*(a+b+c+d+e))}

spp341<-function(a,b,c,d,e) {10-b+(0*(a+b+c+d+e))}
spp342<-function(a,b,c,d,e) {10-b*2+(0*(a+b+c+d+e))}
spp343<-function(a,b,c,d,e) {10-b*3+(0*(a+b+c+d+e))}
spp344<-function(a,b,c,d,e) {10-b*4+(0*(a+b+c+d+e))}
spp345<-function(a,b,c,d,e) {10-b*5+(0*(a+b+c+d+e))}
spp346<-function(a,b,c,d,e) {10-b*6+(0*(a+b+c+d+e))}
spp347<-function(a,b,c,d,e) {10-b*7+(0*(a+b+c+d+e))}
spp348<-function(a,b,c,d,e) {10-b*8+(0*(a+b+c+d+e))}
spp349<-function(a,b,c,d,e) {10-b*9+(0*(a+b+c+d+e))}
spp350<-function(a,b,c,d,e) {10-b*10+(0*(a+b+c+d+e))}

spp351<-function(a,b,c,d,e) {10-b*11+(0*(a+b+c+d+e))}
spp352<-function(a,b,c,d,e) {10-b*20+(0*(a+b+c+d+e))}
spp353<-function(a,b,c,d,e) {10-b*30+(0*(a+b+c+d+e))}
spp354<-function(a,b,c,d,e) {10-b*40+(0*(a+b+c+d+e))}
spp355<-function(a,b,c,d,e) {10-b*50+(0*(a+b+c+d+e))}
spp356<-function(a,b,c,d,e) {10-b*100+(0*(a+b+c+d+e))}
spp357<-function(a,b,c,d,e) {10-b*110+(0*(a+b+c+d+e))}
spp358<-function(a,b,c,d,e) {10-b*120+(0*(a+b+c+d+e))}
spp359<-function(a,b,c,d,e) {10-b*140+(0*(a+b+c+d+e))}
spp360<-function(a,b,c,d,e) {10-b*150+(0*(a+b+c+d+e))}

spp361<-function(a,b,c,d,e) {50+b+(0*(a+b+c+d+e))}
spp362<-function(a,b,c,d,e) {50+b*2+(0*(a+b+c+d+e))}
spp363<-function(a,b,c,d,e) {50+b*3+(0*(a+b+c+d+e))}
spp364<-function(a,b,c,d,e) {50+b*4+(0*(a+b+c+d+e))}
spp365<-function(a,b,c,d,e) {50+b*5+(0*(a+b+c+d+e))}
spp366<-function(a,b,c,d,e) {50+b*6+(0*(a+b+c+d+e))}
spp367<-function(a,b,c,d,e) {50+b*7+(0*(a+b+c+d+e))}
spp368<-function(a,b,c,d,e) {50+b*8+(0*(a+b+c+d+e))}
spp369<-function(a,b,c,d,e) {50+b*9+(0*(a+b+c+d+e))}
spp370<-function(a,b,c,d,e) {50+b*10+(0*(a+b+c+d+e))}

spp371<-function(a,b,c,d,e) {50+b*11+(0*(a+b+c+d+e))}
spp372<-function(a,b,c,d,e) {50+b*20+(0*(a+b+c+d+e))}
spp373<-function(a,b,c,d,e) {50+b*30+(0*(a+b+c+d+e))}
spp374<-function(a,b,c,d,e) {50+b*40+(0*(a+b+c+d+e))}
spp375<-function(a,b,c,d,e) {50+b*50+(0*(a+b+c+d+e))}
spp376<-function(a,b,c,d,e) {50+b*100+(0*(a+b+c+d+e))}
spp377<-function(a,b,c,d,e) {50+b*110+(0*(a+b+c+d+e))}
spp378<-function(a,b,c,d,e) {50+b*120+(0*(a+b+c+d+e))}
spp379<-function(a,b,c,d,e) {50+b*140+(0*(a+b+c+d+e))}
spp380<-function(a,b,c,d,e) {50+b*150+(0*(a+b+c+d+e))}

spp381<-function(a,b,c,d,e) {50-b+(0*(a+b+c+d+e))}
spp382<-function(a,b,c,d,e) {50-b*2+(0*(a+b+c+d+e))}
spp383<-function(a,b,c,d,e) {50-b*3+(0*(a+b+c+d+e))}
spp384<-function(a,b,c,d,e) {50-b*4+(0*(a+b+c+d+e))}
spp385<-function(a,b,c,d,e) {50-b*5+(0*(a+b+c+d+e))}
spp386<-function(a,b,c,d,e) {50-b*6+(0*(a+b+c+d+e))}
spp387<-function(a,b,c,d,e) {50-b*7+(0*(a+b+c+d+e))}
spp388<-function(a,b,c,d,e) {50-b*8+(0*(a+b+c+d+e))}
spp389<-function(a,b,c,d,e) {50-b*9+(0*(a+b+c+d+e))}
spp390<-function(a,b,c,d,e) {50-b*10+(0*(a+b+c+d+e))}

spp391<-function(a,b,c,d,e) {50-b*11+(0*(a+b+c+d+e))}
spp392<-function(a,b,c,d,e) {50-b*20+(0*(a+b+c+d+e))}
spp393<-function(a,b,c,d,e) {50-b*30+(0*(a+b+c+d+e))}
spp394<-function(a,b,c,d,e) {50-b*40+(0*(a+b+c+d+e))}
spp395<-function(a,b,c,d,e) {50-b*50+(0*(a+b+c+d+e))}
spp396<-function(a,b,c,d,e) {50-b*100+(0*(a+b+c+d+e))}
spp397<-function(a,b,c,d,e) {50-b*110+(0*(a+b+c+d+e))}
spp398<-function(a,b,c,d,e) {50-b*120+(0*(a+b+c+d+e))}
spp399<-function(a,b,c,d,e) {50-b*140+(0*(a+b+c+d+e))}
spp400<-function(a,b,c,d,e) {50-b*150+(0*(a+b+c+d+e))}

##########################################
##                                      ##
##  linear medium responses factor B    ##
##                                      ##
##########################################

spp401<-function(a,b,c,d,e) {100+b+(0*(a+b+c+d+e))}
spp402<-function(a,b,c,d,e) {100+b*2+(0*(a+b+c+d+e))}
spp403<-function(a,b,c,d,e) {100+b*3+(0*(a+b+c+d+e))}
spp404<-function(a,b,c,d,e) {100+b*4+(0*(a+b+c+d+e))}
spp405<-function(a,b,c,d,e) {100+b*5+(0*(a+b+c+d+e))}
spp406<-function(a,b,c,d,e) {100+b*6+(0*(a+b+c+d+e))}
spp407<-function(a,b,c,d,e) {100+b*7+(0*(a+b+c+d+e))}
spp408<-function(a,b,c,d,e) {100+b*8+(0*(a+b+c+d+e))}
spp409<-function(a,b,c,d,e) {100+b*9+(0*(a+b+c+d+e))}
spp410<-function(a,b,c,d,e) {100+b*10+(0*(a+b+c+d+e))}

spp411<-function(a,b,c,d,e) {100+b*11+(0*(a+b+c+d+e))}
spp412<-function(a,b,c,d,e) {100+b*20+(0*(a+b+c+d+e))}
spp413<-function(a,b,c,d,e) {100+b*30+(0*(a+b+c+d+e))}
spp414<-function(a,b,c,d,e) {100+b*40+(0*(a+b+c+d+e))}
spp415<-function(a,b,c,d,e) {100+b*50+(0*(a+b+c+d+e))}
spp416<-function(a,b,c,d,e) {100+b*100+(0*(a+b+c+d+e))}
spp417<-function(a,b,c,d,e) {100+b*110+(0*(a+b+c+d+e))}
spp418<-function(a,b,c,d,e) {100+b*120+(0*(a+b+c+d+e))}
spp419<-function(a,b,c,d,e) {100+b*140+(0*(a+b+c+d+e))}
spp420<-function(a,b,c,d,e) {100+b*150+(0*(a+b+c+d+e))}

spp421<-function(a,b,c,d,e) {100-b+(0*(a+b+c+d+e))}
spp422<-function(a,b,c,d,e) {100-b*2+(0*(a+b+c+d+e))}
spp423<-function(a,b,c,d,e) {100-b*3+(0*(a+b+c+d+e))}
spp424<-function(a,b,c,d,e) {100-b*4+(0*(a+b+c+d+e))}
spp425<-function(a,b,c,d,e) {100-b*5+(0*(a+b+c+d+e))}
spp426<-function(a,b,c,d,e) {100-b*6+(0*(a+b+c+d+e))}
spp427<-function(a,b,c,d,e) {100-b*7+(0*(a+b+c+d+e))}
spp428<-function(a,b,c,d,e) {100-b*8+(0*(a+b+c+d+e))}
spp429<-function(a,b,c,d,e) {100-b*9+(0*(a+b+c+d+e))}
spp430<-function(a,b,c,d,e) {100-b*10+(0*(a+b+c+d+e))}

spp431<-function(a,b,c,d,e) {100-b*11+(0*(a+b+c+d+e))}
spp432<-function(a,b,c,d,e) {100-b*20+(0*(a+b+c+d+e))}
spp433<-function(a,b,c,d,e) {100-b*30+(0*(a+b+c+d+e))}
spp434<-function(a,b,c,d,e) {100-b*40+(0*(a+b+c+d+e))}
spp435<-function(a,b,c,d,e) {100-b*50+(0*(a+b+c+d+e))}
spp436<-function(a,b,c,d,e) {100-b*100+(0*(a+b+c+d+e))}
spp437<-function(a,b,c,d,e) {100-b*110+(0*(a+b+c+d+e))}
spp438<-function(a,b,c,d,e) {100-b*120+(0*(a+b+c+d+e))}
spp439<-function(a,b,c,d,e) {100-b*140+(0*(a+b+c+d+e))}
spp440<-function(a,b,c,d,e) {100-b*150+(0*(a+b+c+d+e))}

spp441<-function(a,b,c,d,e) {500+b+(0*(a+b+c+d+e))}
spp442<-function(a,b,c,d,e) {500+b*2+(0*(a+b+c+d+e))}
spp443<-function(a,b,c,d,e) {500+b*3+(0*(a+b+c+d+e))}
spp444<-function(a,b,c,d,e) {500+b*4+(0*(a+b+c+d+e))}
spp445<-function(a,b,c,d,e) {500+b*5+(0*(a+b+c+d+e))}
spp446<-function(a,b,c,d,e) {500+b*6+(0*(a+b+c+d+e))}
spp447<-function(a,b,c,d,e) {500+b*7+(0*(a+b+c+d+e))}
spp448<-function(a,b,c,d,e) {500+b*8+(0*(a+b+c+d+e))}
spp449<-function(a,b,c,d,e) {500+b*9+(0*(a+b+c+d+e))}
spp450<-function(a,b,c,d,e) {500+b*10+(0*(a+b+c+d+e))}

spp451<-function(a,b,c,d,e) {500+b*11+(0*(a+b+c+d+e))}
spp452<-function(a,b,c,d,e) {500+b*20+(0*(a+b+c+d+e))}
spp453<-function(a,b,c,d,e) {500+b*30+(0*(a+b+c+d+e))}
spp454<-function(a,b,c,d,e) {500+b*40+(0*(a+b+c+d+e))}
spp455<-function(a,b,c,d,e) {500+b*50+(0*(a+b+c+d+e))}
spp456<-function(a,b,c,d,e) {500+b*100+(0*(a+b+c+d+e))}
spp457<-function(a,b,c,d,e) {500+b*110+(0*(a+b+c+d+e))}
spp458<-function(a,b,c,d,e) {500+b*120+(0*(a+b+c+d+e))}
spp459<-function(a,b,c,d,e) {500+b*140+(0*(a+b+c+d+e))}
spp460<-function(a,b,c,d,e) {500+b*150+(0*(a+b+c+d+e))}

spp461<-function(a,b,c,d,e) {500-b+(0*(a+b+c+d+e))}
spp462<-function(a,b,c,d,e) {500-b*2+(0*(a+b+c+d+e))}
spp463<-function(a,b,c,d,e) {500-b*3+(0*(a+b+c+d+e))}
spp464<-function(a,b,c,d,e) {500-b*4+(0*(a+b+c+d+e))}
spp465<-function(a,b,c,d,e) {500-b*5+(0*(a+b+c+d+e))}
spp466<-function(a,b,c,d,e) {500-b*6+(0*(a+b+c+d+e))}
spp467<-function(a,b,c,d,e) {500-b*7+(0*(a+b+c+d+e))}
spp468<-function(a,b,c,d,e) {500-b*8+(0*(a+b+c+d+e))}
spp469<-function(a,b,c,d,e) {500-b*9+(0*(a+b+c+d+e))}
spp470<-function(a,b,c,d,e) {500-b*10+(0*(a+b+c+d+e))}

spp471<-function(a,b,c,d,e) {500-b*11+(0*(a+b+c+d+e))}
spp472<-function(a,b,c,d,e) {500-b*20+(0*(a+b+c+d+e))}
spp473<-function(a,b,c,d,e) {500-b*30+(0*(a+b+c+d+e))}
spp474<-function(a,b,c,d,e) {500-b*40+(0*(a+b+c+d+e))}
spp475<-function(a,b,c,d,e) {500-b*50+(0*(a+b+c+d+e))}
spp476<-function(a,b,c,d,e) {500-b*100+(0*(a+b+c+d+e))}
spp477<-function(a,b,c,d,e) {500-b*110+(0*(a+b+c+d+e))}
spp478<-function(a,b,c,d,e) {500-b*120+(0*(a+b+c+d+e))}
spp479<-function(a,b,c,d,e) {500-b*140+(0*(a+b+c+d+e))}
spp480<-function(a,b,c,d,e) {500-b*150+(0*(a+b+c+d+e))}

##########################################
##                                      ##
##  linear Hi responses factor C        ##
##                                      ##
##########################################

spp481<-function(a,b,c,d,e) {1000+c*1000+(0*(a+b+c+d+e))}
spp482<-function(a,b,c,d,e) {5000+c*1000+(0*(a+b+c+d+e))}
spp483<-function(a,b,c,d,e) {10000+c*1000+(0*(a+b+c+d+e))}
spp484<-function(a,b,c,d,e) {15000+c*1000+(0*(a+b+c+d+e))}
spp485<-function(a,b,c,d,e) {20000+c*1000+(0*(a+b+c+d+e))}
spp486<-function(a,b,c,d,e) {25000+c*1000+(0*(a+b+c+d+e))}
spp487<-function(a,b,c,d,e) {30000+c*1000+(0*(a+b+c+d+e))}
spp488<-function(a,b,c,d,e) {35000+c*1000+(0*(a+b+c+d+e))}
spp489<-function(a,b,c,d,e) {40000+c*1000+(0*(a+b+c+d+e))}
spp490<-function(a,b,c,d,e) {50000+c*1000+(0*(a+b+c+d+e))}

spp491<-function(a,b,c,d,e) {1000-c*1000+(0*(a+b+c+d+e))}
spp492<-function(a,b,c,d,e) {5000-c*1000+(0*(a+b+c+d+e))}
spp493<-function(a,b,c,d,e) {10000-c*1000+(0*(a+b+c+d+e))}
spp494<-function(a,b,c,d,e) {15000-c*1000+(0*(a+b+c+d+e))}
spp495<-function(a,b,c,d,e) {20000-c*1000+(0*(a+b+c+d+e))}
spp496<-function(a,b,c,d,e) {25000-c*1000+(0*(a+b+c+d+e))}
spp497<-function(a,b,c,d,e) {30000-c*1000+(0*(a+b+c+d+e))}
spp498<-function(a,b,c,d,e) {35000-c*1000+(0*(a+b+c+d+e))}
spp499<-function(a,b,c,d,e) {40000-c*1000+(0*(a+b+c+d+e))}
spp500<-function(a,b,c,d,e) {50000-c*1000+(0*(a+b+c+d+e))}

spp501<-function(a,b,c,d,e) {10000+c*0.01+(0*(a+b+c+d+e))}
spp502<-function(a,b,c,d,e) {10000+c*0.1+(0*(a+b+c+d+e))}
spp503<-function(a,b,c,d,e) {10000+c+(0*(a+b+c+d+e))}
spp504<-function(a,b,c,d,e) {10000+c*10+(0*(a+b+c+d+e))}
spp505<-function(a,b,c,d,e) {10000+c*50+(0*(a+b+c+d+e))}
spp506<-function(a,b,c,d,e) {10000+c*100+(0*(a+b+c+d+e))}
spp507<-function(a,b,c,d,e) {10000+c*200+(0*(a+b+c+d+e))}
spp508<-function(a,b,c,d,e) {10000+c*500+(0*(a+b+c+d+e))}
spp509<-function(a,b,c,d,e) {10000+c*2000+(0*(a+b+c+d+e))}
spp510<-function(a,b,c,d,e) {10000+c*10000+(0*(a+b+c+d+e))}

spp511<-function(a,b,c,d,e) {10000-c*0.01+(0*(a+b+c+d+e))}
spp512<-function(a,b,c,d,e) {10000-c*0.1+(0*(a+b+c+d+e))}
spp513<-function(a,b,c,d,e) {10000-c+(0*(a+b+c+d+e))}
spp514<-function(a,b,c,d,e) {10000-c*10+(0*(a+b+c+d+e))}
spp515<-function(a,b,c,d,e) {10000-c*50+(0*(a+b+c+d+e))}
spp516<-function(a,b,c,d,e) {10000-c*100+(0*(a+b+c+d+e))}
spp517<-function(a,b,c,d,e) {10000-c*200+(0*(a+b+c+d+e))}
spp518<-function(a,b,c,d,e) {10000-c*500+(0*(a+b+c+d+e))}
spp519<-function(a,b,c,d,e) {10000-c*2000+(0*(a+b+c+d+e))}
spp520<-function(a,b,c,d,e) {10000-c*10000+(0*(a+b+c+d+e))}

spp521<-function(a,b,c,d,e) {1000+c*0.01+(0*(a+b+c+d+e))}
spp522<-function(a,b,c,d,e) {2000+c*0.1+(0*(a+b+c+d+e))}
spp523<-function(a,b,c,d,e) {5000+c+(0*(a+b+c+d+e))}
spp524<-function(a,b,c,d,e) {10000+c*10+(0*(a+b+c+d+e))}
spp525<-function(a,b,c,d,e) {15000+c*50+(0*(a+b+c+d+e))}
spp526<-function(a,b,c,d,e) {20000+c*100+(0*(a+b+c+d+e))}
spp527<-function(a,b,c,d,e) {25000+c*200+(0*(a+b+c+d+e))}
spp528<-function(a,b,c,d,e) {30000+c*500+(0*(a+b+c+d+e))}
spp529<-function(a,b,c,d,e) {40000+c*2000+(0*(a+b+c+d+e))}
spp530<-function(a,b,c,d,e) {50000+c*10000+(0*(a+b+c+d+e))}

spp531<-function(a,b,c,d,e) {50000-c*0.01+(0*(a+b+c+d+e))}
spp532<-function(a,b,c,d,e) {40000-c*0.1+(0*(a+b+c+d+e))}
spp533<-function(a,b,c,d,e) {30000-c+(0*(a+b+c+d+e))}
spp534<-function(a,b,c,d,e) {25000-c*10+(0*(a+b+c+d+e))}
spp535<-function(a,b,c,d,e) {20000-c*50+(0*(a+b+c+d+e))}
spp536<-function(a,b,c,d,e) {15000-c*100+(0*(a+b+c+d+e))}
spp537<-function(a,b,c,d,e) {10000-c*200+(0*(a+b+c+d+e))}
spp538<-function(a,b,c,d,e) {5000-c*500+(0*(a+b+c+d+e))}
spp539<-function(a,b,c,d,e) {2000-c*2000+(0*(a+b+c+d+e))}
spp540<-function(a,b,c,d,e) {1000-c*10000+(0*(a+b+c+d+e))}

spp541<-function(a,b,c,d,e) {1000-c*0.01+(0*(a+b+c+d+e))}
spp542<-function(a,b,c,d,e) {2000-c*0.1+(0*(a+b+c+d+e))}
spp543<-function(a,b,c,d,e) {5000-c+(0*(a+b+c+d+e))}
spp544<-function(a,b,c,d,e) {10000-c*10+(0*(a+b+c+d+e))}
spp545<-function(a,b,c,d,e) {15000-c*50+(0*(a+b+c+d+e))}
spp546<-function(a,b,c,d,e) {20000-c*100+(0*(a+b+c+d+e))}
spp547<-function(a,b,c,d,e) {25000-c*200+(0*(a+b+c+d+e))}
spp548<-function(a,b,c,d,e) {30000-c*500+(0*(a+b+c+d+e))}
spp549<-function(a,b,c,d,e) {40000-c*2000+(0*(a+b+c+d+e))}
spp550<-function(a,b,c,d,e) {50000-c*10000+(0*(a+b+c+d+e))}

spp551<-function(a,b,c,d,e) {50000+c*0.01+(0*(a+b+c+d+e))}
spp552<-function(a,b,c,d,e) {40000+c*0.1+(0*(a+b+c+d+e))}
spp553<-function(a,b,c,d,e) {30000+c+(0*(a+b+c+d+e))}
spp554<-function(a,b,c,d,e) {25000+c*10+(0*(a+b+c+d+e))}
spp555<-function(a,b,c,d,e) {20000+c*50+(0*(a+b+c+d+e))}
spp556<-function(a,b,c,d,e) {15000+c*100+(0*(a+b+c+d+e))}
spp557<-function(a,b,c,d,e) {10000+c*200+(0*(a+b+c+d+e))}
spp558<-function(a,b,c,d,e) {5000+c*500+(0*(a+b+c+d+e))}
spp559<-function(a,b,c,d,e) {2000+c*2000+(0*(a+b+c+d+e))}
spp560<-function(a,b,c,d,e) {1000+c*10000+(0*(a+b+c+d+e))}


##########################################
##                                      ##
##  linear low responses factor C       ##
##                                      ##
##########################################

spp561<-function(a,b,c,d,e) {10+c+(0*(a+b+c+d+e))}
spp562<-function(a,b,c,d,e) {10+c*2+(0*(a+b+c+d+e))}
spp563<-function(a,b,c,d,e) {10+c*3+(0*(a+b+c+d+e))}
spp564<-function(a,b,c,d,e) {10+c*4+(0*(a+b+c+d+e))}
spp565<-function(a,b,c,d,e) {10+c*5+(0*(a+b+c+d+e))}
spp566<-function(a,b,c,d,e) {10+c*6+(0*(a+b+c+d+e))}
spp567<-function(a,b,c,d,e) {10+c*7+(0*(a+b+c+d+e))}
spp568<-function(a,b,c,d,e) {10+c*8+(0*(a+b+c+d+e))}
spp569<-function(a,b,c,d,e) {10+c*9+(0*(a+b+c+d+e))}
spp570<-function(a,b,c,d,e) {10+c*10+(0*(a+b+c+d+e))}

spp571<-function(a,b,c,d,e) {10+c*11+(0*(a+b+c+d+e))}
spp572<-function(a,b,c,d,e) {10+c*20+(0*(a+b+c+d+e))}
spp573<-function(a,b,c,d,e) {10+c*30+(0*(a+b+c+d+e))}
spp574<-function(a,b,c,d,e) {10+c*40+(0*(a+b+c+d+e))}
spp575<-function(a,b,c,d,e) {10+c*50+(0*(a+b+c+d+e))}
spp576<-function(a,b,c,d,e) {10+c*100+(0*(a+b+c+d+e))}
spp577<-function(a,b,c,d,e) {10+c*110+(0*(a+b+c+d+e))}
spp578<-function(a,b,c,d,e) {10+c*120+(0*(a+b+c+d+e))}
spp579<-function(a,b,c,d,e) {10+c*140+(0*(a+b+c+d+e))}
spp580<-function(a,b,c,d,e) {10+c*150+(0*(a+b+c+d+e))}

spp581<-function(a,b,c,d,e) {10-c+(0*(a+b+c+d+e))}
spp582<-function(a,b,c,d,e) {10-c*2+(0*(a+b+c+d+e))}
spp583<-function(a,b,c,d,e) {10-c*3+(0*(a+b+c+d+e))}
spp584<-function(a,b,c,d,e) {10-c*4+(0*(a+b+c+d+e))}
spp585<-function(a,b,c,d,e) {10-c*5+(0*(a+b+c+d+e))}
spp586<-function(a,b,c,d,e) {10-c*6+(0*(a+b+c+d+e))}
spp587<-function(a,b,c,d,e) {10-c*7+(0*(a+b+c+d+e))}
spp588<-function(a,b,c,d,e) {10-c*8+(0*(a+b+c+d+e))}
spp589<-function(a,b,c,d,e) {10-c*9+(0*(a+b+c+d+e))}
spp590<-function(a,b,c,d,e) {10-c*10+(0*(a+b+c+d+e))}

spp591<-function(a,b,c,d,e) {10-c*11+(0*(a+b+c+d+e))}
spp592<-function(a,b,c,d,e) {10-c*20+(0*(a+b+c+d+e))}
spp593<-function(a,b,c,d,e) {10-c*30+(0*(a+b+c+d+e))}
spp594<-function(a,b,c,d,e) {10-c*40+(0*(a+b+c+d+e))}
spp595<-function(a,b,c,d,e) {10-c*50+(0*(a+b+c+d+e))}
spp596<-function(a,b,c,d,e) {10-c*100+(0*(a+b+c+d+e))}
spp597<-function(a,b,c,d,e) {10-c*110+(0*(a+b+c+d+e))}
spp598<-function(a,b,c,d,e) {10-c*120+(0*(a+b+c+d+e))}
spp599<-function(a,b,c,d,e) {10-c*140+(0*(a+b+c+d+e))}
spp600<-function(a,b,c,d,e) {10-c*150+(0*(a+b+c+d+e))}

spp601<-function(a,b,c,d,e) {50+c+(0*(a+b+c+d+e))}
spp602<-function(a,b,c,d,e) {50+c*2+(0*(a+b+c+d+e))}
spp603<-function(a,b,c,d,e) {50+c*3+(0*(a+b+c+d+e))}
spp604<-function(a,b,c,d,e) {50+c*4+(0*(a+b+c+d+e))}
spp605<-function(a,b,c,d,e) {50+c*5+(0*(a+b+c+d+e))}
spp606<-function(a,b,c,d,e) {50+c*6+(0*(a+b+c+d+e))}
spp607<-function(a,b,c,d,e) {50+c*7+(0*(a+b+c+d+e))}
spp608<-function(a,b,c,d,e) {50+c*8+(0*(a+b+c+d+e))}
spp609<-function(a,b,c,d,e) {50+c*9+(0*(a+b+c+d+e))}
spp610<-function(a,b,c,d,e) {50+c*10+(0*(a+b+c+d+e))}

spp611<-function(a,b,c,d,e) {50+c*11+(0*(a+b+c+d+e))}
spp612<-function(a,b,c,d,e) {50+c*20+(0*(a+b+c+d+e))}
spp613<-function(a,b,c,d,e) {50+c*30+(0*(a+b+c+d+e))}
spp614<-function(a,b,c,d,e) {50+c*40+(0*(a+b+c+d+e))}
spp615<-function(a,b,c,d,e) {50+c*50+(0*(a+b+c+d+e))}
spp616<-function(a,b,c,d,e) {50+c*100+(0*(a+b+c+d+e))}
spp617<-function(a,b,c,d,e) {50+c*110+(0*(a+b+c+d+e))}
spp618<-function(a,b,c,d,e) {50+c*120+(0*(a+b+c+d+e))}
spp619<-function(a,b,c,d,e) {50+c*140+(0*(a+b+c+d+e))}
spp620<-function(a,b,c,d,e) {50+c*150+(0*(a+b+c+d+e))}

spp621<-function(a,b,c,d,e) {50-c+(0*(a+b+c+d+e))}
spp622<-function(a,b,c,d,e) {50-c*2+(0*(a+b+c+d+e))}
spp623<-function(a,b,c,d,e) {50-c*3+(0*(a+b+c+d+e))}
spp624<-function(a,b,c,d,e) {50-c*4+(0*(a+b+c+d+e))}
spp625<-function(a,b,c,d,e) {50-c*5+(0*(a+b+c+d+e))}
spp626<-function(a,b,c,d,e) {50-c*6+(0*(a+b+c+d+e))}
spp627<-function(a,b,c,d,e) {50-c*7+(0*(a+b+c+d+e))}
spp628<-function(a,b,c,d,e) {50-c*8+(0*(a+b+c+d+e))}
spp629<-function(a,b,c,d,e) {50-c*9+(0*(a+b+c+d+e))}
spp630<-function(a,b,c,d,e) {50-c*10+(0*(a+b+c+d+e))}

spp631<-function(a,b,c,d,e) {50-c*11+(0*(a+b+c+d+e))}
spp632<-function(a,b,c,d,e) {50-c*20+(0*(a+b+c+d+e))}
spp633<-function(a,b,c,d,e) {50-c*30+(0*(a+b+c+d+e))}
spp634<-function(a,b,c,d,e) {50-c*40+(0*(a+b+c+d+e))}
spp635<-function(a,b,c,d,e) {50-c*50+(0*(a+b+c+d+e))}
spp636<-function(a,b,c,d,e) {50-c*100+(0*(a+b+c+d+e))}
spp637<-function(a,b,c,d,e) {50-c*110+(0*(a+b+c+d+e))}
spp638<-function(a,b,c,d,e) {50-c*120+(0*(a+b+c+d+e))}
spp639<-function(a,b,c,d,e) {50-c*140+(0*(a+b+c+d+e))}
spp640<-function(a,b,c,d,e) {50-c*150+(0*(a+b+c+d+e))}

##########################################
##                                      ##
##  linear medium responses factor C    ##
##                                      ##
##########################################

spp641<-function(a,b,c,d,e) {100+c+(0*(a+b+c+d+e))}
spp642<-function(a,b,c,d,e) {100+c*2+(0*(a+b+c+d+e))}
spp643<-function(a,b,c,d,e) {100+c*3+(0*(a+b+c+d+e))}
spp644<-function(a,b,c,d,e) {100+c*4+(0*(a+b+c+d+e))}
spp645<-function(a,b,c,d,e) {100+c*5+(0*(a+b+c+d+e))}
spp646<-function(a,b,c,d,e) {100+c*6+(0*(a+b+c+d+e))}
spp647<-function(a,b,c,d,e) {100+c*7+(0*(a+b+c+d+e))}
spp648<-function(a,b,c,d,e) {100+c*8+(0*(a+b+c+d+e))}
spp649<-function(a,b,c,d,e) {100+c*9+(0*(a+b+c+d+e))}
spp650<-function(a,b,c,d,e) {100+c*10+(0*(a+b+c+d+e))}

spp651<-function(a,b,c,d,e) {100+c*11+(0*(a+b+c+d+e))}
spp652<-function(a,b,c,d,e) {100+c*20+(0*(a+b+c+d+e))}
spp653<-function(a,b,c,d,e) {100+c*30+(0*(a+b+c+d+e))}
spp654<-function(a,b,c,d,e) {100+c*40+(0*(a+b+c+d+e))}
spp655<-function(a,b,c,d,e) {100+c*50+(0*(a+b+c+d+e))}
spp656<-function(a,b,c,d,e) {100+c*100+(0*(a+b+c+d+e))}
spp657<-function(a,b,c,d,e) {100+c*110+(0*(a+b+c+d+e))}
spp658<-function(a,b,c,d,e) {100+c*120+(0*(a+b+c+d+e))}
spp659<-function(a,b,c,d,e) {100+c*140+(0*(a+b+c+d+e))}
spp660<-function(a,b,c,d,e) {100+c*150+(0*(a+b+c+d+e))}

spp661<-function(a,b,c,d,e) {100-c+(0*(a+b+c+d+e))}
spp662<-function(a,b,c,d,e) {100-c*2+(0*(a+b+c+d+e))}
spp663<-function(a,b,c,d,e) {100-c*3+(0*(a+b+c+d+e))}
spp664<-function(a,b,c,d,e) {100-c*4+(0*(a+b+c+d+e))}
spp665<-function(a,b,c,d,e) {100-c*5+(0*(a+b+c+d+e))}
spp666<-function(a,b,c,d,e) {100-c*6+(0*(a+b+c+d+e))}
spp667<-function(a,b,c,d,e) {100-c*7+(0*(a+b+c+d+e))}
spp668<-function(a,b,c,d,e) {100-c*8+(0*(a+b+c+d+e))}
spp669<-function(a,b,c,d,e) {100-c*9+(0*(a+b+c+d+e))}
spp670<-function(a,b,c,d,e) {100-c*10+(0*(a+b+c+d+e))}

spp671<-function(a,b,c,d,e) {100-c*11+(0*(a+b+c+d+e))}
spp672<-function(a,b,c,d,e) {100-c*20+(0*(a+b+c+d+e))}
spp673<-function(a,b,c,d,e) {100-c*30+(0*(a+b+c+d+e))}
spp674<-function(a,b,c,d,e) {100-c*40+(0*(a+b+c+d+e))}
spp675<-function(a,b,c,d,e) {100-c*50+(0*(a+b+c+d+e))}
spp676<-function(a,b,c,d,e) {100-c*100+(0*(a+b+c+d+e))}
spp677<-function(a,b,c,d,e) {100-c*110+(0*(a+b+c+d+e))}
spp678<-function(a,b,c,d,e) {100-c*120+(0*(a+b+c+d+e))}
spp679<-function(a,b,c,d,e) {100-c*140+(0*(a+b+c+d+e))}
spp680<-function(a,b,c,d,e) {100-c*150+(0*(a+b+c+d+e))}

spp681<-function(a,b,c,d,e) {500+c+(0*(a+b+c+d+e))}
spp682<-function(a,b,c,d,e) {500+c*2+(0*(a+b+c+d+e))}
spp683<-function(a,b,c,d,e) {500+c*3+(0*(a+b+c+d+e))}
spp684<-function(a,b,c,d,e) {500+c*4+(0*(a+b+c+d+e))}
spp685<-function(a,b,c,d,e) {500+c*5+(0*(a+b+c+d+e))}
spp686<-function(a,b,c,d,e) {500+c*6+(0*(a+b+c+d+e))}
spp687<-function(a,b,c,d,e) {500+c*7+(0*(a+b+c+d+e))}
spp688<-function(a,b,c,d,e) {500+c*8+(0*(a+b+c+d+e))}
spp689<-function(a,b,c,d,e) {500+c*9+(0*(a+b+c+d+e))}
spp690<-function(a,b,c,d,e) {500+c*10+(0*(a+b+c+d+e))}

spp691<-function(a,b,c,d,e) {500+c*11+(0*(a+b+c+d+e))}
spp692<-function(a,b,c,d,e) {500+c*20+(0*(a+b+c+d+e))}
spp693<-function(a,b,c,d,e) {500+c*30+(0*(a+b+c+d+e))}
spp694<-function(a,b,c,d,e) {500+c*40+(0*(a+b+c+d+e))}
spp695<-function(a,b,c,d,e) {500+c*50+(0*(a+b+c+d+e))}
spp696<-function(a,b,c,d,e) {500+c*100+(0*(a+b+c+d+e))}
spp697<-function(a,b,c,d,e) {500+c*110+(0*(a+b+c+d+e))}
spp698<-function(a,b,c,d,e) {500+c*120+(0*(a+b+c+d+e))}
spp699<-function(a,b,c,d,e) {500+c*140+(0*(a+b+c+d+e))}
spp700<-function(a,b,c,d,e) {500+c*150+(0*(a+b+c+d+e))}

spp701<-function(a,b,c,d,e) {500-c+(0*(a+b+c+d+e))}
spp702<-function(a,b,c,d,e) {500-c*2+(0*(a+b+c+d+e))}
spp703<-function(a,b,c,d,e) {500-c*3+(0*(a+b+c+d+e))}
spp704<-function(a,b,c,d,e) {500-c*4+(0*(a+b+c+d+e))}
spp705<-function(a,b,c,d,e) {500-c*5+(0*(a+b+c+d+e))}
spp706<-function(a,b,c,d,e) {500-c*6+(0*(a+b+c+d+e))}
spp707<-function(a,b,c,d,e) {500-c*7+(0*(a+b+c+d+e))}
spp708<-function(a,b,c,d,e) {500-c*8+(0*(a+b+c+d+e))}
spp709<-function(a,b,c,d,e) {500-c*9+(0*(a+b+c+d+e))}
spp710<-function(a,b,c,d,e) {500-c*10+(0*(a+b+c+d+e))}

spp711<-function(a,b,c,d,e) {500-c*11+(0*(a+b+c+d+e))}
spp712<-function(a,b,c,d,e) {500-c*20+(0*(a+b+c+d+e))}
spp713<-function(a,b,c,d,e) {500-c*30+(0*(a+b+c+d+e))}
spp714<-function(a,b,c,d,e) {500-c*40+(0*(a+b+c+d+e))}
spp715<-function(a,b,c,d,e) {500-c*50+(0*(a+b+c+d+e))}
spp716<-function(a,b,c,d,e) {500-c*100+(0*(a+b+c+d+e))}
spp717<-function(a,b,c,d,e) {500-c*110+(0*(a+b+c+d+e))}
spp718<-function(a,b,c,d,e) {500-c*120+(0*(a+b+c+d+e))}
spp719<-function(a,b,c,d,e) {500-c*140+(0*(a+b+c+d+e))}
spp720<-function(a,b,c,d,e) {500-c*150+(0*(a+b+c+d+e))}

##########################################
##                                      ##
##  linear Hi responses factor D        ##
##                                      ##
##########################################

spp721<-function(a,b,c,d,e) {1000+d*1000+(0*(a+b+c+d+e))}
spp722<-function(a,b,c,d,e) {5000+d*1000+(0*(a+b+c+d+e))}
spp723<-function(a,b,c,d,e) {10000+d*1000+(0*(a+b+c+d+e))}
spp724<-function(a,b,c,d,e) {15000+d*1000+(0*(a+b+c+d+e))}
spp725<-function(a,b,c,d,e) {20000+d*1000+(0*(a+b+c+d+e))}
spp726<-function(a,b,c,d,e) {25000+d*1000+(0*(a+b+c+d+e))}
spp727<-function(a,b,c,d,e) {30000+d*1000+(0*(a+b+c+d+e))}
spp728<-function(a,b,c,d,e) {35000+d*1000+(0*(a+b+c+d+e))}
spp729<-function(a,b,c,d,e) {40000+d*1000+(0*(a+b+c+d+e))}
spp730<-function(a,b,c,d,e) {50000+d*1000+(0*(a+b+c+d+e))}

spp731<-function(a,b,c,d,e) {1000-d*1000+(0*(a+b+c+d+e))}
spp732<-function(a,b,c,d,e) {5000-d*1000+(0*(a+b+c+d+e))}
spp733<-function(a,b,c,d,e) {10000-d*1000+(0*(a+b+c+d+e))}
spp734<-function(a,b,c,d,e) {15000-d*1000+(0*(a+b+c+d+e))}
spp735<-function(a,b,c,d,e) {20000-d*1000+(0*(a+b+c+d+e))}
spp736<-function(a,b,c,d,e) {25000-d*1000+(0*(a+b+c+d+e))}
spp737<-function(a,b,c,d,e) {30000-d*1000+(0*(a+b+c+d+e))}
spp738<-function(a,b,c,d,e) {35000-d*1000+(0*(a+b+c+d+e))}
spp739<-function(a,b,c,d,e) {40000-d*1000+(0*(a+b+c+d+e))}
spp740<-function(a,b,c,d,e) {50000-d*1000+(0*(a+b+c+d+e))}

spp741<-function(a,b,c,d,e) {10000+d*0.01+(0*(a+b+c+d+e))}
spp742<-function(a,b,c,d,e) {10000+d*0.1+(0*(a+b+c+d+e))}
spp743<-function(a,b,c,d,e) {10000+d+(0*(a+b+c+d+e))}
spp744<-function(a,b,c,d,e) {10000+d*10+(0*(a+b+c+d+e))}
spp745<-function(a,b,c,d,e) {10000+d*50+(0*(a+b+c+d+e))}
spp746<-function(a,b,c,d,e) {10000+d*100+(0*(a+b+c+d+e))}
spp747<-function(a,b,c,d,e) {10000+d*200+(0*(a+b+c+d+e))}
spp748<-function(a,b,c,d,e) {10000+d*500+(0*(a+b+c+d+e))}
spp749<-function(a,b,c,d,e) {10000+d*2000+(0*(a+b+c+d+e))}
spp750<-function(a,b,c,d,e) {10000+d*10000+(0*(a+b+c+d+e))}

spp751<-function(a,b,c,d,e) {10000-d*0.01+(0*(a+b+c+d+e))}
spp752<-function(a,b,c,d,e) {10000-d*0.1+(0*(a+b+c+d+e))}
spp753<-function(a,b,c,d,e) {10000-d+(0*(a+b+c+d+e))}
spp754<-function(a,b,c,d,e) {10000-d*10+(0*(a+b+c+d+e))}
spp755<-function(a,b,c,d,e) {10000-d*50+(0*(a+b+c+d+e))}
spp756<-function(a,b,c,d,e) {10000-d*100+(0*(a+b+c+d+e))}
spp757<-function(a,b,c,d,e) {10000-d*200+(0*(a+b+c+d+e))}
spp758<-function(a,b,c,d,e) {10000-d*500+(0*(a+b+c+d+e))}
spp759<-function(a,b,c,d,e) {10000-d*2000+(0*(a+b+c+d+e))}
spp760<-function(a,b,c,d,e) {10000-d*10000+(0*(a+b+c+d+e))}

spp761<-function(a,b,c,d,e) {1000+d*0.01+(0*(a+b+c+d+e))}
spp762<-function(a,b,c,d,e) {2000+d*0.1+(0*(a+b+c+d+e))}
spp763<-function(a,b,c,d,e) {5000+d+(0*(a+b+c+d+e))}
spp764<-function(a,b,c,d,e) {10000+d*10+(0*(a+b+c+d+e))}
spp765<-function(a,b,c,d,e) {15000+d*50+(0*(a+b+c+d+e))}
spp766<-function(a,b,c,d,e) {20000+d*100+(0*(a+b+c+d+e))}
spp767<-function(a,b,c,d,e) {25000+d*200+(0*(a+b+c+d+e))}
spp768<-function(a,b,c,d,e) {30000+d*500+(0*(a+b+c+d+e))}
spp769<-function(a,b,c,d,e) {40000+d*2000+(0*(a+b+c+d+e))}
spp770<-function(a,b,c,d,e) {50000+d*10000+(0*(a+b+c+d+e))}

spp771<-function(a,b,c,d,e) {50000-d*0.01+(0*(a+b+c+d+e))}
spp772<-function(a,b,c,d,e) {40000-d*0.1+(0*(a+b+c+d+e))}
spp773<-function(a,b,c,d,e) {30000-d+(0*(a+b+c+d+e))}
spp774<-function(a,b,c,d,e) {25000-d*10+(0*(a+b+c+d+e))}
spp775<-function(a,b,c,d,e) {20000-d*50+(0*(a+b+c+d+e))}
spp776<-function(a,b,c,d,e) {15000-d*100+(0*(a+b+c+d+e))}
spp777<-function(a,b,c,d,e) {10000-d*200+(0*(a+b+c+d+e))}
spp778<-function(a,b,c,d,e) {5000-d*500+(0*(a+b+c+d+e))}
spp779<-function(a,b,c,d,e) {2000-d*2000+(0*(a+b+c+d+e))}
spp780<-function(a,b,c,d,e) {1000-d*10000+(0*(a+b+c+d+e))}

spp781<-function(a,b,c,d,e) {1000-d*0.01+(0*(a+b+c+d+e))}
spp782<-function(a,b,c,d,e) {2000-d*0.1+(0*(a+b+c+d+e))}
spp783<-function(a,b,c,d,e) {5000-d+(0*(a+b+c+d+e))}
spp784<-function(a,b,c,d,e) {10000-d*10+(0*(a+b+c+d+e))}
spp785<-function(a,b,c,d,e) {15000-d*50+(0*(a+b+c+d+e))}
spp786<-function(a,b,c,d,e) {20000-d*100+(0*(a+b+c+d+e))}
spp787<-function(a,b,c,d,e) {25000-d*200+(0*(a+b+c+d+e))}
spp788<-function(a,b,c,d,e) {30000-d*500+(0*(a+b+c+d+e))}
spp789<-function(a,b,c,d,e) {40000-d*2000+(0*(a+b+c+d+e))}
spp790<-function(a,b,c,d,e) {50000-d*10000+(0*(a+b+c+d+e))}

spp791<-function(a,b,c,d,e) {50000+d*0.01+(0*(a+b+c+d+e))}
spp792<-function(a,b,c,d,e) {40000+d*0.1+(0*(a+b+c+d+e))}
spp793<-function(a,b,c,d,e) {30000+d+(0*(a+b+c+d+e))}
spp794<-function(a,b,c,d,e) {25000+d*10+(0*(a+b+c+d+e))}
spp795<-function(a,b,c,d,e) {20000+d*50+(0*(a+b+c+d+e))}
spp796<-function(a,b,c,d,e) {15000+d*100+(0*(a+b+c+d+e))}
spp797<-function(a,b,c,d,e) {10000+d*200+(0*(a+b+c+d+e))}
spp798<-function(a,b,c,d,e) {5000+d*500+(0*(a+b+c+d+e))}
spp799<-function(a,b,c,d,e) {2000+d*2000+(0*(a+b+c+d+e))}
spp800<-function(a,b,c,d,e) {1000+d*10000+(0*(a+b+c+d+e))}


##########################################
##                                      ##
##  linear low responses factor D       ##
##                                      ##
##########################################

spp801<-function(a,b,c,d,e) {10+d+(0*(a+b+c+d+e))}
spp802<-function(a,b,c,d,e) {10+d*2+(0*(a+b+c+d+e))}
spp803<-function(a,b,c,d,e) {10+d*3+(0*(a+b+c+d+e))}
spp804<-function(a,b,c,d,e) {10+d*4+(0*(a+b+c+d+e))}
spp805<-function(a,b,c,d,e) {10+d*5+(0*(a+b+c+d+e))}
spp806<-function(a,b,c,d,e) {10+d*6+(0*(a+b+c+d+e))}
spp807<-function(a,b,c,d,e) {10+d*7+(0*(a+b+c+d+e))}
spp808<-function(a,b,c,d,e) {10+d*8+(0*(a+b+c+d+e))}
spp809<-function(a,b,c,d,e) {10+d*9+(0*(a+b+c+d+e))}
spp810<-function(a,b,c,d,e) {10+d*10+(0*(a+b+c+d+e))}

spp811<-function(a,b,c,d,e) {10+d*11+(0*(a+b+c+d+e))}
spp812<-function(a,b,c,d,e) {10+d*20+(0*(a+b+c+d+e))}
spp813<-function(a,b,c,d,e) {10+d*30+(0*(a+b+c+d+e))}
spp814<-function(a,b,c,d,e) {10+d*40+(0*(a+b+c+d+e))}
spp815<-function(a,b,c,d,e) {10+d*50+(0*(a+b+c+d+e))}
spp816<-function(a,b,c,d,e) {10+d*100+(0*(a+b+c+d+e))}
spp817<-function(a,b,c,d,e) {10+d*110+(0*(a+b+c+d+e))}
spp818<-function(a,b,c,d,e) {10+d*120+(0*(a+b+c+d+e))}
spp819<-function(a,b,c,d,e) {10+d*140+(0*(a+b+c+d+e))}
spp820<-function(a,b,c,d,e) {10+d*150+(0*(a+b+c+d+e))}

spp821<-function(a,b,c,d,e) {10-d+(0*(a+b+c+d+e))}
spp822<-function(a,b,c,d,e) {10-d*2+(0*(a+b+c+d+e))}
spp823<-function(a,b,c,d,e) {10-d*3+(0*(a+b+c+d+e))}
spp824<-function(a,b,c,d,e) {10-d*4+(0*(a+b+c+d+e))}
spp825<-function(a,b,c,d,e) {10-d*5+(0*(a+b+c+d+e))}
spp826<-function(a,b,c,d,e) {10-d*6+(0*(a+b+c+d+e))}
spp827<-function(a,b,c,d,e) {10-d*7+(0*(a+b+c+d+e))}
spp828<-function(a,b,c,d,e) {10-d*8+(0*(a+b+c+d+e))}
spp829<-function(a,b,c,d,e) {10-d*9+(0*(a+b+c+d+e))}
spp830<-function(a,b,c,d,e) {10-d*10+(0*(a+b+c+d+e))}

spp831<-function(a,b,c,d,e) {10-d*11+(0*(a+b+c+d+e))}
spp832<-function(a,b,c,d,e) {10-d*20+(0*(a+b+c+d+e))}
spp833<-function(a,b,c,d,e) {10-d*30+(0*(a+b+c+d+e))}
spp834<-function(a,b,c,d,e) {10-d*40+(0*(a+b+c+d+e))}
spp835<-function(a,b,c,d,e) {10-d*50+(0*(a+b+c+d+e))}
spp836<-function(a,b,c,d,e) {10-d*100+(0*(a+b+c+d+e))}
spp837<-function(a,b,c,d,e) {10-d*110+(0*(a+b+c+d+e))}
spp838<-function(a,b,c,d,e) {10-d*120+(0*(a+b+c+d+e))}
spp839<-function(a,b,c,d,e) {10-d*140+(0*(a+b+c+d+e))}
spp840<-function(a,b,c,d,e) {10-d*150+(0*(a+b+c+d+e))}

spp841<-function(a,b,c,d,e) {50+d+(0*(a+b+c+d+e))}
spp842<-function(a,b,c,d,e) {50+d*2+(0*(a+b+c+d+e))}
spp843<-function(a,b,c,d,e) {50+d*3+(0*(a+b+c+d+e))}
spp844<-function(a,b,c,d,e) {50+d*4+(0*(a+b+c+d+e))}
spp845<-function(a,b,c,d,e) {50+d*5+(0*(a+b+c+d+e))}
spp846<-function(a,b,c,d,e) {50+d*6+(0*(a+b+c+d+e))}
spp847<-function(a,b,c,d,e) {50+d*7+(0*(a+b+c+d+e))}
spp848<-function(a,b,c,d,e) {50+d*8+(0*(a+b+c+d+e))}
spp849<-function(a,b,c,d,e) {50+d*9+(0*(a+b+c+d+e))}
spp850<-function(a,b,c,d,e) {50+d*10+(0*(a+b+c+d+e))}

spp851<-function(a,b,c,d,e) {50+d*11+(0*(a+b+c+d+e))}
spp852<-function(a,b,c,d,e) {50+d*20+(0*(a+b+c+d+e))}
spp853<-function(a,b,c,d,e) {50+d*30+(0*(a+b+c+d+e))}
spp854<-function(a,b,c,d,e) {50+d*40+(0*(a+b+c+d+e))}
spp855<-function(a,b,c,d,e) {50+d*50+(0*(a+b+c+d+e))}
spp856<-function(a,b,c,d,e) {50+d*100+(0*(a+b+c+d+e))}
spp857<-function(a,b,c,d,e) {50+d*110+(0*(a+b+c+d+e))}
spp858<-function(a,b,c,d,e) {50+d*120+(0*(a+b+c+d+e))}
spp859<-function(a,b,c,d,e) {50+d*140+(0*(a+b+c+d+e))}
spp860<-function(a,b,c,d,e) {50+d*150+(0*(a+b+c+d+e))}

spp861<-function(a,b,c,d,e) {50-d+(0*(a+b+c+d+e))}
spp862<-function(a,b,c,d,e) {50-d*2+(0*(a+b+c+d+e))}
spp863<-function(a,b,c,d,e) {50-d*3+(0*(a+b+c+d+e))}
spp864<-function(a,b,c,d,e) {50-d*4+(0*(a+b+c+d+e))}
spp865<-function(a,b,c,d,e) {50-d*5+(0*(a+b+c+d+e))}
spp866<-function(a,b,c,d,e) {50-d*6+(0*(a+b+c+d+e))}
spp867<-function(a,b,c,d,e) {50-d*7+(0*(a+b+c+d+e))}
spp868<-function(a,b,c,d,e) {50-d*8+(0*(a+b+c+d+e))}
spp869<-function(a,b,c,d,e) {50-d*9+(0*(a+b+c+d+e))}
spp870<-function(a,b,c,d,e) {50-d*10+(0*(a+b+c+d+e))}

spp871<-function(a,b,c,d,e) {50-d*11+(0*(a+b+c+d+e))}
spp872<-function(a,b,c,d,e) {50-d*20+(0*(a+b+c+d+e))}
spp873<-function(a,b,c,d,e) {50-d*30+(0*(a+b+c+d+e))}
spp874<-function(a,b,c,d,e) {50-d*40+(0*(a+b+c+d+e))}
spp875<-function(a,b,c,d,e) {50-d*50+(0*(a+b+c+d+e))}
spp876<-function(a,b,c,d,e) {50-d*100+(0*(a+b+c+d+e))}
spp877<-function(a,b,c,d,e) {50-d*110+(0*(a+b+c+d+e))}
spp878<-function(a,b,c,d,e) {50-d*120+(0*(a+b+c+d+e))}
spp879<-function(a,b,c,d,e) {50-d*140+(0*(a+b+c+d+e))}
spp880<-function(a,b,c,d,e) {50-d*150+(0*(a+b+c+d+e))}

##########################################
##                                      ##
##  linear medium responses factor D    ##
##                                      ##
##########################################

spp881<-function(a,b,c,d,e) {100+d+(0*(a+b+c+d+e))}
spp882<-function(a,b,c,d,e) {100+d*2+(0*(a+b+c+d+e))}
spp883<-function(a,b,c,d,e) {100+d*3+(0*(a+b+c+d+e))}
spp884<-function(a,b,c,d,e) {100+d*4+(0*(a+b+c+d+e))}
spp885<-function(a,b,c,d,e) {100+d*5+(0*(a+b+c+d+e))}
spp886<-function(a,b,c,d,e) {100+d*6+(0*(a+b+c+d+e))}
spp887<-function(a,b,c,d,e) {100+d*7+(0*(a+b+c+d+e))}
spp888<-function(a,b,c,d,e) {100+d*8+(0*(a+b+c+d+e))}
spp889<-function(a,b,c,d,e) {100+d*9+(0*(a+b+c+d+e))}
spp890<-function(a,b,c,d,e) {100+d*10+(0*(a+b+c+d+e))}

spp891<-function(a,b,c,d,e) {100+d*11+(0*(a+b+c+d+e))}
spp892<-function(a,b,c,d,e) {100+d*20+(0*(a+b+c+d+e))}
spp893<-function(a,b,c,d,e) {100+d*30+(0*(a+b+c+d+e))}
spp894<-function(a,b,c,d,e) {100+d*40+(0*(a+b+c+d+e))}
spp895<-function(a,b,c,d,e) {100+d*50+(0*(a+b+c+d+e))}
spp896<-function(a,b,c,d,e) {100+d*100+(0*(a+b+c+d+e))}
spp897<-function(a,b,c,d,e) {100+d*110+(0*(a+b+c+d+e))}
spp898<-function(a,b,c,d,e) {100+d*120+(0*(a+b+c+d+e))}
spp899<-function(a,b,c,d,e) {100+d*140+(0*(a+b+c+d+e))}
spp900<-function(a,b,c,d,e) {100+d*150+(0*(a+b+c+d+e))}

spp901<-function(a,b,c,d,e) {100-d+(0*(a+b+c+d+e))}
spp902<-function(a,b,c,d,e) {100-d*2+(0*(a+b+c+d+e))}
spp903<-function(a,b,c,d,e) {100-d*3+(0*(a+b+c+d+e))}
spp904<-function(a,b,c,d,e) {100-d*4+(0*(a+b+c+d+e))}
spp905<-function(a,b,c,d,e) {100-d*5+(0*(a+b+c+d+e))}
spp906<-function(a,b,c,d,e) {100-d*6+(0*(a+b+c+d+e))}
spp907<-function(a,b,c,d,e) {100-d*7+(0*(a+b+c+d+e))}
spp908<-function(a,b,c,d,e) {100-d*8+(0*(a+b+c+d+e))}
spp909<-function(a,b,c,d,e) {100-d*9+(0*(a+b+c+d+e))}
spp910<-function(a,b,c,d,e) {100-d*10+(0*(a+b+c+d+e))}

spp911<-function(a,b,c,d,e) {100-d*11+(0*(a+b+c+d+e))}
spp912<-function(a,b,c,d,e) {100-d*20+(0*(a+b+c+d+e))}
spp913<-function(a,b,c,d,e) {100-d*30+(0*(a+b+c+d+e))}
spp914<-function(a,b,c,d,e) {100-d*40+(0*(a+b+c+d+e))}
spp915<-function(a,b,c,d,e) {100-d*50+(0*(a+b+c+d+e))}
spp916<-function(a,b,c,d,e) {100-d*100+(0*(a+b+c+d+e))}
spp917<-function(a,b,c,d,e) {100-d*110+(0*(a+b+c+d+e))}
spp918<-function(a,b,c,d,e) {100-d*120+(0*(a+b+c+d+e))}
spp919<-function(a,b,c,d,e) {100-d*140+(0*(a+b+c+d+e))}
spp920<-function(a,b,c,d,e) {100-d*150+(0*(a+b+c+d+e))}

spp921<-function(a,b,c,d,e) {500+d+(0*(a+b+c+d+e))}
spp922<-function(a,b,c,d,e) {500+d*2+(0*(a+b+c+d+e))}
spp923<-function(a,b,c,d,e) {500+d*3+(0*(a+b+c+d+e))}
spp924<-function(a,b,c,d,e) {500+d*4+(0*(a+b+c+d+e))}
spp925<-function(a,b,c,d,e) {500+d*5+(0*(a+b+c+d+e))}
spp926<-function(a,b,c,d,e) {500+d*6+(0*(a+b+c+d+e))}
spp927<-function(a,b,c,d,e) {500+d*7+(0*(a+b+c+d+e))}
spp928<-function(a,b,c,d,e) {500+d*8+(0*(a+b+c+d+e))}
spp929<-function(a,b,c,d,e) {500+d*9+(0*(a+b+c+d+e))}
spp930<-function(a,b,c,d,e) {500+d*10+(0*(a+b+c+d+e))}

spp931<-function(a,b,c,d,e) {500+d*11+(0*(a+b+c+d+e))}
spp932<-function(a,b,c,d,e) {500+d*20+(0*(a+b+c+d+e))}
spp933<-function(a,b,c,d,e) {500+d*30+(0*(a+b+c+d+e))}
spp934<-function(a,b,c,d,e) {500+d*40+(0*(a+b+c+d+e))}
spp935<-function(a,b,c,d,e) {500+d*50+(0*(a+b+c+d+e))}
spp936<-function(a,b,c,d,e) {500+d*100+(0*(a+b+c+d+e))}
spp937<-function(a,b,c,d,e) {500+d*110+(0*(a+b+c+d+e))}
spp938<-function(a,b,c,d,e) {500+d*120+(0*(a+b+c+d+e))}
spp939<-function(a,b,c,d,e) {500+d*140+(0*(a+b+c+d+e))}
spp940<-function(a,b,c,d,e) {500+d*150+(0*(a+b+c+d+e))}

spp941<-function(a,b,c,d,e) {500-d+(0*(a+b+c+d+e))}
spp942<-function(a,b,c,d,e) {500-d*2+(0*(a+b+c+d+e))}
spp943<-function(a,b,c,d,e) {500-d*3+(0*(a+b+c+d+e))}
spp944<-function(a,b,c,d,e) {500-d*4+(0*(a+b+c+d+e))}
spp945<-function(a,b,c,d,e) {500-d*5+(0*(a+b+c+d+e))}
spp946<-function(a,b,c,d,e) {500-d*6+(0*(a+b+c+d+e))}
spp947<-function(a,b,c,d,e) {500-d*7+(0*(a+b+c+d+e))}
spp948<-function(a,b,c,d,e) {500-d*8+(0*(a+b+c+d+e))}
spp949<-function(a,b,c,d,e) {500-d*9+(0*(a+b+c+d+e))}
spp950<-function(a,b,c,d,e) {500-d*10+(0*(a+b+c+d+e))}

spp951<-function(a,b,c,d,e) {500-d*11+(0*(a+b+c+d+e))}
spp952<-function(a,b,c,d,e) {500-d*20+(0*(a+b+c+d+e))}
spp953<-function(a,b,c,d,e) {500-d*30+(0*(a+b+c+d+e))}
spp954<-function(a,b,c,d,e) {500-d*40+(0*(a+b+c+d+e))}
spp955<-function(a,b,c,d,e) {500-d*50+(0*(a+b+c+d+e))}
spp956<-function(a,b,c,d,e) {500-d*100+(0*(a+b+c+d+e))}
spp957<-function(a,b,c,d,e) {500-d*110+(0*(a+b+c+d+e))}
spp958<-function(a,b,c,d,e) {500-d*120+(0*(a+b+c+d+e))}
spp959<-function(a,b,c,d,e) {500-d*140+(0*(a+b+c+d+e))}
spp960<-function(a,b,c,d,e) {500-d*150+(0*(a+b+c+d+e))}

##########################################
##                                      ##
##  linear Hi responses factor E        ##
##                                      ##
##########################################

spp961<-function(a,b,c,d,e) {1000+e*1000+(0*(a+b+c+d+e))}
spp962<-function(a,b,c,d,e) {5000+e*1000+(0*(a+b+c+d+e))}
spp963<-function(a,b,c,d,e) {10000+e*1000+(0*(a+b+c+d+e))}
spp964<-function(a,b,c,d,e) {15000+e*1000+(0*(a+b+c+d+e))}
spp965<-function(a,b,c,d,e) {20000+e*1000+(0*(a+b+c+d+e))}
spp966<-function(a,b,c,d,e) {25000+e*1000+(0*(a+b+c+d+e))}
spp967<-function(a,b,c,d,e) {30000+e*1000+(0*(a+b+c+d+e))}
spp968<-function(a,b,c,d,e) {35000+e*1000+(0*(a+b+c+d+e))}
spp969<-function(a,b,c,d,e) {40000+e*1000+(0*(a+b+c+d+e))}
spp970<-function(a,b,c,d,e) {50000+e*1000+(0*(a+b+c+d+e))}

spp971<-function(a,b,c,d,e) {1000-e*1000+(0*(a+b+c+d+e))}
spp972<-function(a,b,c,d,e) {5000-e*1000+(0*(a+b+c+d+e))}
spp973<-function(a,b,c,d,e) {10000-e*1000+(0*(a+b+c+d+e))}
spp974<-function(a,b,c,d,e) {15000-e*1000+(0*(a+b+c+d+e))}
spp975<-function(a,b,c,d,e) {20000-e*1000+(0*(a+b+c+d+e))}
spp976<-function(a,b,c,d,e) {25000-e*1000+(0*(a+b+c+d+e))}
spp977<-function(a,b,c,d,e) {30000-e*1000+(0*(a+b+c+d+e))}
spp978<-function(a,b,c,d,e) {35000-e*1000+(0*(a+b+c+d+e))}
spp979<-function(a,b,c,d,e) {40000-e*1000+(0*(a+b+c+d+e))}
spp980<-function(a,b,c,d,e) {50000-e*1000+(0*(a+b+c+d+e))}

spp981<-function(a,b,c,d,e) {10000+e*0.01+(0*(a+b+c+d+e))}
spp982<-function(a,b,c,d,e) {10000+e*0.1+(0*(a+b+c+d+e))}
spp983<-function(a,b,c,d,e) {10000+e+(0*(a+b+c+d+e))}
spp984<-function(a,b,c,d,e) {10000+e*10+(0*(a+b+c+d+e))}
spp985<-function(a,b,c,d,e) {10000+e*50+(0*(a+b+c+d+e))}
spp986<-function(a,b,c,d,e) {10000+e*100+(0*(a+b+c+d+e))}
spp987<-function(a,b,c,d,e) {10000+e*200+(0*(a+b+c+d+e))}
spp988<-function(a,b,c,d,e) {10000+e*500+(0*(a+b+c+d+e))}
spp989<-function(a,b,c,d,e) {10000+e*2000+(0*(a+b+c+d+e))}
spp990<-function(a,b,c,d,e) {10000+e*10000+(0*(a+b+c+d+e))}

spp991<-function(a,b,c,d,e) {10000-e*0.01+(0*(a+b+c+d+e))}
spp992<-function(a,b,c,d,e) {10000-e*0.1+(0*(a+b+c+d+e))}
spp993<-function(a,b,c,d,e) {10000-e+(0*(a+b+c+d+e))}
spp994<-function(a,b,c,d,e) {10000-e*10+(0*(a+b+c+d+e))}
spp995<-function(a,b,c,d,e) {10000-e*50+(0*(a+b+c+d+e))}
spp996<-function(a,b,c,d,e) {10000-e*100+(0*(a+b+c+d+e))}
spp997<-function(a,b,c,d,e) {10000-e*200+(0*(a+b+c+d+e))}
spp998<-function(a,b,c,d,e) {10000-e*500+(0*(a+b+c+d+e))}
spp999<-function(a,b,c,d,e) {10000-e*2000+(0*(a+b+c+d+e))}
spp1000<-function(a,b,c,d,e) {10000-e*10000+(0*(a+b+c+d+e))}

spp1001<-function(a,b,c,d,e) {1000+e*0.01+(0*(a+b+c+d+e))}
spp1002<-function(a,b,c,d,e) {2000+e*0.1+(0*(a+b+c+d+e))}
spp1003<-function(a,b,c,d,e) {5000+e+(0*(a+b+c+d+e))}
spp1004<-function(a,b,c,d,e) {10000+e*10+(0*(a+b+c+d+e))}
spp1005<-function(a,b,c,d,e) {15000+e*50+(0*(a+b+c+d+e))}
spp1006<-function(a,b,c,d,e) {20000+e*100+(0*(a+b+c+d+e))}
spp1007<-function(a,b,c,d,e) {25000+e*200+(0*(a+b+c+d+e))}
spp1008<-function(a,b,c,d,e) {30000+e*500+(0*(a+b+c+d+e))}
spp1009<-function(a,b,c,d,e) {40000+e*2000+(0*(a+b+c+d+e))}
spp1010<-function(a,b,c,d,e) {50000+e*10000+(0*(a+b+c+d+e))}

spp1011<-function(a,b,c,d,e) {50000-e*0.01+(0*(a+b+c+d+e))}
spp1012<-function(a,b,c,d,e) {40000-e*0.1+(0*(a+b+c+d+e))}
spp1013<-function(a,b,c,d,e) {30000-e+(0*(a+b+c+d+e))}
spp1014<-function(a,b,c,d,e) {25000-e*10+(0*(a+b+c+d+e))}
spp1015<-function(a,b,c,d,e) {20000-e*50+(0*(a+b+c+d+e))}
spp1016<-function(a,b,c,d,e) {15000-e*100+(0*(a+b+c+d+e))}
spp1017<-function(a,b,c,d,e) {10000-e*200+(0*(a+b+c+d+e))}
spp1018<-function(a,b,c,d,e) {5000-e*500+(0*(a+b+c+d+e))}
spp1019<-function(a,b,c,d,e) {2000-e*2000+(0*(a+b+c+d+e))}
spp1020<-function(a,b,c,d,e) {1000-e*10000+(0*(a+b+c+d+e))}

spp1021<-function(a,b,c,d,e) {1000-e*0.01+(0*(a+b+c+d+e))}
spp1022<-function(a,b,c,d,e) {2000-e*0.1+(0*(a+b+c+d+e))}
spp1023<-function(a,b,c,d,e) {5000-e+(0*(a+b+c+d+e))}
spp1024<-function(a,b,c,d,e) {10000-e*10+(0*(a+b+c+d+e))}
spp1025<-function(a,b,c,d,e) {15000-e*50+(0*(a+b+c+d+e))}
spp1026<-function(a,b,c,d,e) {20000-e*100+(0*(a+b+c+d+e))}
spp1027<-function(a,b,c,d,e) {25000-e*200+(0*(a+b+c+d+e))}
spp1028<-function(a,b,c,d,e) {30000-e*500+(0*(a+b+c+d+e))}
spp1029<-function(a,b,c,d,e) {40000-e*2000+(0*(a+b+c+d+e))}
spp1030<-function(a,b,c,d,e) {50000-e*10000+(0*(a+b+c+d+e))}

spp1031<-function(a,b,c,d,e) {50000+e*0.01+(0*(a+b+c+d+e))}
spp1032<-function(a,b,c,d,e) {40000+e*0.1+(0*(a+b+c+d+e))}
spp1033<-function(a,b,c,d,e) {30000+e+(0*(a+b+c+d+e))}
spp1034<-function(a,b,c,d,e) {25000+e*10+(0*(a+b+c+d+e))}
spp1035<-function(a,b,c,d,e) {20000+e*50+(0*(a+b+c+d+e))}
spp1036<-function(a,b,c,d,e) {15000+e*100+(0*(a+b+c+d+e))}
spp1037<-function(a,b,c,d,e) {10000+e*200+(0*(a+b+c+d+e))}
spp1038<-function(a,b,c,d,e) {5000+e*500+(0*(a+b+c+d+e))}
spp1039<-function(a,b,c,d,e) {2000+e*2000+(0*(a+b+c+d+e))}
spp1040<-function(a,b,c,d,e) {1000+e*10000+(0*(a+b+c+d+e))}


##########################################
##                                      ##
##  linear low responses factor E       ##
##                                      ##
##########################################

spp1041<-function(a,b,c,d,e) {10+e+(0*(a+b+c+d+e))}
spp1042<-function(a,b,c,d,e) {10+e*2+(0*(a+b+c+d+e))}
spp1043<-function(a,b,c,d,e) {10+e*3+(0*(a+b+c+d+e))}
spp1044<-function(a,b,c,d,e) {10+e*4+(0*(a+b+c+d+e))}
spp1045<-function(a,b,c,d,e) {10+e*5+(0*(a+b+c+d+e))}
spp1046<-function(a,b,c,d,e) {10+e*6+(0*(a+b+c+d+e))}
spp1047<-function(a,b,c,d,e) {10+e*7+(0*(a+b+c+d+e))}
spp1048<-function(a,b,c,d,e) {10+e*8+(0*(a+b+c+d+e))}
spp1049<-function(a,b,c,d,e) {10+e*9+(0*(a+b+c+d+e))}
spp1050<-function(a,b,c,d,e) {10+e*10+(0*(a+b+c+d+e))}

spp1051<-function(a,b,c,d,e) {10+e*11+(0*(a+b+c+d+e))}
spp1052<-function(a,b,c,d,e) {10+e*20+(0*(a+b+c+d+e))}
spp1053<-function(a,b,c,d,e) {10+e*30+(0*(a+b+c+d+e))}
spp1054<-function(a,b,c,d,e) {10+e*40+(0*(a+b+c+d+e))}
spp1055<-function(a,b,c,d,e) {10+e*50+(0*(a+b+c+d+e))}
spp1056<-function(a,b,c,d,e) {10+e*100+(0*(a+b+c+d+e))}
spp1057<-function(a,b,c,d,e) {10+e*110+(0*(a+b+c+d+e))}
spp1058<-function(a,b,c,d,e) {10+e*120+(0*(a+b+c+d+e))}
spp1059<-function(a,b,c,d,e) {10+e*140+(0*(a+b+c+d+e))}
spp1060<-function(a,b,c,d,e) {10+e*150+(0*(a+b+c+d+e))}

spp1061<-function(a,b,c,d,e) {10-e+(0*(a+b+c+d+e))}
spp1062<-function(a,b,c,d,e) {10-e*2+(0*(a+b+c+d+e))}
spp1063<-function(a,b,c,d,e) {10-e*3+(0*(a+b+c+d+e))}
spp1064<-function(a,b,c,d,e) {10-e*4+(0*(a+b+c+d+e))}
spp1065<-function(a,b,c,d,e) {10-e*5+(0*(a+b+c+d+e))}
spp1066<-function(a,b,c,d,e) {10-e*6+(0*(a+b+c+d+e))}
spp1067<-function(a,b,c,d,e) {10-e*7+(0*(a+b+c+d+e))}
spp1068<-function(a,b,c,d,e) {10-e*8+(0*(a+b+c+d+e))}
spp1069<-function(a,b,c,d,e) {10-e*9+(0*(a+b+c+d+e))}
spp1070<-function(a,b,c,d,e) {10-e*10+(0*(a+b+c+d+e))}

spp1071<-function(a,b,c,d,e) {10-e*11+(0*(a+b+c+d+e))}
spp1072<-function(a,b,c,d,e) {10-e*20+(0*(a+b+c+d+e))}
spp1073<-function(a,b,c,d,e) {10-e*30+(0*(a+b+c+d+e))}
spp1074<-function(a,b,c,d,e) {10-e*40+(0*(a+b+c+d+e))}
spp1075<-function(a,b,c,d,e) {10-e*50+(0*(a+b+c+d+e))}
spp1076<-function(a,b,c,d,e) {10-e*100+(0*(a+b+c+d+e))}
spp1077<-function(a,b,c,d,e) {10-e*110+(0*(a+b+c+d+e))}
spp1078<-function(a,b,c,d,e) {10-e*120+(0*(a+b+c+d+e))}
spp1079<-function(a,b,c,d,e) {10-e*140+(0*(a+b+c+d+e))}
spp1080<-function(a,b,c,d,e) {10-e*150+(0*(a+b+c+d+e))}

spp1081<-function(a,b,c,d,e) {50+e+(0*(a+b+c+d+e))}
spp1082<-function(a,b,c,d,e) {50+e*2+(0*(a+b+c+d+e))}
spp1083<-function(a,b,c,d,e) {50+e*3+(0*(a+b+c+d+e))}
spp1084<-function(a,b,c,d,e) {50+e*4+(0*(a+b+c+d+e))}
spp1085<-function(a,b,c,d,e) {50+e*5+(0*(a+b+c+d+e))}
spp1086<-function(a,b,c,d,e) {50+e*6+(0*(a+b+c+d+e))}
spp1087<-function(a,b,c,d,e) {50+e*7+(0*(a+b+c+d+e))}
spp1088<-function(a,b,c,d,e) {50+e*8+(0*(a+b+c+d+e))}
spp1089<-function(a,b,c,d,e) {50+e*9+(0*(a+b+c+d+e))}
spp1090<-function(a,b,c,d,e) {50+e*10+(0*(a+b+c+d+e))}

spp1091<-function(a,b,c,d,e) {50+e*11+(0*(a+b+c+d+e))}
spp1092<-function(a,b,c,d,e) {50+e*20+(0*(a+b+c+d+e))}
spp1093<-function(a,b,c,d,e) {50+e*30+(0*(a+b+c+d+e))}
spp1094<-function(a,b,c,d,e) {50+e*40+(0*(a+b+c+d+e))}
spp1095<-function(a,b,c,d,e) {50+e*50+(0*(a+b+c+d+e))}
spp1096<-function(a,b,c,d,e) {50+e*100+(0*(a+b+c+d+e))}
spp1097<-function(a,b,c,d,e) {50+e*110+(0*(a+b+c+d+e))}
spp1098<-function(a,b,c,d,e) {50+e*120+(0*(a+b+c+d+e))}
spp1099<-function(a,b,c,d,e) {50+e*140+(0*(a+b+c+d+e))}
spp1100<-function(a,b,c,d,e) {50+e*150+(0*(a+b+c+d+e))}

spp1101<-function(a,b,c,d,e) {50-e+(0*(a+b+c+d+e))}
spp1102<-function(a,b,c,d,e) {50-e*2+(0*(a+b+c+d+e))}
spp1103<-function(a,b,c,d,e) {50-e*3+(0*(a+b+c+d+e))}
spp1104<-function(a,b,c,d,e) {50-e*4+(0*(a+b+c+d+e))}
spp1105<-function(a,b,c,d,e) {50-e*5+(0*(a+b+c+d+e))}
spp1106<-function(a,b,c,d,e) {50-e*6+(0*(a+b+c+d+e))}
spp1107<-function(a,b,c,d,e) {50-e*7+(0*(a+b+c+d+e))}
spp1108<-function(a,b,c,d,e) {50-e*8+(0*(a+b+c+d+e))}
spp1109<-function(a,b,c,d,e) {50-e*9+(0*(a+b+c+d+e))}
spp1110<-function(a,b,c,d,e) {50-e*10+(0*(a+b+c+d+e))}

spp1111<-function(a,b,c,d,e) {50-e*11+(0*(a+b+c+d+e))}
spp1112<-function(a,b,c,d,e) {50-e*20+(0*(a+b+c+d+e))}
spp1113<-function(a,b,c,d,e) {50-e*30+(0*(a+b+c+d+e))}
spp1114<-function(a,b,c,d,e) {50-e*40+(0*(a+b+c+d+e))}
spp1115<-function(a,b,c,d,e) {50-e*50+(0*(a+b+c+d+e))}
spp1116<-function(a,b,c,d,e) {50-e*100+(0*(a+b+c+d+e))}
spp1117<-function(a,b,c,d,e) {50-e*110+(0*(a+b+c+d+e))}
spp1118<-function(a,b,c,d,e) {50-e*120+(0*(a+b+c+d+e))}
spp1119<-function(a,b,c,d,e) {50-e*140+(0*(a+b+c+d+e))}
spp1120<-function(a,b,c,d,e) {50-e*150+(0*(a+b+c+d+e))}

##########################################
##                                      ##
##  linear medium responses factor E    ##
##                                      ##
##########################################

spp1121<-function(a,b,c,d,e) {100+e+(0*(a+b+c+d+e))}
spp1122<-function(a,b,c,d,e) {100+e*2+(0*(a+b+c+d+e))}
spp1123<-function(a,b,c,d,e) {100+e*3+(0*(a+b+c+d+e))}
spp1124<-function(a,b,c,d,e) {100+e*4+(0*(a+b+c+d+e))}
spp1125<-function(a,b,c,d,e) {100+e*5+(0*(a+b+c+d+e))}
spp1126<-function(a,b,c,d,e) {100+e*6+(0*(a+b+c+d+e))}
spp1127<-function(a,b,c,d,e) {100+e*7+(0*(a+b+c+d+e))}
spp1128<-function(a,b,c,d,e) {100+e*8+(0*(a+b+c+d+e))}
spp1129<-function(a,b,c,d,e) {100+e*9+(0*(a+b+c+d+e))}
spp1130<-function(a,b,c,d,e) {100+e*10+(0*(a+b+c+d+e))}

spp1131<-function(a,b,c,d,e) {100+e*11+(0*(a+b+c+d+e))}
spp1132<-function(a,b,c,d,e) {100+e*20+(0*(a+b+c+d+e))}
spp1133<-function(a,b,c,d,e) {100+e*30+(0*(a+b+c+d+e))}
spp1134<-function(a,b,c,d,e) {100+e*40+(0*(a+b+c+d+e))}
spp1135<-function(a,b,c,d,e) {100+e*50+(0*(a+b+c+d+e))}
spp1136<-function(a,b,c,d,e) {100+e*100+(0*(a+b+c+d+e))}
spp1137<-function(a,b,c,d,e) {100+e*110+(0*(a+b+c+d+e))}
spp1138<-function(a,b,c,d,e) {100+e*120+(0*(a+b+c+d+e))}
spp1139<-function(a,b,c,d,e) {100+e*140+(0*(a+b+c+d+e))}
spp1140<-function(a,b,c,d,e) {100+e*150+(0*(a+b+c+d+e))}

spp1141<-function(a,b,c,d,e) {100-e+(0*(a+b+c+d+e))}
spp1142<-function(a,b,c,d,e) {100-e*2+(0*(a+b+c+d+e))}
spp1143<-function(a,b,c,d,e) {100-e*3+(0*(a+b+c+d+e))}
spp1144<-function(a,b,c,d,e) {100-e*4+(0*(a+b+c+d+e))}
spp1145<-function(a,b,c,d,e) {100-e*5+(0*(a+b+c+d+e))}
spp1146<-function(a,b,c,d,e) {100-e*6+(0*(a+b+c+d+e))}
spp1147<-function(a,b,c,d,e) {100-e*7+(0*(a+b+c+d+e))}
spp1148<-function(a,b,c,d,e) {100-e*8+(0*(a+b+c+d+e))}
spp1149<-function(a,b,c,d,e) {100-e*9+(0*(a+b+c+d+e))}
spp1150<-function(a,b,c,d,e) {100-e*10+(0*(a+b+c+d+e))}

spp1151<-function(a,b,c,d,e) {100-e*11+(0*(a+b+c+d+e))}
spp1152<-function(a,b,c,d,e) {100-e*20+(0*(a+b+c+d+e))}
spp1153<-function(a,b,c,d,e) {100-e*30+(0*(a+b+c+d+e))}
spp1154<-function(a,b,c,d,e) {100-e*40+(0*(a+b+c+d+e))}
spp1155<-function(a,b,c,d,e) {100-e*50+(0*(a+b+c+d+e))}
spp1156<-function(a,b,c,d,e) {100-e*100+(0*(a+b+c+d+e))}
spp1157<-function(a,b,c,d,e) {100-e*110+(0*(a+b+c+d+e))}
spp1158<-function(a,b,c,d,e) {100-e*120+(0*(a+b+c+d+e))}
spp1159<-function(a,b,c,d,e) {100-e*140+(0*(a+b+c+d+e))}
spp1160<-function(a,b,c,d,e) {100-e*150+(0*(a+b+c+d+e))}

spp1161<-function(a,b,c,d,e) {500+e+(0*(a+b+c+d+e))}
spp1162<-function(a,b,c,d,e) {500+e*2+(0*(a+b+c+d+e))}
spp1163<-function(a,b,c,d,e) {500+e*3+(0*(a+b+c+d+e))}
spp1164<-function(a,b,c,d,e) {500+e*4+(0*(a+b+c+d+e))}
spp1165<-function(a,b,c,d,e) {500+e*5+(0*(a+b+c+d+e))}
spp1166<-function(a,b,c,d,e) {500+e*6+(0*(a+b+c+d+e))}
spp1167<-function(a,b,c,d,e) {500+e*7+(0*(a+b+c+d+e))}
spp1168<-function(a,b,c,d,e) {500+e*8+(0*(a+b+c+d+e))}
spp1169<-function(a,b,c,d,e) {500+e*9+(0*(a+b+c+d+e))}
spp1170<-function(a,b,c,d,e) {500+e*10+(0*(a+b+c+d+e))}

spp1171<-function(a,b,c,d,e) {500+e*11+(0*(a+b+c+d+e))}
spp1172<-function(a,b,c,d,e) {500+e*20+(0*(a+b+c+d+e))}
spp1173<-function(a,b,c,d,e) {500+e*30+(0*(a+b+c+d+e))}
spp1174<-function(a,b,c,d,e) {500+e*40+(0*(a+b+c+d+e))}
spp1175<-function(a,b,c,d,e) {500+e*50+(0*(a+b+c+d+e))}
spp1176<-function(a,b,c,d,e) {500+e*100+(0*(a+b+c+d+e))}
spp1177<-function(a,b,c,d,e) {500+e*110+(0*(a+b+c+d+e))}
spp1178<-function(a,b,c,d,e) {500+e*120+(0*(a+b+c+d+e))}
spp1179<-function(a,b,c,d,e) {500+e*140+(0*(a+b+c+d+e))}
spp1180<-function(a,b,c,d,e) {500+e*150+(0*(a+b+c+d+e))}

spp1181<-function(a,b,c,d,e) {500-e+(0*(a+b+c+d+e))}
spp1182<-function(a,b,c,d,e) {500-e*2+(0*(a+b+c+d+e))}
spp1183<-function(a,b,c,d,e) {500-e*3+(0*(a+b+c+d+e))}
spp1184<-function(a,b,c,d,e) {500-e*4+(0*(a+b+c+d+e))}
spp1185<-function(a,b,c,d,e) {500-e*5+(0*(a+b+c+d+e))}
spp1186<-function(a,b,c,d,e) {500-e*6+(0*(a+b+c+d+e))}
spp1187<-function(a,b,c,d,e) {500-e*7+(0*(a+b+c+d+e))}
spp1188<-function(a,b,c,d,e) {500-e*8+(0*(a+b+c+d+e))}
spp1189<-function(a,b,c,d,e) {500-e*9+(0*(a+b+c+d+e))}
spp1190<-function(a,b,c,d,e) {500-e*10+(0*(a+b+c+d+e))}

spp1191<-function(a,b,c,d,e) {500-e*11+(0*(a+b+c+d+e))}
spp1192<-function(a,b,c,d,e) {500-e*20+(0*(a+b+c+d+e))}
spp1193<-function(a,b,c,d,e) {500-e*30+(0*(a+b+c+d+e))}
spp1194<-function(a,b,c,d,e) {500-e*40+(0*(a+b+c+d+e))}
spp1195<-function(a,b,c,d,e) {500-e*50+(0*(a+b+c+d+e))}
spp1196<-function(a,b,c,d,e) {500-e*100+(0*(a+b+c+d+e))}
spp1197<-function(a,b,c,d,e) {500-e*110+(0*(a+b+c+d+e))}
spp1198<-function(a,b,c,d,e) {500-e*120+(0*(a+b+c+d+e))}
spp1199<-function(a,b,c,d,e) {500-e*140+(0*(a+b+c+d+e))}
spp1200<-function(a,b,c,d,e) {500-e*150+(0*(a+b+c+d+e))}

################################################
##                                            ##
##  multiple linear low responses             ##
##                                            ##
################################################
spp1202<-function(a,b,c,d,e) {10*a+b+c+d+e}
spp1203<-function(a,b,c,d,e) {a+10*b+c+d+e}
spp1204<-function(a,b,c,d,e) {a+b+10*c+d+e}
spp1205<-function(a,b,c,d,e) {a+b+c+10*d+e}
spp1206<-function(a,b,c,d,e) {a+b+c+d+10*e}
spp1207<-function(a,b,c,d,e) {100-(10*a+b+c+d+e)}
spp1208<-function(a,b,c,d,e) {100-(a+10*b+c+d+e)}
spp1209<-function(a,b,c,d,e) {100-(a+b+10*c+d+e)}
spp1210<-function(a,b,c,d,e) {100-(a+b+c+10*d+e)}
spp1201<-function(a,b,c,d,e) {100-(a+b+c+d+10*e)}


spp1211<-function(a,b,c,d,e) {50+(a*b+c+d+e)}
spp1212<-function(a,b,c,d,e) {50+(a*b-c+d+e)}
spp1213<-function(a,b,c,d,e) {50-(a*b-c+d+e)}
spp1214<-function(a,b,c,d,e) {50-(a*b+c+d+e)}
spp1215<-function(a,b,c,d,e) {50-(a+b+c+d*e)}
spp1216<-function(a,b,c,d,e) {50-(a+b+c-d*e)}
spp1217<-function(a,b,c,d,e) {50+(a+b+c+d*e)}
spp1218<-function(a,b,c,d,e) {50+(a+b+c-d*e)}
spp1219<-function(a,b,c,d,e) {50-(a+b*c+d*e)}
spp1220<-function(a,b,c,d,e) {50+(a+b*c+d*e)}

spp1221<-function(a,b,c,d,e) {50+(a*b+10*c+d+e)}
spp1222<-function(a,b,c,d,e) {50+(a*b-10*c+d+e)}
spp1223<-function(a,b,c,d,e) {50-(a*b-c+10*d+e)}
spp1224<-function(a,b,c,d,e) {50-(a*b+c+10*d+e)}
spp1225<-function(a,b,c,d,e) {50-(a+10*b+c+d*e)}
spp1226<-function(a,b,c,d,e) {50-(a+10*b+c-d*e)}
spp1227<-function(a,b,c,d,e) {50+(a+b+c+10*d*e)}
spp1228<-function(a,b,c,d,e) {50+(a+b+c-10*d*e)}
spp1229<-function(a,b,c,d,e) {round(50-(a+b*c+d*e/10))}
spp1230<-function(a,b,c,d,e) {round(50+(a+b*c+d*e/10))}

spp1231<-function(a,b,c,d,e) {round(1+(a*b+c+d+e))}
spp1232<-function(a,b,c,d,e) {round(1-(a*b+c+d+e))}
spp1233<-function(a,b,c,d,e) {round(1+(a*b*c+d+e))}
spp1234<-function(a,b,c,d,e) {round(1-(a*b*c+d+e))}
spp1235<-function(a,b,c,d,e) {round(1+(a*b+c+d+e)/10)}
spp1236<-function(a,b,c,d,e) {round(1+(a*b+c+d+e)/20)}
spp1237<-function(a,b,c,d,e) {round(1+(a+b+c+d+e)/10)}
spp1238<-function(a,b,c,d,e) {round(1-(a+b+c+d+e)/10)}
spp1239<-function(a,b,c,d,e) {round(1+(a*b*c*d*e)/10)}
spp1240<-function(a,b,c,d,e) {round(1-(a*b*c*d*e)/10)}

spp1241<-function(a,b,c,d,e) {round(80+(a*b+c+d+e))}
spp1242<-function(a,b,c,d,e) {round(80-(a*b+c+d+e))}
spp1243<-function(a,b,c,d,e) {round(80+(a*b*c+d+e))}
spp1244<-function(a,b,c,d,e) {round(80-(a*b*c+d+e))}
spp1245<-function(a,b,c,d,e) {round(80+(a*b+c+d+e)/10)}
spp1246<-function(a,b,c,d,e) {round(80+(a*b+c+d+e)/20)}
spp1247<-function(a,b,c,d,e) {round(80+(a+b+c+d+e)/10)}
spp1248<-function(a,b,c,d,e) {round(80-(a+b+c+d+e)/10)}
spp1249<-function(a,b,c,d,e) {round(80+(a*b*c*d*e)/10)}
spp1250<-function(a,b,c,d,e) {round(80-(a*b*c*d*e)/10)}

spp1251<-function(a,b,c,d,e) {round(1+(a*b-c+d+e))}
spp1252<-function(a,b,c,d,e) {round(1-(a*b-c+d+e))}
spp1253<-function(a,b,c,d,e) {round(1+(a*b*c-d+e))}
spp1254<-function(a,b,c,d,e) {round(1-(a*b*c-d+e))}
spp1255<-function(a,b,c,d,e) {round(1-(a-b-c-d-e)/10)}
spp1256<-function(a,b,c,d,e) {round(1+(a-b-c-d-e)/20)}
spp1257<-function(a,b,c,d,e) {round(1+(a-b-c-d-e)/10)}
spp1258<-function(a,b,c,d,e) {round(1+(a-10*b+c+d+e)/10)}
spp1259<-function(a,b,c,d,e) {round(1+(a-10*b*c*d*e)/10)}
spp1260<-function(a,b,c,d,e) {round(1+(a*b*c+d-10*e)/10)}

spp1261<-function(a,b,c,d,e) {round(1+(a+b+c+d+10*e))}
spp1262<-function(a,b,c,d,e) {round(1+(a+b+c+d+50*e))}
spp1263<-function(a,b,c,d,e) {round(1+(a+b+c+d+100*e))}
spp1264<-function(a,b,c,d,e) {round(1+(a+b+c+10*d+e))}
spp1265<-function(a,b,c,d,e) {round(1+(a+b+c+50*d+e))}
spp1266<-function(a,b,c,d,e) {round(1+(a+b+c+100*d+e))}
spp1267<-function(a,b,c,d,e) {round(1+(-10*a+b+c+d+10*e))}
spp1268<-function(a,b,c,d,e) {round(1+(-50*a+b+c+d+50*e))}
spp1269<-function(a,b,c,d,e) {round(1+(a+b-10*c+10*d+e))}
spp1270<-function(a,b,c,d,e) {round(1+(a+b-50*c+50*d+e))}

spp1271<-function(a,b,c,d,e) {10*a*b*c*d*e}
spp1272<-function(a,b,c,d,e) {1+(a+b-10*c+50*d+e)}
spp1273<-function(a,b,c,d,e) {1+(a+b-10*c+100*d+e)}
spp1274<-function(a,b,c,d,e) {1+(a+b-50*c+10*d+e)}
spp1275<-function(a,b,c,d,e) {1+(a+b-100*c+50*d+e)}
spp1276<-function(a,b,c,d,e) {1+(a+10*b+c-50*d+e)}
spp1277<-function(a,b,c,d,e) {1+(a+50*b+c-100*d+e)}
spp1278<-function(a,b,c,d,e) {1+(a+100*b+c-10*d+e)}
spp1279<-function(a,b,c,d,e) {1+(a+b-50*c+d+100*e)}
spp1280<-function(a,b,c,d,e) {1+(a+b-50*c+d+10*e)}

spp1281<-function(a,b,c,d,e) {100-(a+b-10*c+50*d+e)}
spp1282<-function(a,b,c,d,e) {100-(a+b-10*c+100*d+e)}
spp1283<-function(a,b,c,d,e) {100-(a+b-50*c+10*d+e)}
spp1284<-function(a,b,c,d,e) {100-(a+b-100*c+50*d+e)}
spp1285<-function(a,b,c,d,e) {100-(a+10*b+c-50*d+e)}
spp1286<-function(a,b,c,d,e) {100-(a+50*b+c-100*d+e)}
spp1287<-function(a,b,c,d,e) {100-(a+100*b+c-10*d+e)}
spp1288<-function(a,b,c,d,e) {100-(a+b-50*c+d+100*e)}
spp1289<-function(a,b,c,d,e) {100-(a+b-50*c+d+10*e)}
spp1290<-function(a,b,c,d,e) {100-(a+b-50*c*d+10*e)}
################################################
##                                            ##
##  multiple linear med responses             ##
##                                            ##
################################################
spp1292<-function(a,b,c,d,e) {150*a+b+c+d+e}
spp1293<-function(a,b,c,d,e) {a+150*b+c+d+e}
spp1294<-function(a,b,c,d,e) {a+b+150*c+d+e}
spp1295<-function(a,b,c,d,e) {a+b+c+150*d+e}
spp1296<-function(a,b,c,d,e) {a+b+c+d+150*e}
spp1297<-function(a,b,c,d,e) {500-(10*a+b+c+d+e)}
spp1298<-function(a,b,c,d,e) {500-(a+10*b+c+d+e)}
spp1299<-function(a,b,c,d,e) {500-(a+b+10*c+d+e)}
spp1300<-function(a,b,c,d,e) {500-(a+b+c+10*d+e)}
spp1291<-function(a,b,c,d,e) {500-(a+b+c+d+10*e)}


spp1301<-function(a,b,c,d,e) {500+(a*b+c+d+e)}
spp1302<-function(a,b,c,d,e) {500+(a*b-c+d+e)}
spp1303<-function(a,b,c,d,e) {500-(a*b-c+d+e)}
spp1304<-function(a,b,c,d,e) {500-(a*b+c+d+e)}
spp1305<-function(a,b,c,d,e) {500-(a+b+c+d*e)}
spp1306<-function(a,b,c,d,e) {500-(a+b+c-d*e)}
spp1307<-function(a,b,c,d,e) {500+(a+b+c+d*e)}
spp1308<-function(a,b,c,d,e) {500+(a+b+c-d*e)}
spp1309<-function(a,b,c,d,e) {500-(a+b*c+d*e)}
spp1310<-function(a,b,c,d,e) {500+(a+b*c+d*e)}

spp1311<-function(a,b,c,d,e) {500+(a*b+10*c+d+e)}
spp1312<-function(a,b,c,d,e) {500+(a*b-10*c+d+e)}
spp1313<-function(a,b,c,d,e) {500-(a*b-c+10*d+e)}
spp1314<-function(a,b,c,d,e) {500-(a*b+c+10*d+e)}
spp1315<-function(a,b,c,d,e) {500-(a+10*b+c+d*e)}
spp1316<-function(a,b,c,d,e) {500-(a+10*b+c-d*e)}
spp1317<-function(a,b,c,d,e) {500+(a+b+c+10*d*e)}
spp1318<-function(a,b,c,d,e) {500+(a+b+c-10*d*e)}
spp1319<-function(a,b,c,d,e) {round(500-(a+b*c+d*e/10))}
spp1320<-function(a,b,c,d,e) {round(500+(a+b*c+d*e/10))}

spp1321<-function(a,b,c,d,e) {round(100+(a*b+c+d+e))}
spp1322<-function(a,b,c,d,e) {round(100-(a*b+c+d+e))}
spp1323<-function(a,b,c,d,e) {round(100+(a*b*c+d+e))}
spp1324<-function(a,b,c,d,e) {round(100-(a*b*c+d+e))}
spp1325<-function(a,b,c,d,e) {round(100+(a*b+c+d+e)/10)}
spp1326<-function(a,b,c,d,e) {round(100+(a*b+c+d+e)/20)}
spp1327<-function(a,b,c,d,e) {round(100+(a+b+c+d+e)/10)}
spp1328<-function(a,b,c,d,e) {round(100-(a+b+c+d+e)/10)}
spp1329<-function(a,b,c,d,e) {round(100+(a*b*c*d*e)/10)}
spp1330<-function(a,b,c,d,e) {round(100-(a*b*c*d*e)/10)}

spp1331<-function(a,b,c,d,e) {round(800+(a*b+c+d+e))}
spp1332<-function(a,b,c,d,e) {round(800-(a*b+c+d+e))}
spp1333<-function(a,b,c,d,e) {round(800+(a*b*c+d+e))}
spp1334<-function(a,b,c,d,e) {round(800-(a*b*c+d+e))}
spp1335<-function(a,b,c,d,e) {round(800+(a*b+c+d+e)/10)}
spp1336<-function(a,b,c,d,e) {round(800+(a*b+c+d+e)/20)}
spp1337<-function(a,b,c,d,e) {round(800+(a+b+c+d+e)/10)}
spp1338<-function(a,b,c,d,e) {round(800-(a+b+c+d+e)/10)}
spp1339<-function(a,b,c,d,e) {round(800+(a*b*c*d*e)/10)}
spp1340<-function(a,b,c,d,e) {round(800-(a*b*c*d*e)/10)}

spp1341<-function(a,b,c,d,e) {round(500+(a*b-c+d+e))}
spp1342<-function(a,b,c,d,e) {round(500-(a*b-c+d+e))}
spp1343<-function(a,b,c,d,e) {round(500+(a*b*c-d+e))}
spp1344<-function(a,b,c,d,e) {round(500-(a*b*c-d+e))}
spp1345<-function(a,b,c,d,e) {round(500-(a-b-c-d-e)/10)}
spp1346<-function(a,b,c,d,e) {round(500+(a-b-c-d-e)/20)}
spp1347<-function(a,b,c,d,e) {round(500+(a-b-c-d-e)/10)}
spp1348<-function(a,b,c,d,e) {round(500+(a-10*b+c+d+e)/10)}
spp1349<-function(a,b,c,d,e) {round(500+(a-10*b*c*d*e)/10)}
spp1350<-function(a,b,c,d,e) {round(500+(a*b*c+d-10*e)/10)}

spp1351<-function(a,b,c,d,e) {round(1000+(a+b+c+d+10*e))}
spp1352<-function(a,b,c,d,e) {round(1000+(a+b+c+d+50*e))}
spp1353<-function(a,b,c,d,e) {round(1000+(a+b+c+d+100*e))}
spp1354<-function(a,b,c,d,e) {round(1000+(a+b+c+10*d+e))}
spp1355<-function(a,b,c,d,e) {round(1000+(a+b+c+50*d+e))}
spp1356<-function(a,b,c,d,e) {round(1000+(a+b+c+100*d+e))}
spp1357<-function(a,b,c,d,e) {round(1000+(-10*a+b+c+d+10*e))}
spp1358<-function(a,b,c,d,e) {round(1000+(-50*a+b+c+d+50*e))}
spp1359<-function(a,b,c,d,e) {round(1000+(a+b-10*c+10*d+e))}
spp1360<-function(a,b,c,d,e) {round(1000+(a+b-50*c+50*d+e))}

spp1361<-function(a,b,c,d,e) {100*a*b*c*d*e}
spp1362<-function(a,b,c,d,e) {600+(a+b-10*c+50*d+e)}
spp1363<-function(a,b,c,d,e) {600+(a+b-10*c+100*d+e)}
spp1364<-function(a,b,c,d,e) {600+(a+b-50*c+10*d+e)}
spp1365<-function(a,b,c,d,e) {600+(a+b-100*c+50*d+e)}
spp1366<-function(a,b,c,d,e) {600+(a+10*b+c-50*d+e)}
spp1367<-function(a,b,c,d,e) {600+(a+50*b+c-100*d+e)}
spp1368<-function(a,b,c,d,e) {600+(a+100*b+c-10*d+e)}
spp1369<-function(a,b,c,d,e) {600+(a+b-50*c+d+100*e)}
spp1370<-function(a,b,c,d,e) {600+(a+b-50*c+d+10*e)}

spp1372<-function(a,b,c,d,e) {1000-(a+b-10*c+50*d+e)}
spp1373<-function(a,b,c,d,e) {1000-(a+b-10*c+100*d+e)}
spp1374<-function(a,b,c,d,e) {1000-(a+b-50*c+10*d+e)}
spp1375<-function(a,b,c,d,e) {1000-(a+b-100*c+50*d+e)}
spp1376<-function(a,b,c,d,e) {1000-(a+10*b+c-50*d+e)}
spp1377<-function(a,b,c,d,e) {1000-(a+50*b+c-100*d+e)}
spp1378<-function(a,b,c,d,e) {1000-(a+100*b+c-10*d+e)}
spp1379<-function(a,b,c,d,e) {1000-(a+b-50*c+d+100*e)}
spp1380<-function(a,b,c,d,e) {1000-(a+b-50*c+d+10*e)}
spp1371<-function(a,b,c,d,e) {1000-(a+b-50*c*d+10*e)}

################################################
##                                            ##
##  multiple linear hi  responses             ##
##                                            ##
################################################
spp1382<-function(a,b,c,d,e) {1500*a+b+c+d+e}
spp1383<-function(a,b,c,d,e) {a+1500*b+c+d+e}
spp1384<-function(a,b,c,d,e) {a+b+1500*c+d+e}
spp1385<-function(a,b,c,d,e) {a+b+c+1500*d+e}
spp1386<-function(a,b,c,d,e) {a+b+c+d+1500*e}
spp1387<-function(a,b,c,d,e) {5000-(100*a+b+c+d+e)}
spp1388<-function(a,b,c,d,e) {5000-(a+100*b+c+d+e)}
spp1389<-function(a,b,c,d,e) {5000-(a+b+100*c+d+e)}
spp1390<-function(a,b,c,d,e) {5000-(a+b+c+100*d+e)}
spp1381<-function(a,b,c,d,e) {5000-(a+b+c+d+100*e)}


spp1391<-function(a,b,c,d,e) {5000+(a*b+c+d+e)}
spp1392<-function(a,b,c,d,e) {5000+(a*b-c+d+e)}
spp1393<-function(a,b,c,d,e) {5000-(a*b-c+d+e)}
spp1394<-function(a,b,c,d,e) {5000-(a*b+c+d+e)}
spp1395<-function(a,b,c,d,e) {5000-(a+b+c+d*e)}
spp1396<-function(a,b,c,d,e) {5000-(a+b+c-d*e)}
spp1397<-function(a,b,c,d,e) {5000+(a+b+c+d*e)}
spp1398<-function(a,b,c,d,e) {5000+(a+b+c-d*e)}
spp1399<-function(a,b,c,d,e) {5000-(a+b*c+d*e)}
spp1400<-function(a,b,c,d,e) {5000+(a+b*c+d*e)}

spp1401<-function(a,b,c,d,e) {5000+(a*b+100*c+d+e)}
spp1402<-function(a,b,c,d,e) {5000+(a*b-100*c+d+e)}
spp1403<-function(a,b,c,d,e) {5000-(a*b-c+100*d+e)}
spp1404<-function(a,b,c,d,e) {5000-(a*b+c+100*d+e)}
spp1405<-function(a,b,c,d,e) {5000-(a+100*b+c+d*e)}
spp1406<-function(a,b,c,d,e) {5000-(a+100*b+c-d*e)}
spp1407<-function(a,b,c,d,e) {5000+(a+b+c+100*d*e)}
spp1408<-function(a,b,c,d,e) {5000+(a+b+c-100*d*e)}
spp1409<-function(a,b,c,d,e) {round(5000-(a+100*b*c+d*e/10))}
spp1410<-function(a,b,c,d,e) {round(5000+(a+100*b*c+d*e/10))}

spp1411<-function(a,b,c,d,e) {round(1000+(a*b+c+100*d+e))}
spp1412<-function(a,b,c,d,e) {round(1000-(a*b+c+100*d+e))}
spp1413<-function(a,b,c,d,e) {round(1000+(a*b*c+d+100*e))}
spp1414<-function(a,b,c,d,e) {round(1000-(a*b*c+d+100*e))}
spp1415<-function(a,b,c,d,e) {round(1000+(a*b+c+d+100*e)/10)}
spp1416<-function(a,b,c,d,e) {round(1000+(a*b+c+d+100*e)/20)}
spp1417<-function(a,b,c,d,e) {round(1000+(a+b+c+d+100*e)/10)}
spp1418<-function(a,b,c,d,e) {round(1000-(a+b+c+d+100*e)/10)}
spp1419<-function(a,b,c,d,e) {round(1000+(100*a*b*c*d*e)/10)}
spp1420<-function(a,b,c,d,e) {round(1000-(100*a*b*c*d*e)/10)}

spp1421<-function(a,b,c,d,e) {round(8000+(a*b+c+d+e))}
spp1422<-function(a,b,c,d,e) {round(8000-(a*b+c+d+e))}
spp1423<-function(a,b,c,d,e) {round(8000+(a*b*c+d+e))}
spp1424<-function(a,b,c,d,e) {round(8000-(a*b*c+d+e))}
spp1425<-function(a,b,c,d,e) {round(8000+(a*b+c+d+e)/10)}
spp1426<-function(a,b,c,d,e) {round(8000+(a*b+c+d+e)/20)}
spp1427<-function(a,b,c,d,e) {round(8000+(a+b+c+d+e)/10)}
spp1428<-function(a,b,c,d,e) {round(8000-(a+b+c+d+e)/10)}
spp1429<-function(a,b,c,d,e) {round(8000+(a*b*c*d*e)/10)}
spp1430<-function(a,b,c,d,e) {round(8000-(a*b*c*d*e)/10)}

spp1431<-function(a,b,c,d,e) {round(5000+(a*b-c+d+e))}
spp1432<-function(a,b,c,d,e) {round(5000-(a*b-c+d+e))}
spp1433<-function(a,b,c,d,e) {round(5000+(a*b*c-d+e))}
spp1434<-function(a,b,c,d,e) {round(5000-(a*b*c-d+e))}
spp1435<-function(a,b,c,d,e) {round(5000-(a-b-c-d-e)/10)}
spp1436<-function(a,b,c,d,e) {round(5000+(a-b-c-d-e)/20)}
spp1437<-function(a,b,c,d,e) {round(5000+(a-b-c-d-e)/10)}
spp1438<-function(a,b,c,d,e) {round(5000+(a-100*b+c+d+e)/10)}
spp1439<-function(a,b,c,d,e) {round(5000+(a-100*b*c*d*e)/10)}
spp1440<-function(a,b,c,d,e) {round(5000+(a*b*c+d-100*e)/10)}

spp1441<-function(a,b,c,d,e) {round(10000+(a+b+c+d+100*e))}
spp1442<-function(a,b,c,d,e) {round(10000+(a+b+c+d+500*e))}
spp1443<-function(a,b,c,d,e) {round(10000+(a+b+c+d+1000*e))}
spp1444<-function(a,b,c,d,e) {round(10000+(a+b+c+100*d+e))}
spp1445<-function(a,b,c,d,e) {round(10000+(a+b+c+500*d+e))}
spp1446<-function(a,b,c,d,e) {round(10000+(a+b+c+1000*d+e))}
spp1447<-function(a,b,c,d,e) {round(10000+(-100*a+b+c+d+100*e))}
spp1448<-function(a,b,c,d,e) {round(10000+(-500*a+b+c+d+500*e))}
spp1449<-function(a,b,c,d,e) {round(10000+(a+b-100*c+100*d+e))}
spp1450<-function(a,b,c,d,e) {round(10000+(a+b-500*c+500*d+e))}

spp1451<-function(a,b,c,d,e) {1000*a*b*c*d*e}
spp1452<-function(a,b,c,d,e) {6000+(a+b-100*c+500*d+e)}
spp1453<-function(a,b,c,d,e) {6000+(a+b-100*c+1000*d+e)}
spp1454<-function(a,b,c,d,e) {6000+(a+b-500*c+100*d+e)}
spp1455<-function(a,b,c,d,e) {6000+(a+b-1000*c+500*d+e)}
spp1456<-function(a,b,c,d,e) {6000+(a+100*b+c-500*d+e)}
spp1457<-function(a,b,c,d,e) {6000+(a+500*b+c-1000*d+e)}
spp1458<-function(a,b,c,d,e) {6000+(a+1000*b+c-100*d+e)}
spp1459<-function(a,b,c,d,e) {6000+(a+b-500*c+d+1000*e)}
spp1460<-function(a,b,c,d,e) {6000+(a+b-500*c+d+100*e)}

spp1462<-function(a,b,c,d,e) {10000-(a+b-100*c+500*d+e)}
spp1463<-function(a,b,c,d,e) {10000-(a+b-100*c+1000*d+e)}
spp1464<-function(a,b,c,d,e) {10000-(a+b-500*c+100*d+e)}
spp1465<-function(a,b,c,d,e) {10000-(a+b-1000*c+500*d+e)}
spp1466<-function(a,b,c,d,e) {10000-(a+100*b+c-500*d+e)}
spp1467<-function(a,b,c,d,e) {10000-(a+500*b+c-1000*d+e)}
spp1468<-function(a,b,c,d,e) {10000-(a+1000*b+c-100*d+e)}
spp1469<-function(a,b,c,d,e) {10000-(a+b-500*c+d+1000*e)}
spp1470<-function(a,b,c,d,e) {10000-(a+b-500*c+d+100*e)}
spp1461<-function(a,b,c,d,e) {10000-(a+b-500*c*d+100*e)}
################################################
##                                            ##
##  non-linear high complex  response         ##
##                                            ##
################################################

spp1471<-function(a,b,c,d,e) {(0.001*(a-50)^3+3)+(10*b)+(-0.1*(c-50)^2+50)+d+e}
spp1472<-function(a,b,c,d,e) {(0.001*(a-50)^3+3)*(10*b)+(-0.1*(c-50)^2+50)+d+e}
spp1473<-function(a,b,c,d,e) {((1/e)*(a-50)^3+3)+(10*b)+(-0.1*(c-50)^2+50)+d}
spp1474<-function(a,b,c,d,e) {((1/e)*(a-50)^3+3)+(10*b)+(-(1/d)*(c-50)^2+50)}
spp1475<-function(a,b,c,d,e) {((1/e)^2*(a-50)^3+3)+(10*b)+(-(1/d)^2*(c-50)^2+50)}
spp1476<-function(a,b,c,d,e) {((1/d)*(a-50)^3+3)+(10*b)+(-(1/e)*(c-50)^2+50)}
spp1477<-function(a,b,c,d,e) {((1/e)*(a-100)^3+3)*(10*b)+(-(1/d)*(c-50)^2+50)}
spp1478<-function(a,b,c,d,e) {((1/e)^2*(a-100)^3+3)*(10*b)+(-(1/d)^2*(c-100)^2+50)}
spp1479<-function(a,b,c,d,e) {(0.001*(a-50)^3+3)*(10*b)*(1/d)+(-0.1*(c-50)^2+50)+d+e}
spp1480<-function(a,b,c,d,e) {(0.001*(a-50)^3+3)*(10*b)*(-0.001*(e-50)^2+50)+(-0.00000001*(e-50)^5+10)+(-0.1*(c-50)^2+50)+d}

################################################
##                                            ##
##  non-linear high simple   response         ##
##                                            ##
################################################
#complex
spp1481<-function(a,b,c,d,e) {100*a^2+(b^2)/10+c-100*d+e}
spp1482<-function(a,b,c,d,e) {-100*a^2+(b^2)/10+c+100*d+100*e}
spp1483<-function(a,b,c,d,e) {-100*a^2+(b^2)/10+c-100*d+e}
spp1484<-function(a,b,c,d,e) {100*a^2-b+c-100*d^2+e^2}
spp1485<-function(a,b,c,d,e) {(a^2)*b*c*(d^2)*(e^2)}
spp1486<-function(a,b,c,d,e) {-100*a^2+b+c+1000*d^2+e^2}
spp1487<-function(a,b,c,d,e) {100*a^2-b+c-d^2+100*e^2}
spp1488<-function(a,b,c,d,e) {100*e^3+a+(b^2)/10+c-100*d}
spp1489<-function(a,b,c,d,e) {(100*a^2-b+c-100*d^2+e^2)/3}
spp1490<-function(a,b,c,d,e) {(100*a^2-b+c-100*d^2+e^2)/5}


spp1491<-function(a,b,c,d,e) {100*(a+b+c)-(d^e)}
spp1492<-function(a,b,c,d,e) {10000+(a+b+c)-(d*e)}
spp1493<-function(a,b,c,d,e) {100*a*b*c-1/(d/(100*e))}
spp1494<-function(a,b,c,d,e) {1000*e^2-a-b-c-d}
spp1495<-function(a,b,c,d,e) {1000*e^2-a-b-c-d^2}
spp1496<-function(a,b,c,d,e) {(1000*e^2)/a*b*c*d^2}
spp1497<-function(a,b,c,d,e) {(1000*e^2)/a*b*c*d}
spp1498<-function(a,b,c,d,e) {100*a^b^c^d^e}
spp1499<-function(a,b,c,d,e) {1000*a^b+c-100*d^e}
spp1500<-function(a,b,c,d,e) {1000*d^e-100*a^b+c}

spp1501<-function(a,b,c,d,e) {10000-(10000/(100*(a+b+c+d+e)))}
spp1502<-function(a,b,c,d,e) {10000-(10000/(a*b*c*d*e))}
spp1503<-function(a,b,c,d,e) {10000-(10000/(a*b*c*d))+e^2}
spp1504<-function(a,b,c,d,e) {10000*e-(10000/(a*b*c*d))}
spp1505<-function(a,b,c,d,e) {10000-(10000/(100*(a+b+c+d+e)))}
spp1506<-function(a,b,c,d,e) {10000-100*a*b*c*d*e}
spp1507<-function(a,b,c,d,e) {10000-100*(a+b+c+d+e)}
spp1508<-function(a,b,c,d,e) {100*a*b*c-d^(10*e)}
spp1509<-function(a,b,c,d,e) {1000+100*a*b*c-d^(10*e)}
spp1510<-function(a,b,c,d,e) {1000+100*e*b*c-d^(10*a)}


#simple (single factor)

spp1511<-function(a,b,c,d,e) {(rnorm(1, 1000, 100)*(a^2)+0*(a+b+c+d+e))}
spp1512<-function(a,b,c,d,e) {10000-(10000/(100*a))+(0*(a+b+c+d+e))}
spp1513<-function(a,b,c,d,e) {10000-(10000/(100*b))+(0*(a+b+c+d+e))}
spp1514<-function(a,b,c,d,e) {10000-(10000/(100*c))+(0*(a+b+c+d+e))}
spp1515<-function(a,b,c,d,e) {10000-(10000/(100*d))+(0*(a+b+c+d+e))}
spp1516<-function(a,b,c,d,e) {10000-(10000/(100*e))+(0*(a+b+c+d+e))}
spp1517<-function(a,b,c,d,e) {500^a+(0*(a+b+c+d+e))}
spp1518<-function(a,b,c,d,e) {500^b+(0*(a+b+c+d+e))}
spp1519<-function(a,b,c,d,e) {500^c+(0*(a+b+c+d+e))}
spp1520<-function(a,b,c,d,e) {500^d+(0*(a+b+c+d+e))}

spp1521<-function(a,b,c,d,e) {500^e+(0*(a+b+c+d+e))}
spp1522<-function(a,b,c,d,e) {rnorm(1, 1000, 100)*(b^2)+(0*(a+b+c+d+e))}
spp1523<-function(a,b,c,d,e) {rnorm(1, 1000, 100)*(c^2)+(0*(a+b+c+d+e))}
spp1524<-function(a,b,c,d,e) {rnorm(1, 1000, 100)*(d^2)+(0*(a+b+c+d+e))}
spp1525<-function(a,b,c,d,e) {rnorm(1, 1000, 100)*(e^2)+(0*(a+b+c+d+e))}
spp1526<-function(a,b,c,d,e) {(10000/(100*a))+(0*(a+b+c+d+e))}
spp1527<-function(a,b,c,d,e) {(10000/(100*b))+(0*(a+b+c+d+e))}
spp1528<-function(a,b,c,d,e) {(10000/(100*c))+(0*(a+b+c+d+e))}
spp1529<-function(a,b,c,d,e) {(10000/(100*d))+(0*(a+b+c+d+e))}
spp1530<-function(a,b,c,d,e) {(10000/(100*e))+(0*(a+b+c+d+e))}


spp1531<-function(a,b,c,d,e) {100*(a^4)+(0*(a+b+c+d+e))}
spp1532<-function(a,b,c,d,e) {100*(b^4)+(0*(a+b+c+d+e))}
spp1533<-function(a,b,c,d,e) {100*(c^4)+(0*(a+b+c+d+e))}
spp1534<-function(a,b,c,d,e) {100*(d^4)+(0*(a+b+c+d+e))}
spp1535<-function(a,b,c,d,e) {100*(e^4)+(0*(a+b+c+d+e))}
spp1536<-function(a,b,c,d,e) {(10000-(100*a^2))+(0*(a+b+c+d+e))}
spp1537<-function(a,b,c,d,e) {(10000-(100*b^2))+(0*(a+b+c+d+e))}
spp1538<-function(a,b,c,d,e) {(10000-(100*c^2))+(0*(a+b+c+d+e))}
spp1539<-function(a,b,c,d,e) {(10000-(100*d^2))+(0*(a+b+c+d+e))}
spp1540<-function(a,b,c,d,e) {(10000-(100*e^2))+(0*(a+b+c+d+e))}
################################################
##                                            ##
##  non-linear med simple   response          ##
##                                            ##
################################################
#complex
spp1541<-function(a,b,c,d,e) {10*a^2+(b^2)/10+c-10*d+e}
spp1542<-function(a,b,c,d,e) {-10*a^2+(b^2)/10+c+10*d+10*e}
spp1543<-function(a,b,c,d,e) {-10*a^2+(b^2)/10+c-10*d+e}
spp1544<-function(a,b,c,d,e) {10*a^2-b+c-10*d^2+e^2}
spp1545<-function(a,b,c,d,e) {(a^2)*b*c*(d^2)*(e^2)}
spp1546<-function(a,b,c,d,e) {-10*a^2+b+c+10*d^2+e^2}
spp1547<-function(a,b,c,d,e) {a^2-b+c-d^2+10*e^2}
spp1548<-function(a,b,c,d,e) {10*e^3+a+(b^2)/10+c-10*d}
spp1549<-function(a,b,c,d,e) {(10*a^2-b+c-10*d^2+e^2)/3}
spp1550<-function(a,b,c,d,e) {(10*a^2-b+c-10*d^2+e^2)/5}


spp1551<-function(a,b,c,d,e) {10*(a+b+c)-(d^e)}
spp1552<-function(a,b,c,d,e) {10+(a+b+c)-(d*e)}
spp1553<-function(a,b,c,d,e) {10*a*b*c-1/(d/(10*e))}
spp1554<-function(a,b,c,d,e) {10*e^2-a-b-c-d}
spp1555<-function(a,b,c,d,e) {100*e^2-a-b-c-d^2}
spp1556<-function(a,b,c,d,e) {(100*e^2)/a*b*c*d^2}
spp1557<-function(a,b,c,d,e) {(100*e^2)/a*b*c*d}
spp1558<-function(a,b,c,d,e) {10*a^b^c^d^e}
spp1559<-function(a,b,c,d,e) {100*a^b+c-10*d^e}
spp1560<-function(a,b,c,d,e) {100*d^e-10*a^b+c}

spp1561<-function(a,b,c,d,e) {1000-(1000/(100*(a+b+c+d+e)))}
spp1562<-function(a,b,c,d,e) {1000-(1000/(a*b*c*d*e))}
spp1563<-function(a,b,c,d,e) {1000-(1000/(a*b*c*d))+e^2}
spp1564<-function(a,b,c,d,e) {1000*e-(100/(a*b*c*d))}
spp1565<-function(a,b,c,d,e) {100-(100/(100*(a+b+c+d+e)))}
spp1566<-function(a,b,c,d,e) {100-10*a*b*c*d*e}
spp1567<-function(a,b,c,d,e) {100-10*(a+b+c+d+e)}
spp1568<-function(a,b,c,d,e) {10*a*b*c-d^(e)}
spp1569<-function(a,b,c,d,e) {100+10*a*b*c-d^(10*e)}
spp1570<-function(a,b,c,d,e) {100+10*e*b*c-d^(10*a)}


#simple (single factor)

spp1571<-function(a,b,c,d,e) {(rnorm(1, 100, 10)*(a^2)+0*(a+b+c+d+e))}
spp1572<-function(a,b,c,d,e) {1000-(1000/(10*a))+(0*(a+b+c+d+e))}
spp1573<-function(a,b,c,d,e) {1000-(1000/(10*b))+(0*(a+b+c+d+e))}
spp1574<-function(a,b,c,d,e) {1000-(1000/(10*c))+(0*(a+b+c+d+e))}
spp1575<-function(a,b,c,d,e) {1000-(1000/(10*d))+(0*(a+b+c+d+e))}
spp1576<-function(a,b,c,d,e) {1000-(1000/(10*e))+(0*(a+b+c+d+e))}
spp1577<-function(a,b,c,d,e) {50^a+(0*(a+b+c+d+e))}
spp1578<-function(a,b,c,d,e) {50^b+(0*(a+b+c+d+e))}
spp1579<-function(a,b,c,d,e) {50^c+(0*(a+b+c+d+e))}
spp1580<-function(a,b,c,d,e) {50^d+(0*(a+b+c+d+e))}

spp1581<-function(a,b,c,d,e) {50^e+(0*(a+b+c+d+e))}
spp1582<-function(a,b,c,d,e) {rnorm(1, 100, 100)*(b^2)+(0*(a+b+c+d+e))}
spp1583<-function(a,b,c,d,e) {rnorm(1, 100, 100)*(c^2)+(0*(a+b+c+d+e))}
spp1584<-function(a,b,c,d,e) {rnorm(1, 100, 100)*(d^2)+(0*(a+b+c+d+e))}
spp1585<-function(a,b,c,d,e) {rnorm(1, 100, 100)*(e^2)+(0*(a+b+c+d+e))}
spp1586<-function(a,b,c,d,e) {(1000/(10*a))+(0*(a+b+c+d+e))}
spp1587<-function(a,b,c,d,e) {(1000/(10*b))+(0*(a+b+c+d+e))}
spp1588<-function(a,b,c,d,e) {(1000/(10*c))+(0*(a+b+c+d+e))}
spp1589<-function(a,b,c,d,e) {(1000/(10*d))+(0*(a+b+c+d+e))}
spp1590<-function(a,b,c,d,e) {(1000/(10*e))+(0*(a+b+c+d+e))}


spp1591<-function(a,b,c,d,e) {10*(a^4)+(0*(a+b+c+d+e))}
spp1592<-function(a,b,c,d,e) {10*(b^4)+(0*(a+b+c+d+e))}
spp1593<-function(a,b,c,d,e) {10*(c^4)+(0*(a+b+c+d+e))}
spp1594<-function(a,b,c,d,e) {10*(d^4)+(0*(a+b+c+d+e))}
spp1595<-function(a,b,c,d,e) {10*(e^4)+(0*(a+b+c+d+e))}
spp1596<-function(a,b,c,d,e) {(1000-(10*a^2))+(0*(a+b+c+d+e))}
spp1597<-function(a,b,c,d,e) {(1000-(10*b^2))+(0*(a+b+c+d+e))}
spp1598<-function(a,b,c,d,e) {(1000-(10*c^2))+(0*(a+b+c+d+e))}
spp1599<-function(a,b,c,d,e) {(1000-(10*d^2))+(0*(a+b+c+d+e))}
spp1600<-function(a,b,c,d,e) {(1000-(10*e^2))+(0*(a+b+c+d+e))}

################################################
##                                            ##
##  non-linear med low response               ##
##                                            ##
################################################
#complex
spp1601<-function(a,b,c,d,e) {a^2+(b^2)/10+c-10*d+e}
spp1602<-function(a,b,c,d,e) {-a^2+(b^2)/10+c+10*d+10*e}
spp1603<-function(a,b,c,d,e) {-a^2+(b^2)/10+c-10*d+e}
spp1604<-function(a,b,c,d,e) {a^2-b+c-10*d^2+e^2}
spp1605<-function(a,b,c,d,e) {(a^2)*b*c*(d^2)*(e^2)}
spp1606<-function(a,b,c,d,e) {a^2+b+c+d^2+e^2}
spp1607<-function(a,b,c,d,e) {a^2-b+c-d^2+e^2}
spp1608<-function(a,b,c,d,e) {e^3+a+(b^2)/10+c-d}
spp1609<-function(a,b,c,d,e) {(a^2-b+c-d^2+e^2)/3}
spp1610<-function(a,b,c,d,e) {(a^2-b+c-d^2+e^2)/5}


spp1611<-function(a,b,c,d,e) {(a+b+c)-(d^e)}
spp1612<-function(a,b,c,d,e) {1+(a+b+c)-(d*e)}
spp1613<-function(a,b,c,d,e) {1*a*b*c-1/(d/(10*e))}
spp1614<-function(a,b,c,d,e) {e^2-a-b-c-d}
spp1615<-function(a,b,c,d,e) {e^2-a-b-c-d^2}
spp1616<-function(a,b,c,d,e) {(10*e^2)/a*b*c*d^2}
spp1617<-function(a,b,c,d,e) {(10*e^2)/a*b*c*d}
spp1618<-function(a,b,c,d,e) {a^b^c^d^e}
spp1619<-function(a,b,c,d,e) {a^b+c-10*d^e}
spp1620<-function(a,b,c,d,e) {d^e-10*a^b+c}

spp1621<-function(a,b,c,d,e) {100-(100/(10*(a+b+c+d+e)))}
spp1622<-function(a,b,c,d,e) {100-(100/(a*b*c*d*e))}
spp1623<-function(a,b,c,d,e) {100-(100/(a*b*c*d))+e^2}
spp1624<-function(a,b,c,d,e) {100*e-(10/(a*b*c*d))}
spp1625<-function(a,b,c,d,e) {10-(10/(10*(a+b+c+d+e)))}
spp1626<-function(a,b,c,d,e) {10-a*b*c*d*e}
spp1627<-function(a,b,c,d,e) {10-(a+b+c+d+e)}
spp1628<-function(a,b,c,d,e) {a*b*c-d^(e)}
spp1629<-function(a,b,c,d,e) {10+a*b*c-d^(e)}
spp1630<-function(a,b,c,d,e) {10+e*b*c-d^(a)}


#simple (single factor)

spp1631<-function(a,b,c,d,e) {(rnorm(1, 10, 1)*(a^2)+0*(a+b+c+d+e))}
spp1632<-function(a,b,c,d,e) {100-(100/(10*a))+(0*(a+b+c+d+e))}
spp1633<-function(a,b,c,d,e) {100-(100/(10*b))+(0*(a+b+c+d+e))}
spp1634<-function(a,b,c,d,e) {100-(100/(10*c))+(0*(a+b+c+d+e))}
spp1635<-function(a,b,c,d,e) {100-(100/(10*d))+(0*(a+b+c+d+e))}
spp1636<-function(a,b,c,d,e) {100-(100/(10*e))+(0*(a+b+c+d+e))}
spp1637<-function(a,b,c,d,e) {5^a+(0*(a+b+c+d+e))}
spp1638<-function(a,b,c,d,e) {5^b+(0*(a+b+c+d+e))}
spp1639<-function(a,b,c,d,e) {5^c+(0*(a+b+c+d+e))}
spp1640<-function(a,b,c,d,e) {5^d+(0*(a+b+c+d+e))}

spp1641<-function(a,b,c,d,e) {5^e+(0*(a+b+c+d+e))}

spp1642<-function(a,b,c,d,e){rnorm(1, 10, 10)*(b^2)+(0*(a+b+c+d+e))}
spp1643<-function(a,b,c,d,e) {rnorm(1, 10, 10)*(c^2)+(0*(a+b+c+d+e))}
spp1644<-function(a,b,c,d,e) {rnorm(1, 10, 10)*(d^2)+(0*(a+b+c+d+e))}
spp1645<-function(a,b,c,d,e) {rnorm(1, 10, 10)*(e^2)+(0*(a+b+c+d+e))}
spp1646<-function(a,b,c,d,e) {(100/(10*a))+(0*(a+b+c+d+e))}
spp1647<-function(a,b,c,d,e) {(100/(10*b))+(0*(a+b+c+d+e))}
spp1648<-function(a,b,c,d,e) {(100/(10*c))+(0*(a+b+c+d+e))}
spp1649<-function(a,b,c,d,e) {(100/(10*d))+(0*(a+b+c+d+e))}
spp1650<-function(a,b,c,d,e) {(100/(10*e))+(0*(a+b+c+d+e))}


spp1651<-function(a,b,c,d,e) {(a^4)+(0*(a+b+c+d+e))}
spp1652<-function(a,b,c,d,e) {(b^4)+(0*(a+b+c+d+e))}
spp1653<-function(a,b,c,d,e) {(c^4)+(0*(a+b+c+d+e))}
spp1654<-function(a,b,c,d,e) {(d^4)+(0*(a+b+c+d+e))}
spp1655<-function(a,b,c,d,e) {(e^4)+(0*(a+b+c+d+e))}
spp1656<-function(a,b,c,d,e) {(100-(a^2))+(0*(a+b+c+d+e))}
spp1657<-function(a,b,c,d,e) {(100-(b^2))+(0*(a+b+c+d+e))}
spp1658<-function(a,b,c,d,e) {(100-(c^2))+(0*(a+b+c+d+e))}
spp1659<-function(a,b,c,d,e) {(100-(d^2))+(0*(a+b+c+d+e))}
spp1660<-function(a,b,c,d,e) {(100-(e^2))+(0*(a+b+c+d+e))}


################################################
##                                            ##
##  Random responses                          ##
##                                            ##
################################################

#low
spp1661<-function(a,b,c,d,e) {rnorm(1, 50, 50)+(0*(a+b+c+d+e))} #generate 1 number with mean=50, stdev=50
spp1662<-function(a,b,c,d,e) {rnorm(1, 50, 40)+(0*(a+b+c+d+e))}
spp1663<-function(a,b,c,d,e) {rnorm(1, 50, 20)+(0*(a+b+c+d+e))}
spp1664<-function(a,b,c,d,e) {rnorm(1, 50, 10)+(0*(a+b+c+d+e))}
spp1665<-function(a,b,c,d,e) {rnorm(1, 50, 100)+(0*(a+b+c+d+e))}
spp1666<-function(a,b,c,d,e) {rnorm(1, 50, 200)+(0*(a+b+c+d+e))}
spp1667<-function(a,b,c,d,e) {rnorm(1, 50, 1)+(0*(a+b+c+d+e))}
spp1668<-function(a,b,c,d,e) {rnorm(1, 5, 100)+(0*(a+b+c+d+e))}
spp1669<-function(a,b,c,d,e) {rnorm(1, 5, 50)+(0*(a+b+c+d+e))}
spp1670<-function(a,b,c,d,e) {rnorm(1, 5, 5)+(0*(a+b+c+d+e))}

spp1671<-function(a,b,c,d,e) {rnorm(1, 1, 50)+(0*(a+b+c+d+e))}
spp1672<-function(a,b,c,d,e) {rnorm(1, 2, 50)+(0*(a+b+c+d+e))}
spp1673<-function(a,b,c,d,e) {rnorm(1, 10, 50)+(0*(a+b+c+d+e))}
spp1674<-function(a,b,c,d,e) {rnorm(1, 20, 50)+(0*(a+b+c+d+e))}
spp1675<-function(a,b,c,d,e) {rnorm(1, 30, 50)+(0*(a+b+c+d+e))}
spp1676<-function(a,b,c,d,e) {rnorm(1, 40, 50)+(0*(a+b+c+d+e))}
spp1677<-function(a,b,c,d,e) {rnorm(1, 60, 50)+(0*(a+b+c+d+e))}
spp1678<-function(a,b,c,d,e) {rnorm(1, 70, 50)+(0*(a+b+c+d+e))}
spp1679<-function(a,b,c,d,e) {rnorm(1, 80, 50)+(0*(a+b+c+d+e))}
spp1680<-function(a,b,c,d,e) {rnorm(1, 90, 50)+(0*(a+b+c+d+e))}
spp1681<-function(a,b,c,d,e) {rnorm(1, 100, 50)+(0*(a+b+c+d+e))}
spp1682<-function(a,b,c,d,e) {rnorm(1, 1, 5)+(0*(a+b+c+d+e))}
spp1683<-function(a,b,c,d,e) {rnorm(1, 1, 50)+(0*(a+b+c+d+e))}
spp1684<-function(a,b,c,d,e) {rnorm(1, 1, 100)+(0*(a+b+c+d+e))}

#medium


spp1725<-function(a,b,c,d,e) {rnorm(1, 500, 500)+(0*(a+b+c+d+e))} #generate 1 number with mean=50, stdev=50
spp1726<-function(a,b,c,d,e) {rnorm(1, 500, 400)+(0*(a+b+c+d+e))}
spp1727<-function(a,b,c,d,e) {rnorm(1, 500, 200)+(0*(a+b+c+d+e))}
spp1728<-function(a,b,c,d,e) {rnorm(1, 500, 100)+(0*(a+b+c+d+e))}
spp1685<-function(a,b,c,d,e) {rnorm(1, 500, 1000)+(0*(a+b+c+d+e))}
spp1686<-function(a,b,c,d,e) {rnorm(1, 500, 2000)+(0*(a+b+c+d+e))}
spp1687<-function(a,b,c,d,e) {rnorm(1, 500, 10)+(0*(a+b+c+d+e))}
spp1688<-function(a,b,c,d,e) {rnorm(1, 50, 1000)+(0*(a+b+c+d+e))}
spp1689<-function(a,b,c,d,e) {rnorm(1, 50, 500)+(0*(a+b+c+d+e))}
spp1690<-function(a,b,c,d,e) {rnorm(1, 50, 50)+(0*(a+b+c+d+e))}

spp1691<-function(a,b,c,d,e) {rnorm(1, 10, 500)+(0*(a+b+c+d+e))}
spp1692<-function(a,b,c,d,e) {rnorm(1, 20, 500)+(0*(a+b+c+d+e))}
spp1693<-function(a,b,c,d,e) {rnorm(1, 100, 500)+(0*(a+b+c+d+e))}
spp1694<-function(a,b,c,d,e) {rnorm(1, 200, 500)+(0*(a+b+c+d+e))}
spp1695<-function(a,b,c,d,e) {rnorm(1, 300, 500)+(0*(a+b+c+d+e))}
spp1696<-function(a,b,c,d,e) {rnorm(1, 400, 500)+(0*(a+b+c+d+e))}
spp1697<-function(a,b,c,d,e) {rnorm(1, 600, 500)+(0*(a+b+c+d+e))}
spp1698<-function(a,b,c,d,e) {rnorm(1, 700, 500)+(0*(a+b+c+d+e))}
spp1699<-function(a,b,c,d,e) {rnorm(1, 800, 500)+(0*(a+b+c+d+e))}
spp1700<-function(a,b,c,d,e) {rnorm(1, 900, 500)+(0*(a+b+c+d+e))}
spp1701<-function(a,b,c,d,e) {rnorm(1, 1000, 500)+(0*(a+b+c+d+e))}
spp1702<-function(a,b,c,d,e) {rnorm(1, 10, 50)+(0*(a+b+c+d+e))}
spp1703<-function(a,b,c,d,e) {rnorm(1, 10, 500)+(0*(a+b+c+d+e))}
spp1704<-function(a,b,c,d,e) {rnorm(1, 10, 1000)+(0*(a+b+c+d+e))}

#high
spp1729<-function(a,b,c,d,e) {rnorm(1, 5000, 5000)+(0*(a+b+c+d+e))} #generate 1 number with mean=50, stdev=50
spp1730<-function(a,b,c,d,e) {rnorm(1, 5000, 4000)+(0*(a+b+c+d+e))}
spp1731<-function(a,b,c,d,e) {rnorm(1, 5000, 2000)+(0*(a+b+c+d+e))}
spp1732<-function(a,b,c,d,e) {rnorm(1, 5000, 1000)+(0*(a+b+c+d+e))}
spp1705<-function(a,b,c,d,e) {rnorm(1, 5000, 10000)+(0*(a+b+c+d+e))}
spp1706<-function(a,b,c,d,e) {rnorm(1, 5000, 20000)+(0*(a+b+c+d+e))}
spp1707<-function(a,b,c,d,e) {rnorm(1, 5000, 100)+(0*(a+b+c+d+e))}
spp1708<-function(a,b,c,d,e) {rnorm(1, 500, 10000)+(0*(a+b+c+d+e))}
spp1709<-function(a,b,c,d,e) {rnorm(1, 500, 5000)+(0*(a+b+c+d+e))}
spp1710<-function(a,b,c,d,e) {rnorm(1, 500, 500)+(0*(a+b+c+d+e))}

spp1711<-function(a,b,c,d,e) {rnorm(1, 100, 5000)+(0*(a+b+c+d+e))}
spp1712<-function(a,b,c,d,e) {rnorm(1, 200, 5000)+(0*(a+b+c+d+e))}
spp1713<-function(a,b,c,d,e) {rnorm(1, 1000, 5000)+(0*(a+b+c+d+e))}
spp1714<-function(a,b,c,d,e) {rnorm(1, 2000, 5000)+(0*(a+b+c+d+e))}
spp1715<-function(a,b,c,d,e) {rnorm(1, 3000, 5000)+(0*(a+b+c+d+e))}
spp1716<-function(a,b,c,d,e) {rnorm(1, 4000, 5000)+(0*(a+b+c+d+e))}
spp1717<-function(a,b,c,d,e) {rnorm(1, 6000, 5000)+(0*(a+b+c+d+e))}
spp1718<-function(a,b,c,d,e) {rnorm(1, 7000, 5000)+(0*(a+b+c+d+e))}
spp1719<-function(a,b,c,d,e) {rnorm(1, 8000, 5000)+(0*(a+b+c+d+e))}
spp1720<-function(a,b,c,d,e) {rnorm(1, 9000, 5000)+(0*(a+b+c+d+e))}
spp1721<-function(a,b,c,d,e) {rnorm(1, 10000, 5000)+(0*(a+b+c+d+e))}
spp1722<-function(a,b,c,d,e) {rnorm(1, 100, 500)+(0*(a+b+c+d+e))}
spp1723<-function(a,b,c,d,e) {rnorm(1, 100, 5000)+(0*(a+b+c+d+e))}
spp1724<-function(a,b,c,d,e) {rnorm(1, 100, 10000)+(0*(a+b+c+d+e))}


rrarefy2<-function(){
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function is meaningful only for integers (counts)")
  x <- as.matrix(x)
  if (ncol(x) == 1) 
    x <- t(x)
  if (length(sample) > 1 && length(sample) != nrow(x)) 
    stop(gettextf("length of 'sample' and number of rows of 'x' do not match"))
  sample <- rep(sample, length = nrow(x))
  colnames(x) <- colnames(x, do.NULL = FALSE)
  nm <- colnames(x)
  if (any(rowSums(x) < sample)) 
    warning("Some row sums < 'sample' and are not rarefied")
  for (i in 1:nrow(x)) {
    if (sum(x[i, ]) <= sample[i]) 
      next 
    row <- sample(rep(nm, times = x[i, ]), sample[i], replace = T)
    row <- table(row)
    ind <- names(row)
    x[i, ] <- 0
    x[i, ind] <- row
    }
  }
  x
}