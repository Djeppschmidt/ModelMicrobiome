#limaBenchmark


deseqF<-funciton(){

}
limmaF<-funciton(){
  counts<-t(as.data.frame(as.matrix(otu_table(x))))
  categories<-data.frame("Factor"=as.factor(sample_data(x)$Factor))

  f.design<-model.matrix(~factor, categories)
  #colnames(design) <- gsub("group3", "", colnames(design))
  contr.matrix <- makeContrasts(
    C1 = Factor2-Factor1,
    C2 = Factor3-Factor1,
    C3 = Factor4-Factor1,
    C4 = Factor5-Factor1,
    C5 = Factor6-Factor1,
    C6 = Factor7-Factor1,
    C7 = Factor8-Factor1,
    C8 = Factor9-Factor1,
    C9 = Factor10-Factor1,
    C10 = Factor11-Factor1,
    C11 = Factor12-Factor1,

    levels = colnames(f.design))
  contr.matrix

  dge <- DGEList(counts=counts)
  keep <- filterByExpr(dge, f.design)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge) #what happens if we don't do this step?
  v <- voom(dge, f.design, plot=TRUE)
  fitV <- lmFit(v, f.design)
  fitV <- contrasts.fit(fitV, contrasts=contr.matrix)
  fitV <- eBayes(fitV, trend=TRUE)
  #topgenes<-topTable(fitV, number=10, coef=1)
  #siggenes<-Muri.topgenes[topgenes$adj.P.Val<0.05,]

  out<-summary(decideTests(fitV))
  out
}

vst.limma<-funciton(ps){
  counts<-t(as.data.frame(as.matrix(otu_table(ps))))
  categories<-data.frame("Factor"=as.factor(sample_data(ps)$Factor))
  dge <- DGEList(counts=counts)
  dge <- calcNormFactors(dge)
  ps<-otu_table(dge)
  ps#make sure to export otu table
}
vst.deseq2<-funciton(ps){

  vst...
  ps<-otu_table(vst)
}

PERMANOVA<-function(x){
  tab<-as.data.frame(as.matrix(otu_table(x)))
  metadata<-data.frame(sample_data(x))
  out <- adonis(tab~factor+f1+F2+F3+F4+F5, data=metadata)
  out #
}

score.perm<-function(){

  env.comm<-readRDS("")
}

score.da<-function(list){

  t<-list("C1"=,"C2"=,"C3"=,"C4"=,"C5"=,"C6"=,"C7"=,"C8"=,"C9"=,"C10"=,"C11"=)#Make a list of taxa names that are true increases and decreases
  list-t # work n this question later (how to compare outputs)
}


#benchmark preprocessing/take a list of input otu tables, and a function/pipeline argument
#outputs for a single function
#use lapply for multiple functions
Benchmark.alg<-function(x, arguments, pipe){
  require(phyloseq)
  require(limma)
  require(deseq2)
  require()...# finish this list!!

  #import permanova community
  env.comm<-readRDS("")
  #import differential abundance community
  da.comm<-readRDS("")
  #benchmark environmental correlations
  #make a list of normalization to benchmark
  bench.list<-list(vist.deseq2(env.comm), ...)
  perm.out<-lapply(bench.list, PERMANOVA)
  out[,1]<-score.perm(benc.list) #score the permanova r^2 values
  out[,2]<-score.da() #score the differential abundance counts
  out<-data.frame()
  for(i in bench.list){
    out[i,1]<-score.da()
    out[i,2]<-score.perm()
  }#
  out #an array with each
}

#functions for extracting plots from an array!
plot_da<-function(){} #colored line graph for different algorithms

plot_perm<-function(){} #heatmap env.conditions*algorithm

plot_da.pp<-function(){} #colored line graph for different preprocessing

plot_perm.pp<-function(){} #heatmap env.conditions*preprocessing (sum of residuals?)
