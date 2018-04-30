### limma-voom ###
#design<-model.matrix(~categories)

limmaVoomPipe<-function(dat, design, contrasts){
OTU<-as.matrix(t(otu_table(dat)))
sData<-sample_data(dat)
attach(sData)
dge<-DGEList(counts=OTU)
dge<-calcNormFactors(dge) #calculate the normalization factors via EdgeR (see EdgeR documentation)

V <- voom(dge, design, plot=TRUE)
Fit<-lmFit(V, design)
Contr<-makeContrasts(categories , levels=colnames(coef(Fit)))
Tmp<-contrasts.fit(Fit,Contr)
Tmp<-eBayes(Tmp)
Tmp2<-topTable(Tmp, coef=1, sort.by="P", n="INF")#add column for bayesian support
Tmp2

}
