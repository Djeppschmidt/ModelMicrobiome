
#' shell script for benchmarking
#' @param reps number of replicate communities
#' @param commonN number of common species
#' @param groupN number of unique taxa to groups
#' @param singleN number of unique taxa to samples
#' @param D average sampling depth
#' @param V variation in sampling depth
#' @keywords benchmark
#' @export
#' @examples test.script<-benchmark.MM(reps=3, commonN=20, groupN=10, singleN=10, D=2000, V=1000)
#'benchmark.MM()

benchmark.MM<-function(reps, commonN, groupN, singleN, D, V){
  require(phyloseq)
  #require(ModelMicrobiome)
  require(limma)
  require(DESeq2)
  require(vegan)
  array<-c(rep(commonN, reps))
  names(array)<-c(paste0("rep", 1:reps, sep=""))
  array<-sapply(array, run.analysis, simplify=F, USE.NAMES = T, groupN, singleN, D, V)
  array
}

#' workhorse function for benchmark.MM
#' @param commonN number of common species
#' @param groupN number of unique taxa to groups
#' @param singleN number of unique taxa to samples
#' @param D average sampling depth
#' @param V variation in sampling depth
#' @keywords benchmark
#' @export
#' @examples
#' run.analysis()
run.analysis<-function(commonN, groupN, singleN, D, V){
    AllSpp<-c(paste0("spp", c(1:700), sep="")) # make a quick list of all species functions
    AllSpp<-lapply(AllSpp, get) # connect function to name
    AllSpp<-unlist(AllSpp)  # format to be read by downstream functions
    names(AllSpp)<-c(paste0("spp", c(1:700)))

    # Define list of 5 species w/ global distribution
    global.spp<-names(sample(AllSpp, commonN, replace=F))

# define list of species w/ regional distribution
    group.spp<-NULL
    group.spp$group1<-names(sample(AllSpp, groupN, replace=F))
    group.spp$group2<-names(sample(AllSpp, groupN, replace=F))
    group.spp$group3<-names(sample(AllSpp, groupN, replace=F))
    group.spp$group4<-names(sample(AllSpp, groupN, replace=F))
    group.spp$group5<-names(sample(AllSpp, groupN, replace=F))
    group.spp$group6<-names(sample(AllSpp, groupN, replace=F))

# define list of species found at each site
    rando.spp<-NULL
    rando.spp$Site1<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
    rando.spp$Site2<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
    rando.spp$Site3<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
    rando.spp$Site4<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
    rando.spp$Site5<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
    rando.spp$Site6<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
    rando.spp$Site7<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
    rando.spp$Site8<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
    rando.spp$Site9<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
    rando.spp$Site10<-unique(c(names(sample(AllSpp,singleN, replace=F)), c(group.spp$group2), global.spp))
    rando.spp$Site11<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
    rando.spp$Site12<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
    rando.spp$Site13<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
    rando.spp$Site14<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
    rando.spp$Site15<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
    rando.spp$Site16<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
    rando.spp$Site17<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
    rando.spp$Site18<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
    rando.spp$Site19<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
    rando.spp$Site20<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
    rando.spp$Site21<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
    rando.spp$Site22<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
    rando.spp$Site23<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
    rando.spp$Site24<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
    rando.spp$Site25<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
    rando.spp$Site26<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
    rando.spp$Site27<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
    rando.spp$Site28<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
    rando.spp$Site29<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
    rando.spp$Site30<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))

# make list of unique species arrays

    library(reshape2)
    f1c1<-c(5,5,5,5,5,5) # number of selections
    f1c2<-c(1,3,10,30,60,15) # mean value of selections
    f1c3<-c(0.5,1,4,10,20,5) # SD of selections
    F1.frame<-mapply(rnorm, f1c1,f1c2,f1c3) # pick Factor 1 value for each site
    F1<-reshape2::melt(F1.frame)

#F2
    f2c1<-c(5,5,5,5,5,5) # number of selections
    f2c2<-c(34,30,10,55,35,60) # mean value of selections
    f2c3<-c(10,10,3,10,1,20) # SD of selections
    F2.frame<-mapply(rnorm, f2c1,f2c2,f2c3) # pick Factor 2 value for each site
    F2<-reshape2::melt(F2.frame)

#F3
    f3c1<-c(5,5,5,5,5,5) # number of selections
    f3c2<-c(1,3,10,15,3,15) # mean value of selections
    f3c3<-c(0.5,1,3,3,1,5) # SD of selections
    F3.frame<-mapply(rnorm, f3c1,f3c2,f3c3) # pick Factor 3 value for each site
    F3<-reshape2::melt(F3.frame)

#F4
    f4c1<-c(5,5,5,5,5,5) # number of selections
    f4c2<-c(1,3,10,15,3,15) # mean value of selections
    f4c3<-c(0.5,0.5,3,3,1,5) # SD of selections
    F4.frame<-mapply(rnorm, f4c1,f4c2,f4c3) # pick Factor 4 value for each site
    F4<-reshape2::melt(F4.frame)

#F5
    f5c1<-c(5,5,5,5,5,5) # number of selections
    f5c2<-c(50,40,30,20,10,15) # mean value of selections
    f5c3<-c(20,10,10,5,3,5) # SD of selections
    F5.frame<-mapply(rnorm, f5c1,f5c2,f5c3) # pick Factor 5 value for each site
    F5<-reshape2::melt(F5.frame)
    Factors<-data.frame(F1$value,F2$value,F3$value,F4$value,F5$value) # combine factors into data table
    Sites<-c(paste0("Site", 1:30))
    rownames(Factors)<-Sites
    colnames(Factors)<-c("F1","F2","F3","F4","F5")

    output<-list("model"=NULL, "spplist"=NULL, "raw"=NULL, "eRare"=NULL, "pRare"=NULL, "scaled"=NULL, "deseqVST"=NULL, "limma"=NULL)

    #output$model<-NULL
    output$spplist<-rando.spp
    output$model$comm<-make.refcomm(rando.spp, Factors) # output a phyloseq object... will make a list of phyloseq objects
    sample_data(output$model$comm)$Density<-sample_sums(output$model$comm)# add sample sums
    sample_data(output$model$comm)$DensityF<-sample_sums(output$model$comm)/mean(sample_sums(output$model$comm))
    sample_data(output$model$comm)$Factor<-as.factor(c(rep("one",5),rep("two",5),rep("three",5),rep("four",5),rep("five",5),rep("six",5)))
    sample_data(output$model$comm)$Factor2<-as.factor(c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5)))

    sample<-set.seqDepth(D,V)
    output$raw$comm<-make.rarefy2(output$model$comm, sample)         # do normalizations!!

    output$eRare$comm<-make.rarefy2(output$raw$comm, min(sample_sums(output$raw$comm)))
    output$pRare$comm<-make.rarefy2(output$raw$comm, D*sample_sums(output$raw$comm)/mean(sample_sums(output$raw$comm)))
    output$scaled$comm<-make.scaled2(output$raw$comm, val=D*sample_sums(output$model$comm)/mean(sample_sums(output$model$comm)), scale=D)
    output$deseqVST$comm<-make.deseqVST(output$raw$comm, "Factor", l=1)
    output$limma$comm<-make.limmaVST(output$raw$comm, "Factor")
    print(output$limma$comm)
  # add in any additional methods we want to add

    output$raw$PERMANOVA<-make.PERMANOVA(output$raw$comm)   # make permanova tables
    output$eRare$PERMANOVA<-make.PERMANOVA(output$eRare$comm)
    output$pRare$PERMANOVA<-make.PERMANOVA(output$pRare$comm)
    output$scaled$PERMANOVA<-make.PERMANOVA(output$scaled$comm)
    output$deseqVST$PERMANOVA<-make.PERMANOVA(output$deseqVST$comm)
    output$limma$PERMANOVA<-make.PERMANOVA(output$limma$comm)
  # add in any additional methods we want to add

    output$raw$LII<-LII(output$model$comm, output$raw$comm)                    # show amount of information retained after corrections
    output$eRare$LII<-LII(output$model$comm, output$eRare$comm)
    output$pRare$LII<-LII(output$model$comm, output$pRare$comm)
    output$scaled$LII<-LII(output$model$comm, output$scaled$comm)
    output$deseqVST$LII<-LII(output$model$comm, output$deseqVST$comm)
    output$limma$LII<-LII(output$model$comm, output$limma$comm)
  # add in any additional methods we want to add

    output$raw$SVI<-SVI(output$model$comm, output$raw$comm)                  # show over-fitting of model
    output$eRare$SVI<-SVI(output$model$comm, output$eRare$comm)
    output$pRare$SVI<-SVI(output$model$comm, output$pRare$comm)
    output$scaled$SVI<-SVI(output$model$comm, output$scaled$comm)
    output$deseqVST$SVI<-SVI(output$model$comm, output$deseqVST$comm)
    output$limma$SVI<-SVI(output$model$comm, output$limma$comm)

    output
}


#' construct reference community
#' @param rando.spp list of species lists
#' @param Factors data frame of "environmental" data
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' make.refcomm()
make.refcomm<-function(rando.spp, Factors){
l1<-NULL
for (i in 1:length(rando.spp[[1]])){
  l1[i]<-do.call(rando.spp[[1]][i], list(Factors[1,1],Factors[1,2],Factors[1,3],Factors[1,4],Factors[1,5]))
  }
#l1<-data.frame("Site1"=l1, "Spp"=rando.spp[[1]])
names(l1)<-rando.spp[[1]]
for (r in 2:nrow(Factors)) # for each site...
{ l2<-NULL
  for (i in 1:length(rando.spp[[r]])){  # for each species in site...
    l2[i]<-do.call(rando.spp[[r]][i], list(Factors[r,1],Factors[r,2],Factors[r,3],Factors[r,4],Factors[r,5]))
    }
  names(l2)<-rando.spp[[r]]
  l1<-merge(as.data.frame(l1),as.data.frame(l2), by=0, all=T)
  rownames(l1)<-l1$Row.names
  colnames(l1)[colnames(l1) == "l1"] <- "Site1"
 colnames(l1)[colnames(l1) == "l2"] <- paste("Site", r, sep="")
 l1<-l1[,-1]
  }
l1<-round(l1)
l1[mapply(is.infinite, l1)] <- NA
l1[is.na(l1)]<-0
l1[l1<0]<-0
otu<-otu_table(l1, taxa_are_rows = T)
Sa<-sample_data(Factors)
out<-phyloseq(otu, Sa)
out
}

#' construct community
#' @param rando.spp list of species lists
#' @param Factors data frame of "environmental" data
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.comm2()
make.comm2<-function(rando.spp, Factors){
l1<-NULL
for (i in 1:length(rando.spp[[1]])){
  l1[i]<-do.call(rando.spp[[1]][i], list(Factors[1,1],Factors[1,2],Factors[1,3],Factors[1,4],Factors[1,5]))
  }
#l1<-data.frame("Site1"=l1, "Spp"=rando.spp[[1]])
names(l1)<-rando.spp[[1]]
for (r in 2:nrow(Factors)) # for each site...
{ l2<-NULL
  for (i in 1:length(rando.spp[[r]])){  # for each species in site...
    l2[i]<-do.call(rando.spp[[r]][i], list(Factors[r,1],Factors[r,2],Factors[r,3],Factors[r,4],Factors[r,5]))
    }
  names(l2)<-rando.spp[[r]]
  l1<-merge(as.data.frame(l1),as.data.frame(l2), by=0, all=T)
  rownames(l1)<-l1$Row.names
  colnames(l1)[colnames(l1) == "l1"] <- "Site1"
 colnames(l1)[colnames(l1) == "l2"] <- paste("Site", r, sep="")
 l1<-l1[,-1]
  }
l1<-round(l1)
l1[mapply(is.infinite, l1)] <- NA
l1[is.na(l1)]<-0
l1[l1<0]<-0
otu<-otu_table(l1, taxa_are_rows = T)
Sa<-sample_data(Factors)
out<-phyloseq(otu, Sa)
out
}

#' construct community (simple)
#' @param Comm1 list of species functions
#' @param Factors data frame of "environmental" data
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.comm()
make.comm<-function(Comm1, Factors){
otu<-matrix(data=NA, nrow=nrow(Factors), ncol = length(Comm1))
Sites<-c(paste0("Site", 1:30))
for(i in 1:length(Comm1)) {
  for(row in 1:nrow(Factors)){
   otu[row,i]<-do.call(Comm1[[i]], list(Factors[row,1],Factors[row,2],Factors[row,3],Factors[row,4],Factors[row,5]))
      }
}

#' subsample community
#' @param x phyloseq object
#' @param level single value or vector specifying sampling depth
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.rarefy2()
make.rarefy2<-function(x, level){
  require(phyloseq)
  require(vegan)

  sample_data(x)$adj<-level

  if (length(level)==1){
     p<-prune_samples(sample_sums(x)>level, x) # define samples we want to keep, discard rest

  if (nsamples(x)>nsamples(p)){warning(as.character(nsamples(x)-nsamples(p)), " samples have been removed because they are lower than rarefaction limit")}

  r<-as.data.frame(as.matrix(otu_table(p)))
  meta<-sample_data(p)
  rr<-rrarefy(t(r), meta$adj)
  ps<-phyloseq(otu_table(t(rr), taxa_are_rows = T), sample_data(meta))
  ps} else {

     p<-prune_samples(sample_sums(x)>level, x)
  if (nsamples(x)>nsamples(p)){warning(as.character(nsamples(x)-nsamples(p)), " samples have been removed because they are lower than rarefaction limit")}

    r<-as.data.frame(as.matrix(otu_table(p)))
  meta<-sample_data(p)
  rr<-rrarefy(t(r), meta$adj)
  ps<-phyloseq(otu_table(t(rr), taxa_are_rows = T), sample_data(meta))
  ps
  }
}

#' normalize routine using scaling
#' @param ps phyloseq object with community to be normalized
#' @param val mean value for sample scaling
#' @param scale vector of sample relative abundances for scaling
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.scaled2()
make.scaled2<-function(ps, val, scale){
  scaled<-data.frame(mapply(`*`, data.frame(as.matrix(otu_table(transform_sample_counts(ps, function(x) x/sum(x))))), scale * val))# sample_data(ps)$val))
  names<-rownames(data.frame(as.matrix(otu_table(ps))))
  rownames(scaled)<-names
  scaled<-round(scaled)

  p2<-ps
  otu_table(p2)<- otu_table(scaled, taxa_are_rows=T)
  p2
}

#' normalize routine using limma
#' @param ps phyloseq object with community to be normalized
#' @param Factor column from metadata for model matrix
#' @keywords limma variance stabilization
#' @export
#' @examples
#' make.limmaVST()
make.limmaVST<-function(ps, Factor){
  require(phyloseq)
  #ps<-filter_taxa(ps, function(x) sum(x)>0, T)
  counts<-as.data.frame(as.matrix(otu_table(ps)))
  factors<-unlist(sample_data(ps)[,Factor])
  design<-model.matrix(~factors)
  dge <- DGEList(counts=counts)
  dge <- calcNormFactors(dge) #what happens if we don't do this step?
  v<-voom(dge, design, plot=F)
  LimmaVST<-ps
  otu_table(LimmaVST)<-otu_table(v$E, taxa_are_rows = T)
  LimmaVST
}

#' normalize routine using deseq
#' @param ps phyloseq object with community to be normalized
#' @param Factor column from metadata for model matrix
#' @param l logical: 1 = log, other = linear
#' @keywords deseq variance stabilization
#' @export
#' @examples
#' make.deseqVST()
make.deseqVST<-function(ps, Factor, l=1){
r<-phyloseq_to_deseq2(ps, ~Factor)
geoMeans = apply(counts(r), 1, gm_mean)
dds = estimateSizeFactors(r, geoMeans = geoMeans)
#dds<-DESeqDataSetFromMatrix(r, sample_data(ps), design=~Factor)
#dds = estimateSizeFactors(dds)
if (l==1){dds = estimateDispersions(dds)} else {
dds <- estimateDispersionsGeneEst(dds)
 dispersions(dds) <- mcols(dds)$dispGeneEst
  }
vst = getVarianceStabilizedData(dds)
deseqVST<-ps
otu_table(deseqVST) <- otu_table(vst, taxa_are_rows = TRUE)
deseqVST
}

#' geometric mean function for deseq functions
#' @param x data table of community values
#' @keywords geometric mean
#' @export
#' @examples
#' gm_mean()
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#' standardizes the taxon abundance by taking counts relative to mean abundance of that taxon
#' @param ps phyloseq object with community to be normalized
#' @param ref reference community to ensure even lost taxa are included
#' @param method logical: 1 applies a exp() correction for log-transformed sample counts
#' @keywords relative abundance taxa
#' @export
#' @examples
#' Delta.sppcount()
Delta.sppcount<-function(ps,ref, method=0){
  tab<-as.data.frame(as.matrix(otu_table(ref)))
  tab[]<-0
  env<-sample_data(ref)
  tab<-phyloseq(otu_table(tab, taxa_are_rows = TRUE), sample_data(env))
  #make phyloseq object, subset phyloseq then merge second phyloseq.

  #merge_phyloseq()
  out<-as.data.frame(t(as.matrix(otu_table(ps))))
  out1<-as.data.frame(t(as.matrix(otu_table(ps))))
  if(method==1){
    out<-exp(out)
    out<-sapply(seq.int(dim(out)[2]), function(i) out[,i]/mean(out[,i]))
    out[is.nan(out)]<-0
    rownames(out)<-rownames(out1)
    colnames(out)<-colnames(out1)
    ps<-otu_table(out, taxa_are_rows = F)
    ps<-merge_phyloseq(ps, tab)
    ps
  } else {
    out<-sapply(seq.int(dim(out)[2]), function(i) out[,i]/mean(out[,i]))
    out[is.nan(out)]<-0
    rownames(out)<-rownames(out1)
    colnames(out)<-colnames(out1)
    ps<-otu_table(out, taxa_are_rows = F)
    ps<-merge_phyloseq(ps, tab)
    ps
  }
}

#' conducts a PERMANOVA specific for this model
#' @param ps phyloseq object
#' @keywords deseq variance stabilization
#' @export
#' @examples
#' make.PERMANOVA()
make.PERMANOVA<-function(ps){
  require(vegan)
  require(phyloseq)
  ps<-filter_taxa(ps, function(x) sum(x)>0, T)
  x<-as.data.frame(t(as.matrix(otu_table(ps))))
  y<-data.frame(as.matrix(sample_data(ps)))
  y$F1<-as.numeric(as.character(y$F1))
  y$F2<-as.numeric(as.character(y$F2))
  y$F3<-as.numeric(as.character(y$F3))
  y$F4<-as.numeric(as.character(y$F4))
  y$F5<-as.numeric(as.character(y$F5))
  data.frame(x,y)
  a<-adonis(x~Factor2+F1+F2+F3+F4+F5, data=y)
  a
}

#' extract PERMANOVA r2 values from a benchmarked object
#' @param tst product object from benchmark.MM()
#' @keywords PERMANOVA extract r squared
#' @export
#' @examples
#' ext.PERMANOVA()
ext.PERMANOVA<-function(tst){
out<-list("raw"=c(rep(NA, length(tst))), "scaled"=c(rep(NA, length(tst))), "pRare"=c(rep(NA, length(tst))), "eRare"=c(rep(NA, length(tst))), "deseqVST"=c(rep(NA, length(tst))), "limmaVST"=c(rep(NA, length(tst))))

  for(i in 1:length(tst)) {out$raw[i]<-test.script[[i]]$raw$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$scaled[i]<-test.script[[i]]$scaled$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$pRare[i]<-test.script[[i]]$pRare$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$eRare[i]<-test.script[[i]]$eRare$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$deseqVST[i]<-test.script[[i]]$deseqVST$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$limmaVST[i]<-test.script[[i]]$limma$PERMANOVA$aov.tab$R2[7]}
  out1<-as.data.frame(out)
  out2<-boxplot(out1)
  out4<-c("table"=out1,"plot"=out2)
  out4
}

#' extract Lost Information Index values from a benchmarked object
#' @param ps1.R reference phyloseq object
#' @param ps2.T phyloseq object with normalized community
#' @keywords LII lost information index extract r squared
#' @export
#' @examples
#' LII()
LII <-function(ps1.R, ps2.T){
  reference<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps1.R, ps1.R, method=0))))))
  reference2<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps1.R, ps1.R, method=0))))))
  treatment<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps2.T, ps1.R, method=0))))))

  Ci<-sapply(seq.int(dim(reference)[1]), function(i) summary(lm(reference[i,] ~ treatment[i,]))$r.squared)
  Ci2<-sapply(seq.int(dim(reference)[1]), function(i) summary(lm(reference[i,] ~ reference2[i,]))$r.squared)
  out<-sum(Ci2-Ci)
  out
}

#' extract Lost Information Index values from a benchmarked object
#' @param tst product object from benchmark.MM()
#' @keywords LII lost information index extract r squared
#' @export
#' @examples
#' ext.LII()
ext.LII<-function(tst){
out<-list("raw"=c(rep(NA, length(tst))), "scaled"=c(rep(NA, length(tst))), "pRare"=c(rep(NA, length(tst))), "eRare"=c(rep(NA, length(tst))), "deseqVST"=c(rep(NA, length(tst))), "limmaVST"=c(rep(NA, length(tst))))

  for(i in 1:length(tst)) {out$raw[i]<-test.script[[i]]$raw$LII}
  for(i in 1:length(tst)) {out$scaled[i]<-test.script[[i]]$scaled$LII}
  for(i in 1:length(tst)) {out$pRare[i]<-test.script[[i]]$pRare$LII}
  for(i in 1:length(tst)) {out$eRare[i]<-test.script[[i]]$eRare$LII}
  for(i in 1:length(tst)) {out$deseqVST[i]<-test.script[[i]]$deseqVST$LII}
  for(i in 1:length(tst)) {out$limmaVST[i]<-test.script[[i]]$limma$LII}
  out1<-as.data.frame(out)
  out1
}

#' test difference of linear model r2 values
#' @param ps1 product object from benchmark.MM()
#' @param ps2 product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' SVI()
SVI<-function(ps1, ps2){
  reference<-lm.test(ps1)
  trt<-lm.test(ps2)
  o<-reference$lm-trt$lm
  o
}

#' make linear model for each taxon
#' @param ps product object from benchmark.MM()
#' @keywords LII lost information index extract r squared
#' @export
#' @examples
#' lm.test()
lm.test<-function(ps){
anova.otu<-as.data.frame(t(as.matrix(otu_table(ps))))
anova.env<-data.frame(as.matrix(sample_data(ps)))
anova.env$F1<-as.numeric(as.character(anova.env$F1))
anova.env$F2<-as.numeric(as.character(anova.env$F2))
anova.env$F3<-as.numeric(as.character(anova.env$F3))
anova.env$F4<-as.numeric(as.character(anova.env$F4))
anova.env$F5<-as.numeric(as.character(anova.env$F5))
}

#' test difference of linear model r2 values
#' @param ps1.R phyloseq object of reference community
#' @param ps2.T phyloseq object of community to be tested
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' SVI2()
SVI2 <-function(ps1.R, ps2.T){
  reference<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps1.R, ps1.R, method=0))))))
  treatment<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps2.T, ps1.R, method=0))))))
  Ci<-sapply(seq.int(dim(reference)[1]), function(i) sum(abs(reference[i,] - treatment[i,])))
  names(Ci)<-rownames(reference)
  Ci
}

#' make indicators from deseq
#' @param ps product object from benchmark.MM()
#' @keywords deseq indicator species
#' @export
#' @examples
#' deseqIndics()
deseqIndics<-function(ps){
  require(phyloseq)
  require(DESeq2)
  r1<-subset_samples(ps, Factor=="one"|Factor=="two")
  r2<-subset_samples(ps, Factor=="one"|Factor=="three")
  r3<-subset_samples(ps, Factor=="one"|Factor=="four")
  r4<-subset_samples(ps, Factor=="one"|Factor=="five")
  r5<-subset_samples(ps, Factor=="one"|Factor=="six")
  tab<-ldply(list(r1,r2,r3,r4,r5), deseq.res)
  tab<-t(tab)
  tab[tab > 0.05]<-NA
  tab[tab < 0.05]<-1
  tab[is.na(tab)]<-0
  #tab[tab > 0.05]<-0
  colnames(tab)<-c("OneVTwo", "OneVThree", "OneVFour", "OneVFive", "OneVSix")
  tab
}

#' p values of deseq indicator species
#' @param x phyloseq object
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' deseq.res()
deseq.res<-function(x){
  sample_data(x)$Factor <- relevel(sample_data(x)$Factor, "one")
  r<-phyloseq_to_deseq2(x, ~Factor)
  geoMeans = apply(counts(r), 1, gm_mean)
  dds = estimateSizeFactors(r, geoMeans = geoMeans)
  dds<-DESeq(dds)
  res <- results(dds)
  res.p<-res$padj
  names(res.p)<-res@rownames
  res.p}

#' make indicator species
#' @param ps phyloseq object with community to be tested
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' indicspp()
indicspp<-function(ps){
  require(phyloseq)
  require(indicspecies)
  #ps<-filter_taxa(ps, function(x) sum(x)>0, T)
  p<-as.data.frame(as.matrix(t(otu_table(ps))))
  e<-sample_data(ps)$Factor

  m<-multipatt(p, e, control=how(nperm=999), duleg=TRUE)
  m.s<-m$sign
  m.s[is.na(m.s$p.value), 9]<-1.000
  m.s[m.s$p.value>0.05,1:6]<-0
  m.s[is.na(m.s)]<-0
  m.s
  }

#' extract summary from indicspp()
#' @param ps phyloseq
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' indicsummary()
indicsummary<-function(ps){
  require(phyloseq)
  require(indicspecies)
  #ps<-filter_taxa(ps, function(x) sum(x)>0, T)
  p<-as.data.frame(as.matrix(t(otu_table(ps))))
  e<-sample_data(ps)$Factor

  m<-multipatt(p, e, control=how(nperm=999))
  summary(m)
}

#' make limma indicators
#' @param ps1 product object from benchmark.MM()
#' @param Factor product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' limma.Indics3()
limma.Indics3<-function(ps, Factor){
    ps<-filter_taxa(ps, function(x) sum(x)>0, T)
    counts<-as.data.frame(as.matrix(otu_table(ps)))
    factors<-sample_data(ps)$Factor
    factors<-factor(factors, levels(factors)[c(3,6,5,2,1,4)])
    design<-model.matrix(~0+factors)
    contr.matrix<- makeContrasts(
    TwoVOne = factorstwo-factorsone,
    ThreeVOne = factorsthree-factorsone,
    FourVOne = factorsfour-factorsone,
    FiveVOne = factorsfive-factorsone,
    SixVOne = factorssix-factorsone,
    levels = colnames(design))
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge) #what happens if we don't do this step?
    v<-voom(dge, design, plot=F)
    fitV <- lmFit(v, design)
    fitV <- contrasts.fit(fitV, contrasts=contr.matrix)
    fitV <- eBayes(fitV, trend=TRUE)
    sig<-decideTests(fitV)
    sig
}

#' make limma indicators
#' @param ps1 product object from benchmark.MM()
#' @param Factor product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' limma.Indics()
limma.Indics<-function(ps, Factor){

    counts<-as.data.frame(as.matrix(otu_table(ps)))
    factors<-sample_data(ps)$Factor
    factors<-factor(factors, levels(factors)[c(6,12,10,4,3,8,7,1,5,9,2,11)])
    design<-model.matrix(~0+factors)
    contr.matrix<- makeContrasts(
    TwoVOne = factorstwo-factorsone,
    ThreeVOne = factorsthree-factorsone,
    FourVOne = factorsfour-factorsone,
    FiveVOne = factorsfive-factorsone,
    SixVOne = factorssix-factorsone,
    SevenVOne = factorsseven-factorsone,
    EightVOne = factorseight-factorsone,
    NineVOne = factorsnine-factorsone,
    TenVOne = factorsten-factorsone,
    ElevenVOne = factorseleven-factorsone,
    TwelveVOne = factorstwelve-factorsone,
    levels = colnames(design))
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge) #what happens if we don't do this step?
    v<-voom(dge, design, plot=F)
    fitV <- lmFit(v, design)
    fitV <- contrasts.fit(fitV, contrasts=contr.matrix)
    fitV <- eBayes(fitV, trend=TRUE)
    sig<-decideTests(fitV)
    sig
}

#' do filter protocol for limma
#' @param ps product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' filter.limma()
filter.limma<-function(ps){
    x<-as.data.frame(as.matrix(otu_table(ps)))
    group=sample_data(ps)$Factor2
    keep.exprs <- filterByExpr(x, group=group)
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    otu_table(ps)<-otu_table(x, taxa_are_rows = TRUE)
    ps
}

#' do filter protocol for limma
#' @param comm product object from benchmark.MM()
#' @param r product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' make.table()
make.table<-function(comm, r){

    m<-as.data.frame(t(table(sample(rownames(comm),rnorm(1, 250, 75), replace=T, prob=comm[,1]/sum(comm[,1])))))
     m<-m[,colnames(m)!="Var1"]
     colnames(m)[colnames(m)=="Freq"]<-paste("site", 1, r, sep=".")
      m2<-as.data.frame(t(table(sample(rownames(comm),rnorm(1, 250, 75), replace=T, prob=comm[,2]/sum(comm[,2])))))
     m<-merge(m, m2, by="Var2", all=T)
      m<-m[,colnames(m)!="Var1"]
      colnames(m)[colnames(m)=="Freq"]<-paste("site", 2, r, sep=".")
      for(i in 3:ncol(comm)){
       a<-as.data.frame(t(table(sample(rownames(comm),rnorm(1, 250, 75), replace=T, prob=comm[,i]/sum(comm[,i])))))
       a<-a[,colnames(a)!="Var1"]
       m<-merge(m,a, by="Var2", all=T)
        m<-m[,colnames(m)!="Var1"]
        colnames(m)[colnames(m)=="Freq"]<-paste("site", i, r, sep=".")

      }
      #colnames(m)<-c("Var1", paste0("site",c(1:12), r, sep="."))
        #names(m)[names(m)=="Var2"]<-paste0("Sample", r,)
        rownames(m)<-m$Var2
        m
  }



#############
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

####
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

# build my own functions for processing data

model.glm<-function(data,model){
  g<-glm(model, family=quasipoisson(link = "log"), data)
}

model.sp<-function(data, model){
  t<-apply(data, model.glm, model)
  }
#Build a script to run the following platforms:

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
