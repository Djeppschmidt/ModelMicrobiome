### Package contents:

# pipeline: ####
#select spp-> model comm -> sample comm -> pre-process -> make inference -> evaluate inference (make plots)

### All Species Functions
 ### eventually both niche-modeled and interactions

# spp key ####

makeKey<-function(spp, type, response){
  df<-data.frame(names(spp), type, response)
  df
}

# "assembly" tool ####
#(started)
makeComm <- function(spp, Factors){

  out<-matrix(data=NA, nrow=nrow(Factors), ncol = length(spp))

  for(i in 1:length(spp)) {
    for(row in 1:nrow(Factors)){
      out[row,i]<-do.call(spp[[i]], list(Factors[row,1],Factors[row,2],Factors[row,3],Factors[row,4],Factors[row,5]))
    }
  }
  out
}

assemblePS<-function(spp, Factors){
  otu<-makeComm(spp, Factors)
  ps<-merge_phyloseq(otu, Factors)
  ps
}

# a rarefaction Function ####

rarefy2<-
  function (x, sample, replace) #added replace argument as input
  {
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
      row <- sample(rep(nm, times = x[i, ]), sample[i], replace) #use replace argument
      row <- table(row)
      ind <- names(row)
      x[i, ] <- 0
      x[i, ind] <- row
    }
    x
  }

# normalization wrappers: CSS, VST, limma, rarefy_even, rarefy_uneven, rarefy_adj ####

NormAdj<-function(OTU, qpcr.val, set){
  require(vegan)
  a<-sum(qpcr.val)/length(qpcr.val)
  b<-qpcr.val/a
  val<-b*set
  otu<-rarefy(OTU, val)
  otu
} #works for PFLA as long as it's in a vector form of amounts (not proportion)
#should be aggregated to all bacteria level... sum? think about this more.

rarefy_uneven<-function(ps, set){
  require(phyloseq)
  OTU<-otu_table(ps)
  qpcr.val<-sample_data(ps)$qpcr
  otu<-NormAdj(OTU, qpcr.val, set)
  out
}#for multiple set points, iterate within the model run
rarefy_even<-function(ps, val){
  require(phyloseq)
  require(vegan)
  OTU<-otu_table(ps)
  otu<-rarefy(OTU, val)
  otu
} # val can be a single value

#etagenomeSeq object####
#make metagenomeSeq object from phyloseq, calculate cumsum
#sourcecode: http://joey711.github.io/waste-not-supplemental/simulation-differential-abundance/simulation-differential-abundance-server.html
make_metagenomeSeq = function(physeq) {
  require("metagenomeSeq")
  require("phyloseq")
  # Enforce orientation
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  OTU = as(otu_table(physeq), "matrix")
  # Convert sample_data to AnnotatedDataFrame
  ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))
  # define dummy 'feature' data for OTUs, using their name Helps with
  # extraction and relating to taxonomy later on.
  TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), row.names = taxa_names(physeq)))
  # Create the metagenomeSeq object
  MGS = newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)
  # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
  MGS = cumNorm(MGS)
  return(MGS)
}#makes a CSS normalized metagenomseq object


# get MGS -> phyloseq

#
CSSTest<-function(ps, mod){
  pd<-make_metagenomeSeq(ps)
  pd <- pData(pd)#not sure if this is the right next step
  mod <- model.matrix(mod, data = pd)
  res = fitFeatureModel(res, mod)
  coef<-MRcoefs(res)
  coef
}

# extract run metrics ####

# evaluate differential abundance ###
evalDiff<- function(x, key){
  if names(x)!= names(key) {return("check order")}
  df<-data.frame(key, ...) #make logical data.frame from p values
  df
}

# gm_mean ####

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# DESeq2 ####

r.Deseq2<-function(data, model, Type){
  require(DESeq2)
  require(edgeR)
  dds = phyloseq_to_deseq2(data, model)

  geoMeans = apply(counts(dds), 1, gm_mean)

  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  dds = DESeq(dds, fitType=Type)

  res = results(dds)
  res
}#include categories/predictor variables in phyloseq; if a factor, specify
#this one works!

summarizeDeseq2<-function(res){
  T.P<-
  F.P<-
  F.N<-
  T.N<-
}
permuteDESeq2<-function(data, model, Type){
  lapply(runDeseq2, levels, ...)
  runDeseq2(data, model, fitType)
}

#limma-Voom ####
#load dependencies
library(limma)
library(Glimma)
library(edgeR)
r.limma<-function(ps, model){
  require(phyloseq)
  otutab<-as.matrix(t(otu_table(ps)))
  sdat<-sample_data(ps)
  dge<-DGEList(counts=otutab)
  dge<-calcNormFactors(dge)
  design<-model.matrix(model)
  v<-voom(dge, design)
  fit<-lmFit(v, design)
  contr<-makeContrasts(categories, levels=colnames(coef(fit)))
  tmp<-contrasts.fit(fit, contr)
  tmp
}

#almost working.. needs to be truthed.


# benchmark tool ####
benchmark<-function(df, array){}

#

### summary diveristy stats function
# species distribution curve
# skewness value?
# alpha diversity
