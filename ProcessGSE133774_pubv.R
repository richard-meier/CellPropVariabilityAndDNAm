##########################################################################################
# Authors: D. Koestler, R. Meier, E. Nissen
# Date: January 1st, 2021
##########################################################################################


##########################################################################################
# 0a: install and load necessary R/Bioconductor packages
##########################################################################################

library(Biobase)
library(GEOquery)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(quadprog)
library(sirt)


##########################################################################################
# 0b: load the source code for deconvolution and Dirichlet maximum likelihood estimation
##########################################################################################

load("...targetDirectory/ObjectsForDeconvolution_pubv.RData")
source("...targetDirectory/DeconvolutionCode_pubv.R")


##########################################################################################
# 1a: read in the GEO data, extract the beta values and covariates
##########################################################################################

gse.id = "GSE133774"
dataset_directory = paste0("...targetDirectory/",gse.id)
if(!dir.exists(dataset_directory)){
  dir.create(dataset_directory)
}
gset <- getGEO(
  gse.id, GSEMatrix =TRUE, getGPL=FALSE,
  destdir = dataset_directory
)
if (length(gset) > 1) idx <- grep("GPL21145", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

betas = exprs(gset)
covariates = pData(gset)


##########################################################################################
# 1b: subset the data appropriately, e.g., remove any non-whole blood samples and save
# the resulting beta values and covariates to file
##########################################################################################

# remove SNP associated CpGs 
annotation = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotation = annotation[is.na(annotation$Probe_rs),]
annotation = annotation[is.na(annotation$CpG_rs),]
annotation = annotation[is.na(annotation$SBE_rs),]
betas = betas[row.names(betas) %in% row.names(annotation),]

# remove CpGs with missing values
betas = na.omit(betas)

# if "covariates" contains samples that are not whole blood exclude them here !

# only retain samples that are included in both "betas" and "covariates"
betas = betas[,colnames(betas) %in% row.names(covariates)]
covariates = covariates[row.names(covariates) %in% colnames(betas),]

# ensure that the sample order in "betas" and "covariates" is identical
covariates = covariates[match(colnames(betas),row.names(covariates)),]


##########################################################################################
# 2: explore and clean the "covariate" table 
##########################################################################################
# In this section, we want to do both:
#    1. identify covariates that could explain variation in beta values
#    2. reformat covariates for analysis if necessary

covariates$age = as.numeric(covariates$"age:ch1")
hist(covariates$age)

covariates$gender=covariates$"gender:ch1"
table(covariates$gender)

covariates$disease_state=covariates$"disease state:ch1"
table(covariates$disease_state)

# store the clean dataset in a separate file
save(betas, covariates, file = paste0(dataset_directory, "/", gse.id, "_data.RData"))


##########################################################################################
# 3a: deconvolute cellular composition 
##########################################################################################

beta_subset = betas[row.names(betas) %in% IDOL.optim.DMRs.EPIC,]
cellprops = PredictCellComposition(Betas=beta_subset, Array = "EPIC")


##########################################################################################
# 3b: check for outliers across each of the six cell types using the IQR technique
##########################################################################################

quartiles = apply(cellprops,MARGIN=2,FUN=quantile,probs=c(0.25,0.75))
IQRs = quartiles[2,] - quartiles[1,]
bounds = quartiles + rbind(-IQRs*1.5,IQRs*1.5)
checkOutlier = function(k,vals,bounds){
  (vals[,k] < bounds[1,k]) | (vals[,k] > bounds[2,k])
}
outlierMatrix = sapply(1:6, FUN=checkOutlier, vals=cellprops, bounds=bounds)


##########################################################################################
# 3c: save to file, the deconvolution estimates and the matrix indicating flags for potential
# outliers
##########################################################################################

save(cellprops, outlierMatrix, file = paste0(dataset_directory, "/", gse.id, "_deconvolution.RData"))


##########################################################################################
# 4: Fit a Dirichlet model to the deconvolution estimates and estimate the Dirichlet
# parameters, estimate summary statistics, and save to file.
##########################################################################################

output_directory = "...targetDirectory/summary_statistics"

is_outlier = apply(outlierMatrix,MARGIN=1,FUN=sum) > 1
cellprops = cellprops[!is_outlier,]

# obtain dirichlet MLEs and save them
dmle = dirichlet.mle(cellprops)
summaryDirichlet = rbind(
  c(dmle$alpha0, dmle$alpha)
)
summaryDirichlet = round(summaryDirichlet,digits=4)
colnames(summaryDirichlet)[1] = "0"
colnames(summaryDirichlet) = paste0("alpha.",colnames(summaryDirichlet))
summaryDirichlet = data.frame(
  data.set = gse.id,
  data.type = "Full Dataset",
  N = c(nrow(cellprops)),
  summaryDirichlet
)
write.csv(summaryDirichlet,file=paste0(output_directory,"/SMRY_Dirichlet_",gse.id,".csv"),row.names=FALSE)

# obtain summary statistics and save them
meanprops = rbind(
  apply(cellprops,MARGIN=2,FUN=mean)
)
meanprops = round(meanprops,digits=4)
colnames(meanprops) = paste0("mean.prop.",colnames(meanprops))
sdprops = rbind(
  apply(cellprops,MARGIN=2,FUN=sd)
)
sdprops = round(sdprops,digits=4)
colnames(sdprops) = paste0("sd.prop.",colnames(sdprops))
summaryCellprops = data.frame(
  data.set = gse.id,
  data.type = "Full Dataset",
  N = c(nrow(cellprops)),
  meanprops,sdprops
)
write.csv(summaryCellprops,file=paste0(output_directory,"/SMRY_Cellprop",gse.id,".csv"),row.names=FALSE)


##########################################################################################
# 5: Obtain noise level estimates and save to file.
##########################################################################################

# create a function that performs noise level estimation
#   -> NOTE: covariates need to be modified according to the dataset !
estimateSD = function(Y,meta){
  return( summary(lm(Y~age+gender+disease_state,data=meta))$sigma )
}

# estimate noise levels (this will take several minutes)
SD_ESTIMATES = apply(betas,MARGIN=1,FUN=estimateSD,meta=covariates)
save(SD_ESTIMATES, file = paste0(dataset_directory, "/", gse.id, "_noise.RData"))

# create and save noise level summary
summaryMatrix = round(matrix(
  quantile(SD_ESTIMATES,probs=c(0.01,0.05,(1:9)/10,0.95,0.99)), nrow=1
),digits=4)
colnames(summaryMatrix) = paste0("Q",c("01","05",10*(1:9),95,99))
summaryMatrix = data.frame(
  data.set = gse.id,
  data.type = "Full Dataset",
  N = c(nrow(covariates)),
  summaryMatrix,
  mean = round(c(mean(SD_ESTIMATES)),digits=4),
  sd = round(c(sd(SD_ESTIMATES)),digits=4)
)
write.csv(summaryMatrix,file=paste0(output_directory,"/SMRY_SD_",gse.id,".csv"),row.names=FALSE)

