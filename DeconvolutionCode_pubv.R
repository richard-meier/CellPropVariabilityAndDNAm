#############################################################################################
# Authors: D. Koestler
# Date: January 1st, 2021
#
# FUNCTION:  PredictCellComposition 
#    This function estimates the proportion of the six normal leukocyte cell types: CD4T, 
#    CD8T, NK, Bcell, Monocytes, and Neutrophils (Granulocytes) based on the IDOL optimized
#    library for the Illumina HumanMethylation450 array (aka 450K array) (Koestler et al., 2016)
#    or the IDOL optimized library for the Illumina HumanMethylationEPIC array 
#    (Salas and Koestler et al., 2018).  
#
# REQUIRES:	quadprog    
#
# ARGUMENTS:
#             
#	Betas:  		 A J x N matrix of whole blood methylation beta-values; J represents 
#                    the number of CpGs (i.e., ~ 450,000 for the Illumina HumanMethylation450
#                    array or ~850,000 for the EPIC array) and N represents the number of samples 
#
#   Array:           A character vector denoting the array technology used to profile DNA methylation
#                    for the whole-blood DNA methylation data set.  Possible values inlcude "450K" for 
#                    the Illumina HumanMethylation450 array and "EPIC" for the Illumina 
#                    HumanMethylationEPIC array
#
# RETURNS:           An N x 6 matrix with the cellular composition estimates.  Each column represents
#                    a different cell type  
#############################################################################################

PredictCellComposition <- function(Betas, Array = NULL) {
	require(quadprog)
	if(is.null(Array)) stop("user must supply the array technology used!")
	
	if(Array == "450K") {
		# check the overlap in the number of probes in the whole-blood data set with the 
		# 450K reference library
		overlap = intersect(rownames(Betas), IDOL.optim.DMRs.450K)
		print(paste(length(overlap), " out of ", length(IDOL.optim.DMRs.450K), " probes used for",
		      " deconvolution", sep = ""))
		cellest = projectWBCnew(Betas[overlap, ], IDOL.optim.coefEsts.450K[overlap,])
	}
		if(Array == "EPIC") {
		# check the overlap in the number of probes in the whole-blood data set with the 
		# EPIC reference library
		overlap = intersect(rownames(Betas), IDOL.optim.DMRs.EPIC)
		print(paste(length(overlap), " out of ", length(IDOL.optim.DMRs.EPIC), " probes used for",
		      " deconvolution", sep = ""))
		cellest = projectWBCnew(Betas[overlap, ], IDOL.optim.coefEsts.EPIC[overlap,])
	}
	
	return(cellest)
	
}
      