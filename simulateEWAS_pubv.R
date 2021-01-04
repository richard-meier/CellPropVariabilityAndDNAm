##########################################################################################
# Authors: D. Koestler, R. Meier, E. Nissen
# Date: January 1st, 2021
##########################################################################################


# === LIBRARIES AND FUNCTIONS === #

library(gtools)
library(EpiDISH)
library(matrixStats)
library(TCA) # version 1.2 or higher required !

create.combinations = function(...){
	x <- list(...)
	return(ccbn(x,1))
}
ccbn = function(xlist,n){
	if(n >= length(xlist)){
		out = data.frame(xlist[[length(xlist)]])
		colnames(out) = names(xlist)[length(xlist)]
		return( out )
	} else{
		base = ccbn(xlist,n+1)
		out = data.frame(xlist[[n]][1], base)
		colnames(out) = c( names(xlist)[n], colnames(base) )
		if( length(xlist[[n]])>1 ){
			for(i in 2:length(xlist[[n]])){
				tmp = data.frame(xlist[[n]][i], base)
				colnames(tmp) = c( names(xlist)[n], colnames(base) )
				out = rbind(out,tmp)
			}
		}
		return(out)
	}
}


# === CONSTANTS === #

# each of the 6 leukocytes will be a DMCT in REPS number of iterations
# so 6*25 = 150 total simulation iterations
REPS = 25
MCPGS = 500
MDIFFS = MCPGS*0.2
mdh = MDIFFS/2
NSMP = 100
npg = NSMP/2
cell_types = c("CD56+ NK cells","CD14+ Monocytes","CD19+ B cells","CD8+ T cells","CD4+ T cells","Neutrophils")
PROP_MEANS = c(0.062,0.077,0.082,0.096,0.178,0.505); names(PROP_MEANS)=cell_types
DMCT_COMBOS = 1:6
DELTA = 0.1

sim_scn = create.combinations(
	a0 = c(18,73,127),
	SDN = c(0.007, 0.02, 0.053),
	esize = c(0.05,0.1,0.2)
)


# === PREPARE OUTPUT TABLE COLUMNS === #

sim_scn$pow.cdmc = NA
sim_scn$t1e.cdmc = NA
sim_scn$derr.cdmc = NA
sim_scn$pow.tca = NA
sim_scn$t1e.tca = NA
sim_scn$derr.tca = NA

sim_scn$pow.cdmc.c1=NA
sim_scn$pow.cdmc.c2=NA
sim_scn$pow.cdmc.c3=NA
sim_scn$pow.cdmc.c4=NA
sim_scn$pow.cdmc.c5=NA
sim_scn$pow.cdmc.c6=NA
sim_scn$t1e.cdmc.c1=NA
sim_scn$t1e.cdmc.c2=NA
sim_scn$t1e.cdmc.c3=NA
sim_scn$t1e.cdmc.c4=NA
sim_scn$t1e.cdmc.c5=NA
sim_scn$t1e.cdmc.c6=NA

sim_scn$derr.cdmc.c1=NA
sim_scn$derr.cdmc.c2=NA
sim_scn$derr.cdmc.c3=NA
sim_scn$derr.cdmc.c4=NA
sim_scn$derr.cdmc.c5=NA
sim_scn$derr.cdmc.c6=NA

sim_scn$pow.tca.c1=NA
sim_scn$pow.tca.c2=NA
sim_scn$pow.tca.c3=NA
sim_scn$pow.tca.c4=NA
sim_scn$pow.tca.c5=NA
sim_scn$pow.tca.c6=NA
sim_scn$t1e.tca.c1=NA
sim_scn$t1e.tca.c2=NA
sim_scn$t1e.tca.c3=NA
sim_scn$t1e.tca.c4=NA
sim_scn$t1e.tca.c5=NA
sim_scn$t1e.tca.c6=NA

sim_scn$derr.tca.c1=NA
sim_scn$derr.tca.c2=NA
sim_scn$derr.tca.c3=NA
sim_scn$derr.tca.c4=NA
sim_scn$derr.tca.c5=NA
sim_scn$derr.tca.c6=NA


# === RUN SIMULATION === #

set.seed(123)
for(sss in 1:nrow(sim_scn)){

ALPHA0 = sim_scn$a0[sss]
t.effsize = sim_scn$esize[sss]
snoise = sim_scn$SDN[sss]

POWS_CDMC = c()
T1ES_CDMC = c()
POWS_TCA = c()
T1ES_TCA = c()
GAM_TRUE_1 = c()
GAM_TRUE_2 = c()
GAM_TRUE_3 = c()
GAM_TRUE_4 = c()
GAM_TRUE_5 = c()
GAM_TRUE_6 = c()
GAM_CDMC_1 = c()
GAM_CDMC_2 = c()
GAM_CDMC_3 = c()
GAM_CDMC_4 = c()
GAM_CDMC_5 = c()
GAM_CDMC_6 = c()
GAM_TCA_1 = c()
GAM_TCA_2 = c()
GAM_TCA_3 = c()
GAM_TCA_4 = c()
GAM_TCA_5 = c()
GAM_TCA_6 = c()

opc_by_ct_cdmc = list(
POWS_CT1 = c(),
POWS_CT2 = c(),
POWS_CT3 = c(),
POWS_CT4 = c(),
POWS_CT5 = c(),
POWS_CT6 = c(),
T1ES_CT1 = c(),
T1ES_CT2 = c(),
T1ES_CT3 = c(),
T1ES_CT4 = c(),
T1ES_CT5 = c(),
T1ES_CT6 = c()
)

opc_by_ct_tca = list(
POWS_CT1 = c(),
POWS_CT2 = c(),
POWS_CT3 = c(),
POWS_CT4 = c(),
POWS_CT5 = c(),
POWS_CT6 = c(),
T1ES_CT1 = c(),
T1ES_CT2 = c(),
T1ES_CT3 = c(),
T1ES_CT4 = c(),
T1ES_CT5 = c(),
T1ES_CT6 = c()
)

startTime = Sys.time()
cat("cprogr: ")
for(ctype in 1:6){
cat(".")

for(rrr in 1:REPS){

dmcts = ctype
non_dmcts = (1:6)[!((1:6) %in% dmcts)]
true_dmct_cpgs = 1:MDIFFS

G2_means = matrix(NA,nrow=MCPGS,ncol=6)
G1_means = matrix(NA,nrow=MCPGS,ncol=6)
for(idx in 1:mdh){
	for(h in 1:6){
		G2_means[idx,h] = runif(n=1,min=0.1,max=0.3)
		G1_means[idx,h] = G2_means[idx,h] + (h %in% dmcts)*t.effsize
	}
}
for(idx in (mdh+1):MDIFFS){
	for(h in 1:6){
		G2_means[idx,h] = runif(n=1,min=0.7,max=0.9)
		G1_means[idx,h] = G2_means[idx,h] - (h %in% dmcts)*t.effsize
	}
}
for(idx in (MDIFFS+1):MCPGS){
	for(h in 1:6){
		G2_means[idx,h] = runif(n=1,min=0.1,max=0.9)
		G1_means[idx,h] = G2_means[idx,h]
	}
}

GAMMA_TRUE = G1_means - G2_means

Z = list()
for(h in 1:6){
	Z[[h]] = matrix(NA,nrow=NSMP,ncol=MCPGS)
	for(i in 1:npg){
		Z[[h]][i,] = rnorm(MCPGS,mean=G1_means[,h],sd=snoise)
	}
	for(i in (npg+1):NSMP){
		Z[[h]][i,] = rnorm(MCPGS,mean=G2_means[,h],sd=snoise)
	}
	Z[[h]][Z[[h]]<0.0001]=0.0001
	Z[[h]][Z[[h]]>0.9999]=0.9999
}

true_props = rdirichlet(n=NSMP,alpha=ALPHA0*PROP_MEANS)
row.names(true_props) = paste0("SMP",1:NSMP)
colnames(true_props) = paste0("CL",1:6)

X = matrix(0,nrow=NSMP,ncol=MCPGS)
for(h in 1:6){
	for(i in 1:NSMP){
		X[i,] = X[i,] + Z[[h]][i,]*true_props[i,h]
	}
}
xbetas = t(X)
row.names(xbetas) = paste0("CP",1:MCPGS)
colnames(xbetas) = paste0("SMP",1:NSMP)


f = file(); sink(file=f,type="message")

cdmc.fit = CellDMC(beta.m=xbetas, pheno.v=rep(c(1,0),each=npg), frac.m=true_props, adjPMethod = "fdr")

sink(type="message"); close(f)

C1 = matrix(rep(c(1,0),each=npg),ncol=1)
colnames(C1)="TREAT"
row.names(C1) = row.names(true_props)

tca.mdl.A <- tca(
	X = xbetas, W = true_props, C1 = C1, C2 = NULL,
	parallel = TRUE, num_cores = 8, verbose = FALSE, 
	log_file = "C:\\Users\\r047m063\\Desktop\\tca_log.txt"
)


pow_cdmc = cdmc.fit$coe[[dmcts[1]]]$p[true_dmct_cpgs]
pow_cdmc = mean(pow_cdmc < 0.05)

t1e_cdmc = c(
	cdmc.fit$coe[[dmcts[1]]]$p[-true_dmct_cpgs], 
	cdmc.fit$coe[[non_dmcts[1]]]$p, cdmc.fit$coe[[non_dmcts[2]]]$p, 
	cdmc.fit$coe[[non_dmcts[3]]]$p, cdmc.fit$coe[[non_dmcts[4]]]$p,
	cdmc.fit$coe[[non_dmcts[5]]]$p
)
t1e_cdmc = mean(t1e_cdmc < 0.05)


pow_tca = tca.mdl.A$gammas_hat_pvals[,dmcts[1]][true_dmct_cpgs]
pow_tca = mean(pow_tca < 0.05)

t1e_tca = c(
	tca.mdl.A$gammas_hat_pvals[,dmcts[1]][-true_dmct_cpgs], 
	tca.mdl.A$gammas_hat_pvals[,non_dmcts[1]], tca.mdl.A$gammas_hat_pvals[,non_dmcts[2]], 
	tca.mdl.A$gammas_hat_pvals[,non_dmcts[3]], tca.mdl.A$gammas_hat_pvals[,non_dmcts[4]],
	tca.mdl.A$gammas_hat_pvals[,non_dmcts[5]]
)
t1e_tca = mean(t1e_tca < 0.05)


for(h in 1:6){
	if(h %in% dmcts){
		opc_by_ct_cdmc[[paste0("POWS_CT",h)]] = c(
			opc_by_ct_cdmc[[paste0("POWS_CT",h)]],
			cdmc.fit$coe[[h]]$p[true_dmct_cpgs] < 0.05
		)
		opc_by_ct_cdmc[[paste0("T1ES_CT",h)]] = c(
			opc_by_ct_cdmc[[paste0("T1ES_CT",h)]],
			cdmc.fit$coe[[h]]$p[-true_dmct_cpgs] < 0.05
		)
	} else{
		opc_by_ct_cdmc[[paste0("T1ES_CT",h)]] = c(
			opc_by_ct_cdmc[[paste0("T1ES_CT",h)]],
			cdmc.fit$coe[[h]]$p < 0.05
		)
	}
}

for(h in 1:6){
	if(h %in% dmcts){
		opc_by_ct_tca[[paste0("POWS_CT",h)]] = c(
			opc_by_ct_tca[[paste0("POWS_CT",h)]],
			tca.mdl.A$gammas_hat_pvals[,h][true_dmct_cpgs] < 0.05
		)
		opc_by_ct_tca[[paste0("T1ES_CT",h)]] = c(
			opc_by_ct_tca[[paste0("T1ES_CT",h)]],
			tca.mdl.A$gammas_hat_pvals[,h][-true_dmct_cpgs] < 0.05
		)
	} else{
		opc_by_ct_tca[[paste0("T1ES_CT",h)]] = c(
			opc_by_ct_tca[[paste0("T1ES_CT",h)]],
			tca.mdl.A$gammas_hat_pvals[,h] < 0.05
		)
	}
}

POWS_CDMC = c(POWS_CDMC,pow_cdmc)
T1ES_CDMC = c(T1ES_CDMC,t1e_cdmc)
POWS_TCA = c(POWS_TCA,pow_tca)
T1ES_TCA = c(T1ES_TCA,t1e_tca)
GAM_TRUE_1 = c(GAM_TRUE_1,GAMMA_TRUE[,1])
GAM_TRUE_2 = c(GAM_TRUE_2,GAMMA_TRUE[,2])
GAM_TRUE_3 = c(GAM_TRUE_3,GAMMA_TRUE[,3])
GAM_TRUE_4 = c(GAM_TRUE_4,GAMMA_TRUE[,4])
GAM_TRUE_5 = c(GAM_TRUE_5,GAMMA_TRUE[,5])
GAM_TRUE_6 = c(GAM_TRUE_6,GAMMA_TRUE[,6])
GAM_CDMC_1 = c(GAM_CDMC_1,cdmc.fit$coe[[1]]$Estimate)
GAM_CDMC_2 = c(GAM_CDMC_2,cdmc.fit$coe[[2]]$Estimate)
GAM_CDMC_3 = c(GAM_CDMC_3,cdmc.fit$coe[[3]]$Estimate)
GAM_CDMC_4 = c(GAM_CDMC_4,cdmc.fit$coe[[4]]$Estimate)
GAM_CDMC_5 = c(GAM_CDMC_5,cdmc.fit$coe[[5]]$Estimate)
GAM_CDMC_6 = c(GAM_CDMC_6,cdmc.fit$coe[[6]]$Estimate)
GAM_TCA_1 = c(GAM_TCA_1,tca.mdl.A$gammas_hat[,1])
GAM_TCA_2 = c(GAM_TCA_2,tca.mdl.A$gammas_hat[,2])
GAM_TCA_3 = c(GAM_TCA_3,tca.mdl.A$gammas_hat[,3])
GAM_TCA_4 = c(GAM_TCA_4,tca.mdl.A$gammas_hat[,4])
GAM_TCA_5 = c(GAM_TCA_5,tca.mdl.A$gammas_hat[,5])
GAM_TCA_6 = c(GAM_TCA_6,tca.mdl.A$gammas_hat[,6])

}

}
cat("\n")

G1 = cbind(GAM_TRUE_1,GAM_TRUE_2,GAM_TRUE_3,GAM_TRUE_4,GAM_TRUE_5,GAM_TRUE_6)
G2 = cbind(GAM_CDMC_1,GAM_CDMC_2,GAM_CDMC_3,GAM_CDMC_4,GAM_CDMC_5,GAM_CDMC_6)
G3 = cbind(GAM_TCA_1,GAM_TCA_2,GAM_TCA_3,GAM_TCA_4,GAM_TCA_5,GAM_TCA_6)
errorMat_cdmc = (G1-G2)^2
errorMat_tca = (G1-G3)^2
errVec_cdmc = sqrt(colMeans(errorMat_cdmc))
errVec_tca = sqrt(colMeans(errorMat_tca))

sim_scn$pow.cdmc[sss] = mean(POWS_CDMC)
sim_scn$t1e.cdmc[sss] = mean(T1ES_CDMC)
sim_scn$derr.cdmc[sss] = mean(errVec_cdmc)

sim_scn$pow.tca[sss] = mean(POWS_TCA)
sim_scn$t1e.tca[sss] = mean(T1ES_TCA)
sim_scn$derr.tca[sss] = mean(errVec_tca)

sim_scn$pow.cdmc.c1[sss] = mean(opc_by_ct_cdmc[["POWS_CT1"]])
sim_scn$pow.cdmc.c2[sss] = mean(opc_by_ct_cdmc[["POWS_CT2"]])
sim_scn$pow.cdmc.c3[sss] = mean(opc_by_ct_cdmc[["POWS_CT3"]])
sim_scn$pow.cdmc.c4[sss] = mean(opc_by_ct_cdmc[["POWS_CT4"]])
sim_scn$pow.cdmc.c5[sss] = mean(opc_by_ct_cdmc[["POWS_CT5"]])
sim_scn$pow.cdmc.c6[sss] = mean(opc_by_ct_cdmc[["POWS_CT6"]])

sim_scn$t1e.cdmc.c1[sss] = mean(opc_by_ct_cdmc[["T1ES_CT1"]])
sim_scn$t1e.cdmc.c2[sss] = mean(opc_by_ct_cdmc[["T1ES_CT2"]])
sim_scn$t1e.cdmc.c3[sss] = mean(opc_by_ct_cdmc[["T1ES_CT3"]])
sim_scn$t1e.cdmc.c4[sss] = mean(opc_by_ct_cdmc[["T1ES_CT4"]])
sim_scn$t1e.cdmc.c5[sss] = mean(opc_by_ct_cdmc[["T1ES_CT5"]])
sim_scn$t1e.cdmc.c6[sss] = mean(opc_by_ct_cdmc[["T1ES_CT6"]])

sim_scn$derr.cdmc.c1[sss] = sqrt(mean((GAM_TRUE_1 - GAM_CDMC_1)^2))
sim_scn$derr.cdmc.c2[sss] = sqrt(mean((GAM_TRUE_2 - GAM_CDMC_2)^2))
sim_scn$derr.cdmc.c3[sss] = sqrt(mean((GAM_TRUE_3 - GAM_CDMC_3)^2))
sim_scn$derr.cdmc.c4[sss] = sqrt(mean((GAM_TRUE_4 - GAM_CDMC_4)^2))
sim_scn$derr.cdmc.c5[sss] = sqrt(mean((GAM_TRUE_5 - GAM_CDMC_5)^2))
sim_scn$derr.cdmc.c6[sss] = sqrt(mean((GAM_TRUE_6 - GAM_CDMC_6)^2))

sim_scn$pow.tca.c1[sss] = mean(opc_by_ct_tca[["POWS_CT1"]])
sim_scn$pow.tca.c2[sss] = mean(opc_by_ct_tca[["POWS_CT2"]])
sim_scn$pow.tca.c3[sss] = mean(opc_by_ct_tca[["POWS_CT3"]])
sim_scn$pow.tca.c4[sss] = mean(opc_by_ct_tca[["POWS_CT4"]])
sim_scn$pow.tca.c5[sss] = mean(opc_by_ct_tca[["POWS_CT5"]])
sim_scn$pow.tca.c6[sss] = mean(opc_by_ct_tca[["POWS_CT6"]])

sim_scn$t1e.tca.c1[sss] = mean(opc_by_ct_tca[["T1ES_CT1"]])
sim_scn$t1e.tca.c2[sss] = mean(opc_by_ct_tca[["T1ES_CT2"]])
sim_scn$t1e.tca.c3[sss] = mean(opc_by_ct_tca[["T1ES_CT3"]])
sim_scn$t1e.tca.c4[sss] = mean(opc_by_ct_tca[["T1ES_CT4"]])
sim_scn$t1e.tca.c5[sss] = mean(opc_by_ct_tca[["T1ES_CT5"]])
sim_scn$t1e.tca.c6[sss] = mean(opc_by_ct_tca[["T1ES_CT6"]])

sim_scn$derr.tca.c1[sss] = sqrt(mean((GAM_TRUE_1 - GAM_TCA_1)^2))
sim_scn$derr.tca.c2[sss] = sqrt(mean((GAM_TRUE_2 - GAM_TCA_2)^2))
sim_scn$derr.tca.c3[sss] = sqrt(mean((GAM_TRUE_3 - GAM_TCA_3)^2))
sim_scn$derr.tca.c4[sss] = sqrt(mean((GAM_TRUE_4 - GAM_TCA_4)^2))
sim_scn$derr.tca.c5[sss] = sqrt(mean((GAM_TRUE_5 - GAM_TCA_5)^2))
sim_scn$derr.tca.c6[sss] = sqrt(mean((GAM_TRUE_6 - GAM_TCA_6)^2))

print(sim_scn[sss,])
print(Sys.time()-startTime)
cat("\n")

}


# === STORE RESULTS IN CSV FILE === #

outfile = "...targetDirectory\\simulation_result_table.csv"
write.csv(sim_scn,file=outfile,row.names=FALSE)

