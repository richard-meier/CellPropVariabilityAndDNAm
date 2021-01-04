##########################################################################################
# Authors: D. Koestler, R. Meier, E. Nissen
# Date: January 1st, 2021
##########################################################################################


GSE_IDS = c(
"GSE145233", "GSE149318", "GSE155426", "GSE123914", 
"GSE125367", "GSE131461", "GSE133774", "GSE140038"
)

files = paste0(
	"...targetDirectory/summary_statistics/SMRY_Dirichlet_",
	GSE_IDS, ".csv"
)

dat1 = read.csv(files[1])
for(i in 2:length(files)){
	dat1 = rbind(dat1,read.csv(files[i]))
}

files = paste0(
	"...targetDirectory/summary_statistics/SMRY_SD_",
	GSE_IDS, ".csv"
)

dat2 = read.csv(files[1])
for(i in 2:length(files)){
	dat2 = rbind(dat2,read.csv(files[i]))
}

eprops = dat1[,-(1:4)] / dat1[,4]

cleanNoiseTable = dat2[,c("data.set","data.type","N","Q10","Q50","Q90")]
cleanNoiseTable$Q10 = round(cleanNoiseTable$Q10,digits=3)
cleanNoiseTable$Q50 = round(cleanNoiseTable$Q50,digits=3)
cleanNoiseTable$Q90 = round(cleanNoiseTable$Q90,digits=3)
cleanNoiseTable = rbind(
	cleanNoiseTable,
	data.frame(data.set="AVERAGE",data.type="----------------",N="---",
	Q10=round(mean(cleanNoiseTable[,4]),digits=3),
	Q50=round(mean(cleanNoiseTable[,5]),digits=3),
	Q90=round(mean(cleanNoiseTable[,6]),digits=3)
))

cleanDirichletTable = dat1
eprops = round( dat1[,-(1:4)] / dat1[,4] , digits=3)
cleanDirichletTable[,-(1:4)] = eprops
cleanDirichletTable$alpha.0 = round(cleanDirichletTable$alpha.0,digits=3)
colnames(cleanDirichletTable) = c("data.set", "data.type", "N", "alpha.0", "W.CD4T", "W.CD8T", "W.Bcell", "W.NK", "W.Monocyte", "W.Neutrophil")
cleanDirichletTable = rbind(
	cleanDirichletTable,
	data.frame(data.set="AVERAGE",data.type="----------------",N="---",
	alpha.0=round(mean(cleanDirichletTable[,4]),digits=3),
	W.CD4T=round(mean(cleanDirichletTable[,5]),digits=3),
	W.CD8T=round(mean(cleanDirichletTable[,6]),digits=3),
	W.Bcell=round(mean(cleanDirichletTable[,7]),digits=3),
	W.NK=round(mean(cleanDirichletTable[,8]),digits=3),
	W.Monocyte=round(mean(cleanDirichletTable[,9]),digits=3),
	W.Neutrophil=round(mean(cleanDirichletTable[,10]),digits=3)
))

write.csv(cleanNoiseTable,file="...targetDirectory/summary_statistics/Overall_Noise_Summary.csv",row.names=FALSE)
write.csv(cleanDirichletTable,file="...targetDirectory/summary_statistics/Overall_Dirichlet_Summary.csv",row.names=FALSE)

