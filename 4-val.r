library(nipals)
library(gge)
source('gge.r')

pheno=read.table('FinalYLD.txt')
NoTraits=ncol(pheno)
valpheno=read.table('val')
estpheno=read.table('est')
GEBVs=matrix(NA,nrow=nrow(valpheno),ncol=NoTraits)
for (j in 1:NoTraits) {
	x=read.csv(paste0("PC",j,".GEBV"),sep='\t',header=F)
	GEBVs[,j]=x[1:nrow(valpheno),2]
}
YLD=as.matrix(read.table("stRefGEBVs.txt"))
YLD=YLD[row.names(estpheno),]
m = gge(YLD,maxiter=1e7)
pc=prcomp(m$x)
mx=t(t(GEBVs)*pc$sdev) %*% solve(t(t(m$locCoord)/pc$sdev))
row.names(mx)=row.names(valpheno)
write.table(mx,'PredictedPheno',quote=F,sep='\t')
write.table(diag(cor(mx,valpheno,use='pairwise.complete.obs')),'accuracy',quote=F,sep='\t')
