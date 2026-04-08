library(nipals)
library(gge)
source('gge.r')

geno=read.table('FinalGeno.txt',header=T)
pheno=read.table('est')
geno2=geno[,row.names(pheno)]
trngeno=t(geno2)
trngeno=scale(x=trngeno,center=T,scale=F)
trngeno[is.na(trngeno)]=0
GEBVs=NULL
for (i in 1:ncol(pheno)) {
	x=read.table(paste0(colnames(pheno)[i],".Effects"))
	GEBV=trngeno%*%as.matrix(x[,2])
	GEBV=(GEBV*lm(pheno[,i]~as.numeric(GEBV))$coefficients[2])+(lm(pheno[,i]~as.numeric(GEBV))$coefficients[1])
	GEBVs=cbind(GEBVs,GEBV[,1])
}
colnames(GEBVs)=colnames(pheno)
GEBVs=GEBVs[row.names(pheno),]
write.table(GEBVs,'stRefGEBVs.txt',quote=F,sep='\t')


YLD=as.matrix(read.table("stRefGEBVs.txt"))
YLD2=as.matrix(read.table('est'))
m = gge(YLD,maxiter=1e7)
pcs=m$genCoord
pcs[is.na(YLD2)]=NA
colnames(pcs)=colnames(YLD)
write.table(pcs,'pcs',quote=F,sep="\t")


