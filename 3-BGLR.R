library(BGLR)
nIter = 10000
burnIn= 5000

pheno=read.table('pcs')
phenoOrig=read.table('FinalYLD.txt')
geno=read.table('FinalGeno.txt',header=T)
geno=geno[,row.names(phenoOrig)]

geno=t(geno)
#row.names(geno)=geno[,1]
#geno= geno[,-1]
for (trait in 1:ncol(pheno)) {
	trngeno=geno[row.names(pheno)[!is.na(pheno[,trait])],]
	trngeno=scale(x=trngeno,center=T,scale=F)
	trngeno[is.na(trngeno)]=0

	trnpheno=pheno[row.names(pheno)[!is.na(pheno[,trait])],]
	trnpheno=trnpheno[,trait]
	ETA=list(MRK=list(X = trngeno, model = "BRR"))
	FF=BGLR (y=trnpheno,ETA=ETA,nIter=nIter,burnIn=burnIn,saveAt=paste0("PC",trait),verbose=F)

	valpheno=read.table('val')
	valgeno=geno[row.names(valpheno),]
	Effects=as.matrix(FF$ETA$MRK$b)
	row.names(Effects)=colnames(trngeno)
	write.table(Effects,paste0("PC",trait,".Effects"),quote=F,sep="\t",col.names=F)
	GEBV2=as.matrix(trngeno)%*%as.matrix(Effects[,1])
	GEBV=as.matrix(valgeno)%*%as.matrix(Effects[,1])
	GEBV=(GEBV*lm(FF$y~as.numeric(GEBV2))$coefficients[2])+(lm(FF$y~as.numeric(GEBV2))$coefficients[1])
	write.table(GEBV,paste0("PC",trait,".GEBV"),quote=F,sep="\t",col.names=F)
}
