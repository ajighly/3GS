library(BGLR)
nIter = 10000
burnIn= 5000

trnpheno=read.table('est')
pheno=read.table('FinalYLD.txt')
geno=read.table('FinalGeno.txt',header=T)
trngeno=geno[,row.names(trnpheno)]
trngeno=t(trngeno)
trngeno=scale(x=trngeno,center=T,scale=F)

valpheno=read.table('FinalYLD.txt')
geno=t(geno)
valgeno=geno[row.names(valpheno),]

ETA=list(MRK=list(X = trngeno, model = "BRR"))
for (i in 1:ncol(pheno)) {
	trnphe=trnpheno[,i]
	FF=BGLR (y=trnphe,ETA=ETA,nIter=nIter,burnIn=burnIn,saveAt=colnames(pheno)[i],verbose=F)
	Effects=as.matrix(FF$ETA$MRK$b)
	row.names(Effects)=colnames(trngeno)
	write.table(Effects,paste0(colnames(pheno)[i],".Effects"),quote=F,sep="\t",col.names=F)
	GEBV2=as.matrix(trngeno)%*%as.matrix(Effects[,1])
	GEBV=as.matrix(valgeno)%*%as.matrix(Effects[,1])
	GEBV=(GEBV*lm(FF$y~as.numeric(GEBV2))$coefficients[2])+(lm(FF$y~as.numeric(GEBV2))$coefficients[1])
	write.table(GEBV,paste0(colnames(pheno)[i],".GEBV"),quote=F,sep="\t",col.names=F)
}
