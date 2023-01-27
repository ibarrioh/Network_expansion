
#############################################################################################
#############################
#############################

##Correlations & distances are calculated both with and without normalization
##In the paper we only use results that are NOT norm

setwd("/tables_expansion/RDS_astro/")

path=paste(	"/tables_expansion/RDS_astro/",list.files(pattern= ".rds"),sep="")

gwas_all=matrix(0,18410,length(path))

colnames(gwas_all)=1:ncol(gwas_all)

for (i in 1:length(path)){
  
  temp=readRDS(path[i])
  gwas_all[,i]=as.numeric(temp[,"page.rank"])
  colnames(gwas_all)[i]=temp[1,"Trait"]
}


rownames(gwas_all)=temp[,"ENSG"]

###Due to size constrains this file is not present in the repository
saveRDS(gwas_all,"/tables_expansion/gwas_astro_pageRank.rds")


Zscore=gwas_all

for (i in 1:nrow(Zscore)){
  
  Zscore[i,]=(gwas_all[i,]-mean(gwas_all[i,]))/sd(gwas_all[i,])
  
}

Zscore=Zscore[!is.na(Zscore[,1]),]


#############################
#############################

####Correlations (not Norm)

###Pearson

corr=cor(gwas_all,method="pearson")

ut <- upper.tri(corr)
corr=data.frame(	i = rownames(corr)[row(corr)[ut]],
                 j = rownames(corr)[col(corr)[ut]],
                 cor=t(corr)[ut])

corr.traits.norm=corr

colnames(corr.traits.norm)[ncol(corr.traits.norm)]="corr.pearson"

###Spearman

corr=cor(gwas_all,method="spearman")

ut <- upper.tri(corr)
corr=data.frame(	i = rownames(corr)[row(corr)[ut]],
                 j = rownames(corr)[col(corr)[ut]],
                 cor=t(corr)[ut])

corr.traits.norm=cbind(corr.traits.norm,corr[,"cor"])

colnames(corr.traits.norm)[ncol(corr.traits.norm)]="corr.spearman"


####Correlations (Zscores)

###Pearson

corr=cor(Zscore,method="pearson")

ut <- upper.tri(corr)
corr=data.frame( i = rownames(corr)[row(corr)[ut]],
                 j = rownames(corr)[col(corr)[ut]],
                 cor=t(corr)[ut])

corr.traits.zsco=corr

colnames(corr.traits.zsco)[ncol(corr.traits.zsco)]="corr.pearson"

###Spearman

corr=cor(Zscore,method="spearman")

ut <- upper.tri(corr)
corr=data.frame(	i = rownames(corr)[row(corr)[ut]],
                 j = rownames(corr)[col(corr)[ut]],
                 cor=t(corr)[ut])

corr.traits.zsco=cbind(corr.traits.zsco,corr[,"cor"])

colnames(corr.traits.zsco)[ncol(corr.traits.zsco)]="corr.spearman"

saveRDS(corr.traits.zsco,"/tables_expansion/corr.traits.zsco.rds")
saveRDS(corr.traits.norm,"/tables_expansion/corr.traits.norm.rds")

#############################


#####Distances matrix (not norm & Zscore)

dist.norm=data.frame(	t(combn(colnames(gwas_all),2)), 
                      as.numeric(dist(t(gwas_all),method="euclidean")),
                      as.numeric(dist(t(gwas_all),method="manhattan")),
                      as.numeric(dist(t(gwas_all),method="canberra")))

names(dist.norm) <- c("i","j","euclidean","manhattan","canberra")

dist.Zsco=data.frame(	t(combn(colnames(Zscore),2)), 
                      as.numeric(dist(t(Zscore),method="euclidean")),
                      as.numeric(dist(t(Zscore),method="manhattan")),
                      as.numeric(dist(t(Zscore),method="canberra")))

names(dist.Zsco) <- c("i","j","euclidean","manhattan","canberra")

saveRDS(dist.Zsco,"/tables_expansion/dist.Zsco.rds")
saveRDS(dist.norm,"/tables_expansion/dist.norm.rds")

#####Cosine distances

cos.sim <- function(ix) 
{
  A = X[ix[1],]
  B = X[ix[2],]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
} 

##not norm

X=t(gwas_all)
n <- nrow(X) 
cmb <- expand.grid(i=1:n, j=1:n) 
cos.norm <- matrix(apply(cmb,1,cos.sim),n,n)

colnames(cos.norm)=colnames(gwas_all)
rownames(cos.norm)=colnames(gwas_all)
ut <- lower.tri(cos.norm)
cos.norm=data.frame(	i = rownames(cos.norm)[row(cos.norm)[ut]],
                     j = rownames(cos.norm)[col(cos.norm)[ut]],
                     cos.dist=t(cos.norm)[ut])

####Zscore

X=t(Zscore)
n <- nrow(X) 
cmb <- expand.grid(i=1:n, j=1:n) 
cos.Zsco<- matrix(apply(cmb,1,cos.sim),n,n)

colnames(cos.Zsco)=colnames(Zscore)
rownames(cos.Zsco)=colnames(Zscore)
ut <- lower.tri(cos.Zsco)
cos.Zsco=data.frame(	i = rownames(cos.Zsco)[row(cos.Zsco)[ut]],
                     j = rownames(cos.Zsco)[col(cos.Zsco)[ut]],
                     cos.dist=t(cos.Zsco)[ut])

saveRDS(cos.Zsco,"/tables_expansion/cos.Zsco.rds")
saveRDS(cos.norm,"/tables_expansion/cos.norm.rds")


