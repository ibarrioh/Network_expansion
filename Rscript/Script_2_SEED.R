
################################################################################################
library(igraph)
library(pROC)


anotation<-function(node.list,genes,disease){
  
  genes=genes[genes[,"disease"]%in%disease,]
  
  for (i in 1:nrow(node.list)){
    if(node.list[i,"ENSG"]%in%genes[,"gene"]){
      node.list[i,"padj"]=max(as.numeric(genes[genes[,"gene"]%in%node.list[i,"ENSG"],"padj"]))
    }else{
      node.list[i,"padj"]=0
    }}
  
  return(node.list)
  
}


################################################################################

string=		readRDS("/tables_expansion/Combined_STRINGv11_OTAR281119_FILTER.rds")
EFO_DOI=as.matrix(read.delim("/tables_expansion/EFO_DOI_long_matrix.csv",sep=","))

EFO_DOI=cbind(EFO_DOI,NA,NA)
colnames(EFO_DOI)[(ncol(EFO_DOI)-1):ncol(EFO_DOI)]=c("Zscore_nodes","Zscore_TP")

setwd("/tables_expansion/RDS_astro")


path=list.files(pattern= ".rds")
path=gsub(".rds","",path)
path=gsub("nodes.","",path)
path=path[path%in%EFO_DOI[,"EFO"]]

path=cbind(path,paste(	"/tables_expansion/RDS_astro/nodes.",path,".rds",sep=""))
colnames(path)=c("EFO","path")

################################################################################
####Loop clusters
i <- as.numeric(commandArgs(trailingOnly = TRUE))

####Loopeando cooo
################################################################################

setwd("/tables_expansion/ROCs")

all=readRDS(path[i,"path"])

results.AUCs=cbind(EFO_DOI[EFO_DOI[,"EFO"]%in%path[i,"EFO"],],"all")

colnames(results.AUCs)[ncol(results.AUCs)]="score"

####Listos para empezar con la tabla, primera posicion (all)




for (j in 1:nrow(results.AUCs)){
  
  
  ####Primero AUCs sin mas
  
  genes=unlist(strsplit(results.AUCs[j,"genes"],";"))
  
  temp.astro=cbind(all[	all[,"padj"]=="0",c("ENSG","page.rank")],0)
  
  colnames(temp.astro)[ncol(temp.astro)]="ROC"
  
  temp.astro[temp.astro[,"ENSG"]%in%genes,"ROC"]=1	
  
  results.AUCs[j,"count"]=sum(temp.astro[,"ENSG"]%in%genes)
  
  if(sum(temp.astro[,"ENSG"]%in%genes)>=10){
    
    results.AUCs[j,"AUC"]=roc(	as.numeric(temp.astro[,"ROC"]),
                               as.numeric(temp.astro[,"page.rank"]),
                               direction="<")$auc
  }else{
    
    results.AUCs[j,"AUC"]=NA    
    
  }
  
  
  #########NODES Randomization
  
  if(sum(temp.astro[,"ENSG"]%in%genes)>=10){
    
    ###Definimos las matrices
    
    
    matriz=all
    rownames(matriz)=matriz[,"ENSG"]
    
    Zscores=rep(0,20)
    
    for (k in 1:20){
      
      matriz[,"padj"]=sample(all[,"padj"])
      
      net.clean=graph_from_data_frame(d=string,vertices=matriz[,"ENSG"],directed=F)
      E(net.clean)$weight=rep(1,nrow(string))
      
      net.clean=igraph::simplify(net.clean,
                                 remove.loops = T,
                                 remove.multiple = T ,
                                 edge.attr.comb = c(weight="max","ignore"))
      
      
      E(net.clean)$weight=rep(1,nrow(string))
      matriz[,"page.rank"]=page_rank(net.clean, personalized=as.numeric(matriz[,"padj"]), weights=E(net.clean)$weight)$vector
      
      
      temp.scrmbl=cbind(matriz[as.numeric(matriz[,"padj"])==0,"page.rank"],0)
      colnames(temp.scrmbl)=c("pagerank","ROC")
      
      temp.scrmbl[rownames(temp.scrmbl)%in%genes,"ROC"]=1 
      
      if(sum(rownames(temp.scrmbl)%in%genes)>=10){
        
        Zscores[k]=roc(as.numeric(temp.scrmbl[,"ROC"]),
                       as.numeric(temp.scrmbl[,"pagerank"]),
                       direction="<")$auc
        
      }else{
        Zscores[k]=NA
      }
    }
    
    results.AUCs[j,"Zscore_nodes"]=   (as.numeric(results.AUCs[j,"AUC"])-mean(Zscores,na.rm=T))/ sd(Zscores,na.rm=T)
    
    
  }else{
    
    results.AUCs[j,"Zscore_nodes"]=NA    
    
  }
  
  #########TP RANDOMIZATION
  
  
  if(sum(temp.astro[,"ENSG"]%in%genes)>=10){
    
    
    
    A=temp.astro[temp.astro[,"ROC"]=="1",]
    B=temp.astro[temp.astro[,"ROC"]=="0",]
    
    A[,"ROC"]=0
    
    sel=c(rep(0,nrow(B)-nrow(A)),
          rep(1,nrow(A)))
    
    temp.TP=rep(0,20)
    
    for (k in 1:length(temp.TP)){
      
      B[,"ROC"]=sample(sel)
      
      temp.AUC=rbind(A,B)
      
      temp.TP[k]=roc(as.numeric(temp.AUC[,3]),
                     as.numeric(temp.AUC[,2]),
                     direction="<")$auc
      
    }
    
    results.AUCs[j,"Zscore_TP"]=(as.numeric(results.AUCs[j,"AUC"])-mean(as.numeric(temp.TP)))/sd(as.numeric(temp.TP))
    
  }else{
    
    results.AUCs[j,"Zscore_TP"]=NA   
    
  }}


saveRDS(results.AUCs,paste("/tables_expansion/ROCs/",path[i,"EFO"],".rds",sep=""))
