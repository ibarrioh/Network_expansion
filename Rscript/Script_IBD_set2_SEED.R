
##libraries

library(igraph)

###FUNCIONES Modificando ASTRO para que no se lie con la ausencia de GAD

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

astro<-function(node.gwas,edge.string,all.nodes){
  
  ##Diffusion
  net=graph_from_data_frame(d=edge.string,vertices=node.gwas,directed=F)
  
  E(net)$weight=as.numeric(as.character(edge.string[,"combined_score"]))
  
  net.clean=igraph::simplify(net,
                             remove.loops = T,
                             remove.multiple = T ,
                             edge.attr.comb = c(weight="max","ignore"))
  
  page.rank=page_rank(net.clean, personalized=as.numeric(node.gwas[,"padj"]), weights=E(net.clean)$weight)
  
  node.gwas=cbind(node.gwas,page.rank$vector)
  colnames(node.gwas)[ncol(node.gwas)]="page.rank"
  
  deg=igraph::degree(net.clean)
  
  ##Network filter
  
  node.filter=node.gwas[as.numeric(node.gwas[,"page.rank"])>quantile(as.numeric(node.gwas[,"page.rank"]))[4],]
  
  colnames(node.filter)[1]="ENSP"
  
  edge.filter=edge.string[	as.character(edge.string[,1])%in%node.filter[,"ENSP"] & 
                             as.character(edge.string[,2])%in%node.filter[,"ENSP"] ,]
  
  node.filter=node.filter[node.filter[,"ENSP"]%in%c(as.character(edge.filter[,1]),as.character(edge.filter[,2])),]
  
  
  edge.string=edge.filter[,1:2]
  
  net=graph_from_data_frame(d=edge.string,vertices=node.filter,directed=F)
  
  E(net)$weight=as.numeric(as.character(edge.filter[,"combined_score"]))
  
  net.clean=igraph::simplify(net,
                             remove.loops = T,
                             remove.multiple = T ,
                             edge.attr.comb = c(weight="max","ignore"))
  
  cwt=cluster_walktrap(	net.clean, 
                        weights = E(net.clean)$weight, 
                        steps = 6,
                        merges = TRUE, 
                        modularity = TRUE, 
                        membership = TRUE)
  
  degree=igraph::degree(net.clean)
  
  node.filter=cbind(node.filter,degree,cwt$membership,cwt$modularity)
  
  colnames(node.filter)[(ncol(node.filter)-2):ncol(node.filter)]=c("node.degree","cluster.walktrap","modularity.walktrap")
  
  clust=as.matrix(as.data.frame(table(cwt$membership)))
  
  ####Recluster
  
  if(sum(as.numeric(clust[,2])>=300)>0){
    
    temp=rbind(c("0","0"),clust[as.numeric(clust[,2])>=300,])
    
    for (i in 1:nrow(temp)){
      
      node.re=node.filter[node.filter[,"cluster.walktrap"]==temp[i,1],]	
      edge.re=edge.filter[	as.character(edge.filter[,1])%in%node.re[,"ENSP"] &
                             as.character(edge.filter[,2])%in%node.re[,"ENSP"],]
      node.re=node.re[	node.re[,"ENSP"]%in%c(as.character(edge.re[,1]),as.character(edge.re[,2])),]
      
      net=graph_from_data_frame(d=as.data.frame(edge.re[,1:2]),vertices=node.re,directed=F)
      
      E(net)$weight=as.numeric(as.character(edge.re[,"combined_score"]))
      
      net.re=igraph::simplify(	net,
                               remove.loops = T,
                               remove.multiple = T ,
                               edge.attr.comb = c(weight="max","ignore"))
      
      cwt.re=cluster_walktrap(net.re, 
                              weights = E(net.re)$weight, 
                              steps = 6,
                              merges = TRUE, 
                              modularity = TRUE, 
                              membership = TRUE)
      
      node.re[,"cluster.walktrap"]=paste(node.re[,"cluster.walktrap"],cwt.re$membership,sep=";")
      node.re[,"modularity.walktrap"]=cwt.re$modularity
      
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"cluster.walktrap"]=node.re[,"cluster.walktrap"]
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"modularity.walktrap"]=node.re[,"modularity.walktrap"]
      
    }
    
  }
  
  ####Re-Recluster
  
  clust=as.matrix(as.data.frame(table(node.filter[,"cluster.walktrap"])))
  
  if(sum(as.numeric(clust[,2])>=300)>0){
    
    temp=rbind(c("0","0"),clust[as.numeric(clust[,2])>=300,])
    
    for (i in 1:nrow(temp)){
      
      node.re=node.filter[node.filter[,"cluster.walktrap"]==temp[i,1],]	
      edge.re=edge.filter[	as.character(edge.filter[,1])%in%node.re[,"ENSP"] &
                             as.character(edge.filter[,2])%in%node.re[,"ENSP"],]
      node.re=node.re[	node.re[,"ENSP"]%in%c(as.character(edge.re[,1]),as.character(edge.re[,2])),]
      
      net=graph_from_data_frame(d=as.data.frame(edge.re[,1:2]),vertices=node.re,directed=F)
      
      E(net)$weight=as.numeric(as.character(edge.re[,"combined_score"]))
      
      net.re=igraph::simplify(	net,
                               remove.loops = T,
                               remove.multiple = T ,
                               edge.attr.comb = c(weight="max","ignore"))
      
      cwt.re=cluster_walktrap(net.re, 
                              weights = E(net.re)$weight, 
                              steps = 6,
                              merges = TRUE, 
                              modularity = TRUE, 
                              membership = TRUE)
      
      node.re[,"cluster.walktrap"]=paste(node.re[,"cluster.walktrap"],cwt.re$membership,sep=";")
      node.re[,"modularity.walktrap"]=cwt.re$modularity
      
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"cluster.walktrap"]=node.re[,"cluster.walktrap"]
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"modularity.walktrap"]=node.re[,"modularity.walktrap"]
      
    }
    
  }
  
  
  ####Re-Re-Recluster
  
  clust=as.matrix(as.data.frame(table(node.filter[,"cluster.walktrap"])))
  
  if(sum(as.numeric(clust[,2])>=300)>0){
    
    temp=rbind(c("0","0"),clust[as.numeric(clust[,2])>=300,])
    
    for (i in 1:nrow(temp)){
      
      node.re=node.filter[node.filter[,"cluster.walktrap"]==temp[i,1],]	
      edge.re=edge.filter[	as.character(edge.filter[,1])%in%node.re[,"ENSP"] &
                             as.character(edge.filter[,2])%in%node.re[,"ENSP"],]
      node.re=node.re[	node.re[,"ENSP"]%in%c(as.character(edge.re[,1]),as.character(edge.re[,2])),]
      
      net=graph_from_data_frame(d=as.data.frame(edge.re[,1:2]),vertices=node.re,directed=F)
      
      E(net)$weight=as.numeric(as.character(edge.re[,"combined_score"]))
      
      net.re=igraph::simplify(	net,
                               remove.loops = T,
                               remove.multiple = T ,
                               edge.attr.comb = c(weight="max","ignore"))
      
      cwt.re=cluster_walktrap(net.re, 
                              weights = E(net.re)$weight, 
                              steps = 6,
                              merges = TRUE, 
                              modularity = TRUE, 
                              membership = TRUE)
      
      node.re[,"cluster.walktrap"]=paste(node.re[,"cluster.walktrap"],cwt.re$membership,sep=";")
      node.re[,"modularity.walktrap"]=cwt.re$modularity
      
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"cluster.walktrap"]=node.re[,"cluster.walktrap"]
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"modularity.walktrap"]=node.re[,"modularity.walktrap"]
      
    }
    
  }
  
  if(all.nodes==F){
    
    return(node.filter)	
    
  }else{
    
    node.gwas=cbind(node.gwas,deg)
    colnames(node.gwas)[ncol(node.gwas)]="degree"
    
    temp=node.filter[,c("ENSP","cluster.walktrap")]
    
    node.gwas=as.matrix(merge(node.gwas,temp,by.x="ENSG",by.y="ENSP",all.x=T))
    
    return(node.gwas)
    
  }
}

########################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################

########################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################

###Para transferir a yoda

all.gene.gwas=	readRDS("/Network_expansion/tables_IBD/set2/all_gene_gwas_FILTER_set2.rds")
all.node.gwas=	readRDS("/tables_expansion/all_node_gwas_FILTER.rds")
string=		readRDS("/tables_expansion/Combined_STRINGv11_OTAR281119_FILTER.rds")


################################################################################################

###Rscript

dis=table(all.gene.gwas[,"disease"])

##Argumento para loopear

i <- as.numeric(commandArgs(trailingOnly = TRUE))

##########
##########
##########
##########NETWORK

node.gwas=anotation(all.node.gwas,all.gene.gwas,disease=names(dis)[i])
astro.all=astro(node.gwas,as.data.frame(string),all.nodes=T)

astro.all=cbind(astro.all,names(dis)[i])
colnames(astro.all)[ncol(astro.all)]="Trait"

rm(node.gwas)
gc()

#####
#####
#####
#####
#####
#####
#####1000Randomizacions 
#####


matriz=matrix(0,nrow(all.node.gwas),1000)

rownames(matriz)=all.node.gwas[,"ENSG"]


node.gwas.first=anotation(all.node.gwas,all.gene.gwas,disease=names(dis)[i])

weights=c(	rep(0,(nrow(node.gwas.first)-(sum(node.gwas.first[,"padj"]!="0")*2))),
           node.gwas.first[node.gwas.first[,"padj"]!="0","padj"])

for (j in 1:1000){
  
  node.gwas=node.gwas.first
  
  node.gwas[node.gwas.first[,"padj"]!="0","padj"]=0
  node.gwas[node.gwas.first[,"padj"]=="0","padj"]=sample(weights) 
  
  net.clean=graph_from_data_frame(d=string,vertices=node.gwas,directed=F)
  E(net.clean)$weight=rep(1,nrow(string))
  
  net.clean=igraph::simplify(net.clean,
                             remove.loops = T,
                             remove.multiple = T ,
                             edge.attr.comb = c(weight="max","ignore"))
  
  E(net.clean)$weight=rep(1,nrow(string))
  matriz[,j]=page_rank(net.clean, personalized=as.numeric(node.gwas[,"padj"]), weights=E(net.clean)$weight)$vector
  
}

###Lines to save randomizations if needed
#path=paste("matrix_node/",names(dis)[i],".rds",sep="")
#saveRDS(matriz,path)

############################################################################################################
#####Tables assembly
############################################################################################################
############################################################################################################

astro.all=cbind(astro.all, "")
colnames(astro.all)[ncol(astro.all)]="Zsco.page.rank.node"

astro.all=cbind(astro.all, "")
colnames(astro.all)[ncol(astro.all)]="rankingIte1000.node"


##Re ordenar la matrizz de las randomizaciones 


matriz=matriz[astro.all[,"ENSG"],]


zsco=matrix(0,nrow(matriz),5)
colnames(zsco)=c("ori","zsco","mean","sd","rank")
rownames(zsco)=rownames(matriz)

zsco[,"ori"]=as.numeric(astro.all[,"page.rank"])

for (j in 1:nrow(matriz)){
  
  zsco[j,"mean"]=mean(matriz[j,])
  zsco[j,"sd"]=sd(matriz[j,])
  
  vec=c(zsco[j,"ori"],matriz[j,2:ncol(matriz)])
  names(vec)=c("ori",paste("ite",1:999))
  vec=vec[order(vec)]
  
  zsco[j,"rank"]= grep("ori",names(vec)) 
}

zsco[,"zsco"]=(zsco[,"ori"]-zsco[,"mean"])/zsco[,"sd"]

astro.all[,"Zsco.page.rank.node"]=zsco[,"zsco"]
astro.all[,"rankingIte1000.node"]=zsco[,"rank"]


path=paste(	"/Network_expansion/tables_IBD/set2/RDS_astro/",
            "ZSCO.",
            names(dis)[i],
            ".rds",	sep="")

saveRDS(astro.all,path)
