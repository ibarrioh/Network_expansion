
##libraries

library(igraph)

###Functions modifying the network expansion method
###There are two, "annotation" (meant as a way to prepare the required files for the method itself) and astro. 
###Astro requires 3 types of file, node.gwas (the GWAS hits, they are the output file from anotation function),edge.string (the interactome)
###all.nodes (the nodes of the interactome with ENSG and gene names, this is also used as template for the output)

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
