
##libraries

library(igraph)

###Functions to run the network expansion method and the scoring of the modules 

anotation<-function(node.list,genes,disease){
  
  genes=genes[genes[,"disease"]%in%disease,]
  
  for (j in 1:nrow(node.list)){
    if(node.list[j,"ENSG"]%in%genes[,"gene"]){
      node.list[j,"padj"]=max(as.numeric(genes[genes[,"gene"]%in%node.list[j,"ENSG"],"padj"]))
    }else{
      node.list[j,"padj"]=0
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
  
  
  ####Re-Re-Re-Recluster
  
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




#########################################################################################################################
#########################################################################################################################
####Load the requested files 

all.gene.gwas=	readRDS("/nfs/leia/research/beltrao/ibarrioh/OTAR44/genetic_portal/score050_270120/all.gene.gwas_filter_GP.rds")
all.node.gwas=	readRDS("/nfs/leia/research/beltrao/ibarrioh/OTAR44/genetic_portal/score020_161219/all_node_gwas_FILTER.rds")
string=		readRDS("/nfs/leia/research/beltrao/ibarrioh/OTAR44/genetic_portal/score020_161219/Combined_STRINGv11_OTAR281119_FILTER.rds")

setwd("/nfs/leia/research/beltrao/ibarrioh/OTAR44/genetic_portal/score050_270120")


########

temp=all.gene.gwas[,c("gene","disease")]
temp=temp[!duplicated(temp),]

dis=table(temp[,"disease"])
dis=dis[dis>1]


##Argumento para loopear

i <- as.numeric(commandArgs(trailingOnly = TRUE))

########Primero el metodo Astro hasta el final


node.gwas=anotation(all.node.gwas,all.gene.gwas,disease=names(dis[i]))

temp.nodes=astro(node.gwas,as.data.frame(string),all.nodes=T)


############################################################################################
####Calculo de modulos significativos


node=cbind(temp.nodes,names(dis[i]))

#######################################################################################################

colnames(node)[7]="Trait"


node=cbind(	node,0,0,0,1,1)

colnames(node)[(ncol(node)-4):ncol(node)]=c("Selected.cluster","Selected.fisher","Selected.KS","padj.fisher","padj.KS")


clusters=as.matrix(as.data.frame(table(node[,"cluster.walktrap"])))
clusters=cbind(clusters,0,0,0)
colnames(clusters)=c("clust","all.nodes","gwa.nodes","fisher","KS")
clusters=clusters[as.numeric(clusters[,"all.nodes"])>=10,]

for (j in 1:nrow(clusters)){
  
  clusters[j,"gwa.nodes"]=sum(node[node[,"cluster.walktrap"]%in%clusters[j,"clust"],"padj"]!="0")
  
}

clusters=clusters[as.numeric(clusters[,"gwa.nodes"])>0,]

if(length(clusters)>5){
  ##KS	
  
  for (j in 1:nrow(clusters)){
    
    x=log10(as.numeric(node[!is.na(node[,"cluster.walktrap"]),"page.rank"]))
    y=log10(as.numeric(node[node[,"cluster.walktrap"]==clusters[j,"clust"],"page.rank"]))
    
    clusters[j,"KS"]=ks.test(x,y,alternative="greater")$p.value
    
  }
  
  ##Fisher
  
  for (j in 1:nrow(clusters)){
    
    matrix=cbind(	c(as.numeric(clusters[j,"gwa.nodes"]),sum(node[,"padj"]!="0")-as.numeric(clusters[j,"gwa.nodes"])),
                  c(as.numeric(clusters[j,"all.nodes"]),sum(!duplicated(node[,1]))-as.numeric(clusters[j,"all.nodes"])))
    
    
    clusters[j,"fisher"]=fisher.test(matrix, alternative="greater")$p.value
  }
  
  clusters[,"fisher"]=p.adjust(as.numeric(clusters[,"fisher"]),method="BH")
  clusters[,"KS"]=p.adjust(as.numeric(clusters[,"KS"]),method="BH")
  
  
  
  node[node[,"cluster.walktrap"]%in%clusters[,"clust"],"Selected.cluster"]=1
  node[node[,"cluster.walktrap"]%in%clusters[as.numeric(clusters[,"fisher"])<=0.05,"clust"],"Selected.fisher"]=1	
  node[node[,"cluster.walktrap"]%in%clusters[as.numeric(clusters[,"KS"])<=0.05,"clust"],"Selected.KS"]=1
  
  for (j in 1:nrow(clusters)){
    
    node[node[,"cluster.walktrap"]%in%clusters[j,"clust"],"padj.KS"]=clusters[j,"KS"]
    node[node[,"cluster.walktrap"]%in%clusters[j,"clust"],"padj.fisher"]=clusters[j,"fisher"]
    
  }
  
  
  
  
}else{	
  
  node[,c("Selected.cluster",
          "Selected.fisher",
          "Selected.KS",
          "padj.fisher",
          "padj.KS")]=NA
  
}



path=paste(	"/nfs/leia/research/beltrao/ibarrioh/OTAR44/genetic_portal/score050_270120/RDS_astro/",
            "nodes.",
            names(dis[i]),
            ".rds",	sep="")

saveRDS(node,path)
