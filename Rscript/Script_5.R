

####Jaccard indexes of trait ancestries (EFO) to use as benchmark for trait to trait distances
####This script OVERWRITE the file with the correlations adding new columns and values, careful, consider changing the paths


jaccard.IB<-function(set1,set2){
  jac=sum(set1%in%set2)/sum(!duplicated(c(set1,set2)))
  return(jac)
}

setwd("/tables_expansion")


corr.traits.norm=readRDS("/tables_expansion/corr.traits.norm.EFOanot.rds")
efo=readRDS("/tables_expansion/EFO191219_table.rds")

corr.traits.norm=as.matrix(corr.traits.norm)
corr.traits.norm=cbind(corr.traits.norm,0)
colnames(corr.traits.norm)[ncol(corr.traits.norm)]="common_ancest_EFO_jaccard"



for (i in 1:nrow(corr.traits.norm)){
  
  corr.traits.norm[i,"common_ancest_EFO_jaccard"]=jaccard.IB(unlist(efo$ancestors[corr.traits.norm[i,"i"]]),
                                                             unlist(efo$ancestors[corr.traits.norm[i,"j"]]))
}


saveRDS(corr.traits.norm,"/tables_expansion/corr.traits.norm.EFOanotROC.rds")

#################################################################################################################################

corr.traits.zsco=readRDS("/tables_expansion/corr.traits.zsco.EFOanot.rds")
efo=readRDS("/tables_expansion/EFO191219_table.rds")

corr.traits.zsco=as.matrix(corr.traits.zsco)
corr.traits.zsco=cbind(corr.traits.zsco,0)
colnames(corr.traits.zsco)[ncol(corr.traits.zsco)]="common_ancest_EFO_jaccard"



for (i in 1:nrow(corr.traits.zsco)){
  
  corr.traits.zsco[i,"common_ancest_EFO_jaccard"]=jaccard.IB(unlist(efo$ancestors[corr.traits.zsco[i,"i"]]),
                                                             unlist(efo$ancestors[corr.traits.zsco[i,"j"]]))
}


saveRDS(corr.traits.zsco,"/tables_expansion/corr.traits.zsco.EFOanotROC.rds")


