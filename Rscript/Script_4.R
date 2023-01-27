
####Jaccard score calculations (genes from significant modules, figures 3 & 4)

jaccard.IB<-function(set1,set2){
  jac=sum(set1%in%set2)/sum(!duplicated(c(set1,set2)))
  return(jac)
}

######

setwd("/tables_expansion")

####This file is compressed using 7z (jaccard_KS.7z) and should be decompressed prior to run
jaccard.KS=readRDS("/tables_expansion/jaccard_KS.rds")
geneList=readRDS("/tables_expansion/genesList_KS.rds")

for(i in 1:nrow(jaccard.KS)){
  
  jaccard.KS[i,"jaccIndx"]=jaccard.IB(unlist(geneList[jaccard.KS[i,"A"]]),
                                      unlist(geneList[jaccard.KS[i,"B"]]))
  
}

saveRDS(jaccard.KS,"/tables_expansion/jaccard_KS_results.rds")

