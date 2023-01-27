
#####Tabla con todas las traits juntas 


setwd("/tables_expansion/RDS_astro/")

path=paste(	"/tables_expansion/RDS_astro/",list.files(pattern= ".rds"),sep="")

all=readRDS(path[1])


for (i in 2:length(path)){
  
  all=rbind(all,readRDS(path[i]))
}


###Due to size constrains this file is not present in the repository
saveRDS(all,"/tables_expansion/all_together_sig.rds")
