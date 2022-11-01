l = list(v16_7min, v16_10min, v16_15min, v16_21min, v16_30min, v16_60min)
#####################
num_sites = c()
num_proteins = c()
num_classI = c()
l_sites = list()
l_classI = list()
l_proteins = list()
for (i in 1:length(l)){
  l_sites[[i]] = l[[i]] %>% distinct(EG.IntPIMID, .keep_all = T)
  l_sites[[i]] = l_sites[[i]][which(l_sites[[i]]$EG.IntPIMID %like% '80'),]
  l_classI[[i]] = l_sites[[i]][which(l_sites[[i]]$EG.PTMAssayProbability <= 0.75),]
  l_proteins[[i]] = l[[i]] %>% distinct(PG.ProteinAccessions, .keep_all = T)
  num_sites[i] = nrow(l_sites[[i]])
  num_classI[i] = nrow(l_classI[[i]])
  num_proteins[i] = nrow(l_proteins[[i]])
}
#####################
tmp = data.frame()
tmp1 = data.frame()
tmp2 = c()
tmp3 = list()
grad = c(7,10,15,21,30,60)
for (i in 1:length(l_sites)){
  tmp = l_sites[[i]]
  for (j in 1:grad[i]){
    tmp1 = tmp[which(tmp$EG.ApexRT < j),]
    tmp2[j] = nrow(tmp1)
  }
  tmp3[[i]] = tmp2
}

