library(stringi)
library(stringr)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

report = list(`7min_IM`,
              `10min_IM`,
              `15min_IM`,
              `21min_IM`,
              `30min_IM`,
              `60min_IM`)
grad = c('7','10','15','21','30','60')
tmp = data.frame()
result = list()
result1 = list()
result2 = list()


for (i in 1:length(report)){
  tmp = report[[i]]
  tmp$EG.IntPIMID = gsub("\\[.16\\]", "", tmp$EG.IntPIMID)
  tmp$EG.IntPIMID = gsub("\\[.42\\]", "", tmp$EG.IntPIMID)
  tmp$EG.IntPIMID = gsub("\\[.57\\]", "", tmp$EG.IntPIMID)
  tmp = tmp%>%distinct(EG.IntPIMID, .keep_all = T)
  tmp = tmp[is.na(tmp$FG.IonMobilityPeakWidth) == F,]
  tmp = tmp[which(tmp$EG.IntPIMID %like% '80'),]
  tmp$EG.IntPIMID = gsub("\\[.80\\]", "i", tmp$EG.IntPIMID)
  tmp = tmp[which(tmp$FG.Charge == 2),]
  #tmp = tmp[which(str_count(tmp$EG.IntPIMID,'i') == 1),]
  tmp = split(tmp, tmp$PEP.StrippedSequence)
  data1 = list()
  x = 0
  y = 0
  for (seq in names(tmp)){
    if (nrow(tmp[[seq]]) > 1){
      data1[[seq]] = tmp[[seq]]
      data1[[seq]] = data.frame('Stripped_Seq' = data1[[seq]]$PEP.StrippedSequence,
                                'Modif' = data1[[seq]]$EG.IntPIMID,
                                'RT' = data1[[seq]]$EG.ApexRT,
                                'Apex_IM' = data1[[seq]]$FG.ApexIonMobility,
                                'Width_IM' = data1[[seq]]$FG.IonMobilityPeakWidth,
                                'Start_RT' = data1[[seq]]$EG.StartRT,
                                'Peak_Width' = data1[[seq]]$EG.PeakWidth)
      data1[[seq]]$Modif = gsub("\\[.42\\]", "", data1[[seq]]$Modif)
      data1[[seq]]$Modif = gsub("\\[.16\\]", "", data1[[seq]]$Modif)
      data1[[seq]]$Modif = gsub("\\[.80\\]", "i", data1[[seq]]$Modif)
      data1[[seq]]$Modif = str_replace_all(data1[[seq]]$Modif, '_','')
      data1[[seq]]$Len = str_length(data1[[seq]]$Stripped_Seq)
      data1[[seq]]$pos1 = NA
      data1[[seq]]$pos2 = NA
      data1[[seq]]$pos3 = NA
      data1[[seq]]$rel1 = NA
      data1[[seq]]$rel2 = NA
      data1[[seq]]$rel3 = NA
      for (j in 1:nrow(data1[[seq]])){
        x = as.numeric(unlist(gregexpr('i',data1[[seq]]$Modif[j])))
        if(length (x) == 1){
          x = x - 1
          data1[[seq]]$pos1[j] = x
          data1[[seq]]$rel1[j] = data1[[seq]]$pos1[j]/data1[[seq]]$Len[j]
        }
        if(length(x) == 2){
          x[1] = x[1] - 1
          x[2] = x[2] - 2
          data1[[seq]]$pos1[j] = x[1]
          data1[[seq]]$pos2[j] = x[2]
          data1[[seq]]$rel1[j] = data1[[seq]]$pos1[j]/data1[[seq]]$Len[j]
          data1[[seq]]$rel2[j] = data1[[seq]]$pos2[j]/data1[[seq]]$Len[j]
        }
        if(length(x) == 3){
          x[1] = x[1] - 1
          x[2] = x[2] - 2
          x[3] = x[3] - 3
          data1[[seq]]$pos1[j] = x[1]
          data1[[seq]]$pos2[j] = x[2]
          data1[[seq]]$pos3[j] = x[3]
          data1[[seq]]$rel1[j] = data1[[seq]]$pos1[j]/data1[[seq]]$Len[j]
          data1[[seq]]$rel2[j] = data1[[seq]]$pos2[j]/data1[[seq]]$Len[j]
          data1[[seq]]$rel3[j] = data1[[seq]]$pos3[j]/data1[[seq]]$Len[j]

        }
      }
    }
  }
  df1 = list()
  for (seq1 in names(data1)){
    if(nrow(data1[[seq1]])==2){
      df1[[seq1]] = data1[[seq1]]
    }
  }
  tmp1 = c()
  tmp2 = c()
  tmp3 = c()
  tmp4 = c()
  df1_new = list()
  for (k in 1:length(df1)){
    df1_new[[k]] = df1[[k]][order(df1[[k]]$Start_RT, decreasing = F),]
    tmp1[k] = (df1_new[[k]]$RT[2] - df1_new[[k]]$RT[1])/(0.5*(sum(df1_new[[k]]$Peak_Width)))
    #df1_new[[k]] = df1_new[[k]][order(df1_new[[k]]$FG.ApexIonMobility, decreasing = F),]
    tmp2[k] = abs(df1_new[[k]]$Apex_IM[2] - df1_new[[k]]$Apex_IM[1])/(0.5*(sum(df1_new[[k]]$Width_IM)))
    tmp3[k] = df1[[k]]$Stripped_Seq[1]
    tmp4[k] = abs(diff(df1[[k]]$rel1))
  }

  d = data.frame('RT' = tmp1, 'IM' = tmp2, 'Seq' = tmp3, 'rel' = tmp4)
  #d = d[which(d$RT > 0),]
  #d = d[which(d$rel != 0),]
  d3 = d[which(d$IM > 1),]
  d1 = d[which(d$IM > 1),]
  d1$ID = 'up'
  d2 = d[which(d$RT > 1),]
  d2$ID = 'down'
  a = d
  a$Grad = paste0(grad[i],'min')
  result[[i]] = d
  result1[[i]] = d2
  result2[[i]] = a
}
