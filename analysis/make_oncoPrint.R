# import list of annotated mutations to eventDataFrame

eventDataFrame <- read_csv("mutations.csv", col_names = TRUE)
eventDataFrame <- eventDataFrame[eventDataFrame$Include==1,]
eventDataFrame <- eventDataFrame[eventDataFrame$Func.refGene %in% c("exonic","splicing"),]
eventDataFrame <- eventDataFrame[!eventDataFrame$ExonicFunc.refGene %in% c("synonymous SNV"),]
eventDataFrame <- eventDataFrame[MyNumeric(eventDataFrame$esp6500si_all)<0.01,]
eventDataFrame <- eventDataFrame[MyNumeric(eventDataFrame$'1000g2015aug_all')<0.01,]
eventDataFrame <- eventDataFrame[MyNumeric(eventDataFrame$ExAC_ALL)<0.01,]
eventDataFrame <- eventDataFrame[!(eventDataFrame$Ref=="C" & eventDataFrame$Alt=="T" & GetCosmicNumber(eventDataFrame$cosmic68)<3 & eventDataFrame$AF<0.1),]
eventDataFrame <- eventDataFrame[!(eventDataFrame$Ref=="G" & eventDataFrame$Alt=="A" & GetCosmicNumber(eventDataFrame$cosmic68)<3 & eventDataFrame$AF<0.1),]

eventDataFrame <- eventDataFrame[, c("ExonicFunc.refGene","Chr","Start","End","Ref","Alt",
                                     "AF","Gene.refGene","study","patient","tissue_pathology","event_type","sampleID",
                                     "esp6500si_all","cosmic68","1000g2015aug_all","ExAC_ALL")]
colnames(eventDataFrame) <- c("ann","chr","pos","stop","ref","alt","AF","gene","study","patient","tissue_pathology",
                              "event_type","sampleID","esp6500si_all","COSMIC","normal_af","exac_all")

eventDataFrame <- eventDataFrame[!eventDataFrame$ann %in% c("synonymous SNV"),]
eventDataFrame$ann[eventDataFrame$ann %in% c(".")] <- "splicing"
eventDataFrame$ann[eventDataFrame$ann %in% c("nonsynonymous SNV")] <- "missense"
eventDataFrame$ann[eventDataFrame$ann %in% c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution")] <- "inframe_indel"
eventDataFrame$ann[eventDataFrame$ann %in% c("frameshift deletion","frameshift insertion","frameshift substitution")] <- "frameshift"
eventDataFrame$ann[eventDataFrame$ann %in% c("stopgain","stoploss","startgain","startloss")] <- "nonsense"
eventDataFrame <- eventDataFrame[eventDataFrame$AF >= 0.02,]

patients <- unique(eventDataFrame$patient)

# import list of genes to GenesPanel 
# in our case we have 45 genes:
GenesPanel <- c("ARID1A","GATA3","PTEN","ATM","KMT2D","RB1","AKT1","CDH1","TP53","MAP2K4","NCOR1","NF1","ERBB2","BRCA1","RUNX1","CHEK2","BAP1","PIK3CA","FBXW7","MAP3K1","PIK3R1","KMT2C","NOTCH1","SF3B1","PBRM1","PDGFRA","CCND3","ESR1","ARID1B","EGFR","BRAF","FGFR1","MYC","CDKN2A","FGFR2","CCND1","ERBB3","MDM2","TBX3","BRCA2","IGF1R","CBFB","SMAD4","STK11","CCNE1") 

geneMatrix <- matrix("",nrow=length(patients),ncol=length(GenesPanel))
rownames(geneMatrix) <- patients
colnames(geneMatrix) <- GenesPanel
ER <- rep("",length(patients));  Her2 <- rep("",length(patients));  Grade <- rep("",length(patients)); 

for (i in 1:nrow(geneMatrix))
{
  Ind <- which(eventDataFrame$patient==rownames(geneMatrix)[i])
  ER[i] <- eventDataFrame$ER[Ind][1]
  Her2[i] <- eventDataFrame$Her2[Ind][1]
  Grade[i] <- eventDataFrame$Grade[Ind][1]
  for (j in 1:ncol(geneMatrix))
  {
    Ind <- which(eventDataFrame$patient==rownames(geneMatrix)[i] & eventDataFrame$gene==colnames(geneMatrix)[j])
    if (length(Ind)==1)
    {
      geneMatrix[i,j] <- eventDataFrame$ann[Ind]
    }
    if (length(Ind)>1)
    {
      if (length(unique(eventDataFrame$sampleID[Ind]))==1)
      {  
        geneMatrix[i,j] <- "multi_hit"
      }
      else
      {
        geneMatrix[i,j] <- eventDataFrame$ann[Ind[1]]
      }
    }
  }
}


write.table(data.frame(Gene="gene",what="what", Positive="Positive",Positive_tot="Positive_tot", Negative="Negative",Negative_Tot="Negative_Tot",Pval="p.value"), 
            file="MutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=FALSE)
for (g in GenesPanel)
{
  Res <- data.frame(Gene=g,what="ER", Positive=sum(geneMatrix[which(ER=="Positive"),g]!=""),Positive_tot=length(which(ER=="Positive")), Negative=sum(geneMatrix[which(ER=="Negative"),g]!=""),Negative_Tot=length(which(ER=="Negative")))
  Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
  write.table(Res, file="MutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
  Res <- data.frame(Gene=g,what="Her2", Positive=sum(geneMatrix[which(Her2=="Positive"),g]!=""),Positive_tot=length(which(Her2=="Positive")), Negative=sum(geneMatrix[which(Her2=="Negative"),g]!=""),Negative_Tot=length(which(Her2=="Negative")))
  Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
  write.table(Res, file="MutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
}
Res <- data.frame(Gene="all genes",what="ER", Positive=sum(geneMatrix[which(ER=="Positive"),]!=""),Positive_tot=length(GenesPanel)*length(which(ER=="Positive")), Negative=sum(geneMatrix[which(ER=="Negative"),]!=""),Negative_Tot=length(GenesPanel)*length(which(ER=="Negative")))
Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
write.table(Res, file="MutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
Res <- data.frame(Gene="all genes",what="Her2", Positive=sum(geneMatrix[which(Her2=="Positive"),]!=""),Positive_tot=length(GenesPanel)*length(which(Her2=="Positive")), Negative=sum(geneMatrix[which(Her2=="Negative"),]!=""),Negative_Tot=length(GenesPanel)*length(which(Her2=="Negative")))
Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
write.table(Res, file="MutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)

write.table(data.frame(Gene="gene",what="what", Positive="Int/Low",Positive_tot="Int/Low_tot", Negative="High",Negative_Tot="High_Tot",Pval="p.value"), 
            file="MutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
for (g in GenesPanel)
{
  Res <- data.frame(Gene=g,what="Grade", Positive=sum(geneMatrix[which(Grade %in% c("Low","Intermediate","1","2")),g]!=""),Positive_tot=length(which(Grade %in% c("Low","Intermediate","1","2"))), Negative=sum(geneMatrix[which(Grade %in% c("High","3")),g]!=""),Negative_Tot=length(which(Grade %in% c("High","3"))))
  Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
  write.table(Res, file="MutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)
}
Res <- data.frame(Gene="all genes",what="Grade", Positive=sum(geneMatrix[which(Grade %in% c("Low","Intermediate","1","2")),]!=""),Positive_tot=length(GenesPanel)*length(which(Grade%in% c("Low","Intermediate","1","2"))), Negative=sum(geneMatrix[which(Grade %in% c("High","3")),]!=""),Negative_Tot=length(GenesPanel)*length(which(Grade %in% c("High","3"))))
Res$pValue <- fisher.test(matrix(c(Res$Positive,Res$Positive_tot,Res$Negative,Res$Negative_Tot),2,2))$p.value
write.table(Res, file="MutationsStatistics.csv", row.names=FALSE, col.names=FALSE,sep="\t", quote = FALSE, append=TRUE)


col = c(missense = "violet", multi_hit = "black", nonsense = "red", UTR_variant = "yellow", splice_site = "green", inframe_indel = "orange",frameshift="blue",splicing="yellow")
alter_fun <- list(
  missense=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["missense"], col=NA))
  },
  UTR_variant=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["UTR_variant"], col=NA))
  },
  splice_site=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["splice_site"], col=NA))
  },
  inframe_indel=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["inframe_indel"], col=NA))
  },
  frameshift=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["frameshift"], col=NA))
  },
  nonsense=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["nonsense"], col=NA))
  },
  multi_hit=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["multi_hit"], col=NA))
  },
  splicing=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["splicing"], col=NA))
  }
)



write.table(cbind(colSums(geneMatrix!="")[roworder],rep(nrow(geneMatrix),ncol(geneMatrix))), file="plots/GenesFreq.tsv", row.names=TRUE,sep="\t", quote = FALSE)

PastelOrange <- rgb(red=100/100, green=70/100, blue=28/100, alpha=1)
PastelRed <- rgb(red=100/100, green=41/100, blue=38/100, alpha=1)
PastelYellow <- rgb(red=99/100, green=99/100, blue=59/100, alpha=1)

PastelGreen1 <- rgb(red=218/255, green=241/255, blue=219/255, alpha=1)
PastelGreen2 <- rgb(red=143/255, green=214/255, blue=148/255, alpha=1)
PastelGreen3 <- rgb(red=60/255, green=165/255, blue=68/255, alpha=1)
PastelGreen4 <- rgb(red=14/255, green=37/255, blue=15/255, alpha=1)

PastelViolet1 <- rgb(red=229/255, green=204/255, blue=228/255, alpha=1)
PastelViolet2 <- rgb(red=177/255, green=102/255, blue=174/255, alpha=1)

PastelBlue1 <- rgb(red=22/255, green=232/255, blue=235/255, alpha=1)
PastelBlue2 <- rgb(red=126/255, green=164/255, blue=179/255, alpha=1)

DarkerTurquoise <- rgb(0/255,131/255,133/255)
Maroon <- rgb(128/255,0,0)

Grade[Grade %in% c("Intermediate","Low")] <- "Low/Int"
Grade[is.na(Grade)] <- "Unknown"
Grade[Grade %in% c("1","2")] <- "Low/Int"
Grade[Grade %in% c("3")] <- "High"


ha = HeatmapAnnotation(df = data.frame(ER=ER,Her2=Her2,Grade=Grade), col = list(
  ER = c("Positive" =  DarkerTurquoise, "Negative" = "paleturquoise","Unknown" = "grey"),
  Her2 = c("Positive" =  DarkerTurquoise, "Negative" = "paleturquoise","Unknown" = "grey"),
  Grade = c("High" = DarkerTurquoise,"Low/Int" =  "paleturquoise","Unknown" = "grey")),
  numbers=anno_text(patients,gp=gpar(fontsize=3)),
  annotation_height =1,
  annotation_width =1,
  annotation_name_gp = gpar(fontsize=8),
  gp = gpar(col = "black",lwd=0.05)
  )



library(ComplexHeatmap)
pdf("oncoPrint.pdf",height=5)
oncoPrint(t(geneMatrix),
          col=col,
          alter_fun=alter_fun,
          show_column_names=TRUE,
          remove_empty_columns = FALSE,
          remove_empty_rows = FALSE,
          row_names_gp = gpar(fontsize = 5),
          column_names_gp = gpar(fontsize = 3),
          show_pct=TRUE,
          column_order=NULL,
          row_order=NULL,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot()),
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE)),
          heatmap_legend_param = list(title = "Alterations"),
          bottom_annotation = ha,
          alter_fun_is_vectorized = FALSE,
          pct_gp = gpar(fontsize = 6)
)
dev.off()
