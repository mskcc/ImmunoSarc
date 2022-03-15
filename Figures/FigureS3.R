library(data.table)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(stringr)


cols <- c("Leiomyosarcoma"="#8DD3C7",
          "UPS/MFH/High Grade MFS"="#FFED6F",
          "Osteosarcoma"="#BEBADA",
          "Chondrosarcoma"="#FB8072",
          "Liposarcoma"="#80B1D3",
          "Other"="gray50",
          "Small Blue Round Cell"='#e6194b',
          "Vascular"='#3cb44b',
          "ASPS"='#9a6324',
          "PR"="mediumpurple3",
          "SD"="violet",
          "PD"="#D9D9D9",
          "Baseline"="#CCEBC5",
          "On-Treatment"="#FDB462"
)

aux <- fread("data/SampleSourceData.txt",data.table=FALSE)
aux2 <- fread("data/GSEASourceData.txt",data.table=FALSE)
supp1 <- fread("TableS1.txt",data.table=FALSE)

supp1$id <- paste0(supp1$`Patient ID`,";",supp1$`Sample Timepoint`)
aux$BestResponse <- factor(aux$BestResponse,levels=c("PR","SD","PD"))


pdf("FigureS3_A.pdf",width=8,height=5)
ggplot(aux[aux$SampleTimepoint!="Progression",],aes(x=BestResponse,y=PDCD1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.7,size=3,aes(color=Cohort),width=0.2) +
  scale_color_manual(values=cols) +
  theme_bw() +
  labs(x="Best Response",y="Expression of PD1") +
  theme(text=element_text(size=18)) +
  facet_grid(. ~ SampleTimepoint,space="free",scale="free")
dev.off()



pdt <- do.call(rbind,lapply(c("PR","SD","PD"),function(x) {
  temp3 <- do.call(rbind,lapply(c("Baseline","On-Treatment"),function(y) {
    temp <- aux[!is.na(aux$PDCD1) & aux$SampleTimepoint==y,]
    t <- t.test(temp$PDCD1[temp$BestResponse==x],temp$PDCD1[temp$BestResponse !=x])
    data.frame(Timepoint=y,Test=x,p=t$p.value,meanx=t$estimate[1],meany=t$estimate[2],stringsAsFactors = FALSE)
  }))
  temp3
}))



pdf("FigureS3_I.pdf",width=8,height=5)
ggplot(aux[aux$SampleTimepoint!="Progression",],aes(y=GLI1,x=BestResponse)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=3,aes(color=Cohort),width = 0.2) +
  scale_color_manual(values=cols) +
  theme_bw() +
  labs(x="BestResponse",y="GLI1") +
  theme(text = element_text(size=20)) +
  facet_grid(.~ SampleTimepoint,space="free",scale="free")
dev.off()
glit <- do.call(rbind,lapply(c("PR","SD","PD"),function(x) {
  temp3 <- do.call(rbind,lapply(c("Baseline","On-Treatment"),function(y) {
    temp <- aux
    temp <- temp[!is.na(temp$GLI1) & temp$SampleTimepoint==y,]
    t <- t.test(temp$GLI1[temp$BestResponse==x],temp$GLI1[temp$BestResponse !=x])
    data.frame(Timepoint=y,Test=x,p=t$p.value,meanx=t$estimate[1],meany=t$estimate[2],stringsAsFactors = FALSE)
  }))
  temp3
}))





pdf("FigureS3_H.pdf",width=14.5,height=6.5)
group <- c(strsplit("CDK5R1/TLE3/L1CAM/HEY1/PTCH1/THY1/VLDLR/SCG2/CELSR1/HEY2/ACHE/UNC5C/GLI1/CNTFR/NKX6-1/AMOT/CRMP1","/")[[1]])
temp.1 <- aux[aux$SampleID %in% supp1$id[!is.na(supp1$RNA)],]
lapply(c("Baseline","On-Treatment"),function(time) {
  temp.1 <- temp.1[temp.1$SampleTimepoint==time,]
  temp.1$BestResponse <- factor(temp.1$BestResponse,levels=c("PR","SD","PD"))
  temp.1$Pat <- paste0(temp.1$PatientID)
  colord <- temp.1[order(temp.1$BestResponse),"Pat"]
  rownames(temp.1) <- temp.1$Pat
  temp.1 <- temp.1[colord,]
  temp <- temp.1
  temp <- temp[,as.character(unique(group))]
  temp <- t(scale(temp,center=TRUE,scale=TRUE))
  colnames(temp) <- temp.1$Pat
  col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
  ha2 = HeatmapAnnotation(show_annotation_name = F,
                          show_legend=F,
                          annotation_name_gp = gpar(fontsize = 12),
                          BestResponse=temp.1[,"BestResponse"],
                          Cohort=temp.1[,"Cohort"],
                          TimePoint=temp.1[,"SampleTimepoint"],
                          col=list(BestResponse=cols,
                                   Cohort=cols,
                                   TimePoint=cols
                          ),
                          gp = gpar(fontsize = 18, col="black"),
                          na_col="white")
  h4=Heatmap(as.matrix((temp)),
             show_heatmap_legend = FALSE,
             column_title = paste0(time),
             column_dend_height = unit(3, "cm"),
             column_dend_reorder =FALSE,
             cluster_columns=FALSE,
             cluster_rows=TRUE,
             col=col_fun,
             row_dend_reorder = FALSE,
             show_row_dend = FALSE,
             gap = unit(2, "mm"),
             top_annotation = ha2
  )
  if(time=="Baseline") {hb <<- h4} else {hc <<- h4}
  cat(paste0(time,"\n"))
})
ht_list = hb + hc
draw(ht_list, column_title_gp = gpar(fontsize = 16))
dev.off()





pdf("FigureS3_DEFG.pdf",width=8,height=8)
temptab <- do.call(rbind,lapply(c("T-cell (CD8)","Hedgehog Signaling"),function(y) {
  temp <- aux
  colnames(temp) <- gsub("Cohort","Group",gsub(gsub("[(]","[(]",gsub("[)]","[)]",y)),"Test",colnames(temp)))
  temptab2 <- do.call(rbind,lapply(c("Baseline","On-Treatment"),function(time) {
    temp <- temp[temp$SampleTimepoint == time & !is.na(temp$SampleTimepoint) & !is.na(temp$Group) & !is.na(temp$Test),]
    tests <- unique(temp$Group[!(temp$Group %in% c("Other"))])
    testcount <- plyr::count(temp$Group)
    tests <- tests[tests %in% testcount$x[testcount$freq >= 3]]
    temptab3 <- do.call(rbind,lapply(tests,function(piece) {
      temp$Piece <- ifelse(temp$Group==piece,1,0)
      t <- t.test(temp$Test[temp$Piece==0],temp$Test[temp$Piece==1])
      data.frame(Piece=piece,Test=y,Time=time,NumTrue=length(temp$Test[temp$Piece==1]),NumFalse=length(temp$Test[temp$Piece==0]),MeanTrue=t$estimate[2],MeanFalse=t$estimate[1],p=t$p.value,stringsAsFactors = FALSE)
    }))
    temptab3$label <- paste0(temptab3$Piece,": ",round(temptab3$p,3))
    title <- paste0("Time: ",time,"\n",paste(temptab3$label,collapse="\n"))
    colnames(temp) <- gsub("BestResponse","color",colnames(temp)) 
    temp$Group <- factor(temp$Group,levels=c("PR","SD","PD","Liposarcoma","Chondrosarcoma","Osteosarcoma","UPS/MFH/High Grade MFS","Leiomyosarcoma","Vascular","Small Blue Round Cell","ASPS","Other"))
    p <- ggplot(temp,aes(x=as.factor(Group),y=Test)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width=0.25,size=5,alpha=0.7,aes(color=color)) +
      labs(x="Cohort",y=y,title=title) +
      theme_bw() +
      scale_color_manual(values=cols,na.value="gray50") +
      theme(text = element_text(size=18),
            title = element_text(size=12),
            axis.title = element_text(size=18),
            axis.text.x = element_text(angle=45,hjust=1))
    print(p)
    temptab3
  }))
  temptab2
}))
dev.off()


pdf("FigureS3_BC.pdf",width=5,height=4)
lapply(c("FullLM","PDvNot"),function(x) {
  temp <- aux2
  colnames(temp) <- gsub(paste0(".",x),"",colnames(temp))
  temp <- temp[temp$p.adjust < 0.05,]
  temp <- temp[order(temp$p.adjust,temp$pvalue),]
  use <- c(temp$ID[temp$NES<0 & temp$p.adjust < 0.05][1:10],temp$ID[temp$NES>0 & temp$p.adjust < 0.05][1:10])
  top <- temp[temp$ID %in% use,]
  top <- top[(order(top$NES)),]
  top$Term <- str_to_title(gsub("_"," ",gsub("HALLMARK_","",top$ID)),locale="en")
  top$Term <- ifelse(top$Term=="Il6 Jak Stat3 Signaling","IL6/JAK/STAT3 Signaling",top$Term)
  top$Term <- ifelse(top$Term=="E2f Targets","E2F Targets",top$Term)
  top$Term <- ifelse(top$Term=="G2m Checkpoint","G2M Checkpoint",top$Term)
  top$Term <- factor(top$Term,levels=unique(top$Term))
  p <- ggplot(top,aes(y=NES,x=Term)) +
    geom_bar(stat="identity") +
    labs(fill="AdjustedP < 0.05") +
    theme_bw() +
    scale_x_discrete(labels=function(x) str_wrap(x, width = 40)) +
    coord_flip()
  print(p)
  cat(1)
})
dev.off()




