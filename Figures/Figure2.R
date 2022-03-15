library(data.table)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(survival)
library(ggfortify)
library(survminer)




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
          "PDbyNTL"="#FCCDE5",
          "PD"="#D9D9D9"
)

quartcol <- c(brewer.pal(5,"Greens")[2:5])
names(quartcol) <- c("1","2","3","4")

batchcol <- rev(c(brewer.pal(3,"Reds")[2:3],brewer.pal(7,"Blues")[c(2,5,7)]))
names(batchcol) <- rev(c("On-Treatment-C","On-Treatment-A","Base-C","Base-B","Base-A"))

aux <- fread("data/SampleSourceData.txt",data.table = FALSE)
supp1 <- fread("TableS1.txt",data.table=FALSE) #download separately

supp1$id <- paste0(supp1$`Patient ID`,";",supp1$`Sample Timepoint`)
rownames(supp1) <- supp1$id
aux$`PD-1, % Positive Cells` <- supp1[aux$SampleID,"PD-1, % Positive Cells"]
aux$`CD8, % Positive Cells` <- supp1[aux$SampleID,"CD8, % Positive Cells"]
colnames(aux) <- gsub("[.]IHCquart","",colnames(aux))
aux$PFS <- aux$PFS.days / 7


lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}



pdf("Figure2_AB.pdf")
do.call(rbind,lapply(c("PD-1, % Positive Cells","CD8, % Positive Cells"),function(x) {
  temp <- aux
  colnames(temp) <- gsub(x,"Feature",colnames(temp))
  temp <- temp[!is.na(temp$SampleTimepoint) & !is.na(temp$Feature),]
  ttab <- do.call(rbind,lapply(c("Baseline","On-Treatment"),function(time) {
    temp$Group <- ifelse(temp$BestResponse=="PR",1,0)
    temp <- temp[temp$SampleTimepoint == time,]
    lmsum <- summary(lm(Feature ~ Cohort + Group,data=temp))
    lmfull <- (lm(Feature ~ Cohort + Group,data=temp))
    p <- lmsum$coefficients[length(rownames(lmsum$coefficients)),4]
    totalp <- lmp(lmfull)
    data.frame(Group=x,Test="BestResponse",Timepoint=time,Modelp=totalp,Testp=p,stringsAsFactors = FALSE)
  }))
  temp$BestResponse <- factor(temp$BestResponse,level=c("PR","SD","PD"))
  p <- ggplot(temp[temp$SampleTimepoint != "Progression",],aes(x=BestResponse,y=Feature)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=3,aes(color=Cohort),width = 0.2) +
    scale_color_manual(values=cols) +
    theme_bw() +
    labs(x="Best Response",y=x) +
    theme(text = element_text(size=20)) +
    facet_grid(. ~ SampleTimepoint,scales = "free",space="free")
  print(p)
  ttab
}))
dev.off()



pdf("Figure2_C.pdf")
do.call(rbind,lapply("PD-1, % Positive Cells",function(x) {
  temp <- aux
  colnames(temp) <- gsub(x,"Feature",colnames(temp))
  temp <- temp[!is.na(temp$SampleTimepoint) & !is.na(temp$Feature),]
  temp2 <- spread(temp[,c("PatientID","BestResponse","Cohort","SampleTimepoint","Feature")],SampleTimepoint,Feature)
  temp2$num <- rowSums(!is.na(temp2[,c("Baseline","On-Treatment")]))
  temp2 <- temp2[temp2$num == 2,]
  temp2$Change <- temp2$`On-Treatment` - temp2$Baseline
  ttab <- do.call(rbind,lapply("PRvNot",function(z) {
    temp2$Group <- ifelse(temp2$BestResponse=="PR",1,0)
    lmsum <- summary(lm(Change ~ Cohort + Group,data=temp2))
    lmfull <- (lm(Change ~ Cohort + Group,data=temp2))
    p <- lmsum$coefficients[length(rownames(lmsum$coefficients)),4]
    totalp <- lmp(lmfull)
    data.frame(Group=x,Test="BestResponse",Time="Subtract",Modelp=totalp,Testp=p,stringsAsFactors = FALSE)
  }))
  temp2$BestResponse <- factor(temp2$BestResponse,level=c("PR","SD","PD"))
  p <- ggplot(temp2,aes(x=BestResponse,y=Change)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=3,aes(color=Cohort),width = 0.2) +
    scale_color_manual(values=cols) +
    theme_bw() +
    labs(x="Best Response",y=paste0("Change in\n",x)) +
    theme(text = element_text(size=20))
  print(p)
  ttab
}))
dev.off()



pdf("Figure2_DE.pdf",width=15,height=8)
lapply(c("Baseline","On-Treatment"),function(time) {
  temp.1 <- aux[aux$SampleTimepoint==time & !is.na(aux$ImmuneGroup),c("PatientID","SampleTimepoint","BestResponse","Cohort","PFSGroup","PDL1","CD8","CD68","FOXP3","PD-1","T-cell","T-cell (CD8)","NK cell","B-cell","Macrophage/Monocyte","Myeloid dendritic cell","Neutrophil","Endothelial cell","Cancer-associated fibroblast","ImmuneGroup","OppositeImmuneGroup")]
  temp <- temp.1
  roworder <- c("T-cell (CD8)","T-cell","B-cell","Macrophage/Monocyte","Neutrophil","Myeloid dendritic cell","NK cell","Endothelial cell","Cancer-associated fibroblast")
  temp.1$BestResponse <- factor(temp.1$BestResponse,levels=c("PR","SD","PD"))
  temp <- temp[,colnames(temp) %in% roworder]
  temp <- t(scale(temp,center=TRUE,scale=TRUE))
  colnames(temp) <- temp.1$PatientID
  temp <- temp[roworder,]
  rownames(temp.1) <- temp.1$PatientID
  column_dend = as.dendrogram(hclust(dist(t(temp), method = "manhattan"), method = "ward.D"))
  column_dendlog.1 = hclust(dist(t(temp), method = "manhattan"), method = "ward.D")
  column_dend = dendextend::rotate(column_dend, temp.1$PatientID[order(temp.1$BestResponse,temp.1$PFSGroup)])
  bgroup <- if(time=="Baseline"){temp.1[,"ImmuneGroup"]}else{temp.1[,"OppositeImmuneGroup"]}
  ogroup <- if(time=="On-Treatment"){temp.1[,"ImmuneGroup"]}else{temp.1[,"OppositeImmuneGroup"]}
  ha = HeatmapAnnotation(show_annotation_name = T,
                         show_legend=T,
                         annotation_name_gp = gpar(fontsize = 14),
                         BestResponse=temp.1[,"BestResponse"],
                         Cohort=temp.1[,"Cohort"],
                         PD1=temp.1[,"PD-1"],
                         PDL1=temp.1[,"PDL1"],
                         FOXP3=temp.1[,"FOXP3"],
                         CD68=temp.1[,"CD68"],
                         CD8=temp.1[,"CD8"],
                         BaselineGroup=bgroup,
                         `On-TreatmentGroup`=ogroup,
                         col=list(BestResponse=cols,
                                  Cohort=cols,
                                  PD1=quartcol,
                                  PDL1=quartcol,
                                  FOXP3=quartcol,
                                  CD68=quartcol,
                                  CD8=quartcol,
                                  BaselineGroup=batchcol,
                                  `On-TreatmentGroup`=batchcol),
                         gp = gpar(fontsize = 18, col="black"),
                         na_col="white",
                         gap = unit(c(0,2,0,0,0,0,2,0,0), "mm"))
  h=Heatmap(as.matrix(temp),
            heatmap_legend_param = list(title = "Zscore"),
            column_title = paste0(time),
            column_dend_height = unit(3, "cm"),
            column_dend_reorder =T,
            cluster_columns=column_dend,
            cluster_rows=FALSE,
            row_dend_reorder = FALSE,
            gap = unit(30, "mm"),
            row_names_gp = gpar(fontsize = 14),
            top_annotation = ha
  )
  print(h)
  cat(paste0(time,"\n"))
})
dev.off()



pdf("Figure2_FG.pdf")
lapply(c("Baseline","On-Treatment"),function(y) {
  temp <- aux[aux$SampleTimepoint==y & !is.na(aux$ImmuneGroup),]
  fit6 = survfit(formula=Surv(as.numeric(PFS),as.numeric(factor(PFSCensor)))~ImmuneGroup,data=temp)
  temp2 <- temp[temp$ImmuneGroup %in% c("Base-A","Base-C","On-Treatment-A","On-Treatment-C"),]
  sd <- survdiff(Surv(as.numeric(PFS), PFSCensor==1)~ImmuneGroup,data=temp2)
  p.val1 <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  ttab <- data.frame(first="A",second="C",p=p.val1,stringsAsFactors = FALSE)
  ttab$label <- paste0(ttab$first," vs ",ttab$second,": ",round(ttab$p,3))
  linklabel <- paste(ttab$label,collapse="\n")
  plot.data=fortify(fit6,surv.connect=TRUE)
  p <- ggplot(plot.data, aes(time, surv,color=strata)) +
    geom_step(size=2) +
    geom_point(data = subset(plot.data, n.censor > 0),size=3,shape=3,aes(color=strata)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values=batchcol) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    labs(color="Group",x="Time (Weeks)",y="PFS (%)",title=paste0("PFS with Immune Groups at ",y," Timepoint\n",linklabel)) +
    theme(legend.position=c(.77,.77),
          legend.key=element_blank(),
          legend.background=element_blank(),
          aspect.ratio=1.0,
          legend.text=element_text(size=12),
          axis.text.y=element_text(size=15),
          axis.text.x=element_text(size=15),
          axis.title=element_text(size=16)) +
    coord_cartesian(xlim=c(0,105),expand=FALSE)
  q <- ggsurvplot(fit6, data = temp, risk.table = TRUE,palette=as.character(batchcol[names(batchcol) %in% unique(temp$ImmuneGroup)]),break.x.by=25,xlim=c(0,105),ylab="PFS")
  print(p)
  print(q)
  cat(y,"\n")
})
dev.off()





