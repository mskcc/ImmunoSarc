library(data.table)
library(msigdbr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ComplexHeatmap)
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
          "PD"="#D9D9D9",
          "Baseline"="#CCEBC5",
          "On-Treatment"="#FDB462"
)

cols2 <- c("purple","blue","red","gray50")
names(cols2) <- c("H/H","H/L","L/H","L/L")


aux <- fread("data/GSEASourceData.txt",data.table = FALSE)
aux2 <- fread("data/SampleSourceData.txt",data.table = FALSE)
aux3 <- fread("data/DEGSourceData.txt",data.table = FALSE)
supp1 <- fread("TableS1.txt",data.table=FALSE)

supp1$id <- paste0(supp1$`Patient ID`,";",supp1$`Sample Timepoint`)
aux2$PFS <- aux2$PFS.days / 7

all_gene_setsh = as.data.frame(msigdbr(species = "Homo sapiens"))
myHall = all_gene_setsh[all_gene_setsh$gs_cat=="H",]
myHall2 <- myHall[,c("gs_name","gene_symbol")]
colnames(myHall2) <- c("ont","gene")


pdf("Figure3_A.pdf",width=5,height=6)
temp <- aux3
temp2.1 <- myHall$gene_symbol[myHall$gs_name=="HALLMARK_HEDGEHOG_SIGNALING"]
temp$Hedgehog <- ifelse(temp$target_id %in% temp2.1,"Yes","No")
temp$label <- ifelse(temp$target_id=="PDCD1","PDCD1",ifelse(temp$target_id=="GLI1","GLI1",""))
temp$logq <- -log(temp$qval,2)
temp <- temp[(order(temp$Hedgehog)),]
p <- ggplot(temp,aes(x=b,y=logq,fill=Hedgehog)) +
  geom_point(size=5,alpha=0.7,shape=21,color="gray50") +
  geom_label_repel(aes(label=label),segment.color="black",label.size = NA,segment.size = 0.5,point.padding = 0.25,max.overlaps=10000,xlim=c(2.5,3)) +
  geom_hline(yintercept= -log(0.05,2),linetype="dashed") +
  scale_fill_manual(values=c("white","red")) +
  labs(y="-Log2(q)",x="beta") +
  theme_bw() +
  theme(text=element_text(size=18))
print(p)
dev.off()


pdf("Figure3_B.pdf",width=5,height=4)
temp <- aux
colnames(temp) <- gsub("[.]PRvNot","",colnames(temp))
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
dev.off()



pdf("Figure3_C.pdf",width=14.5,height=6.5)
lapply(c("Baseline","On-Treatment"),function(time) {
  temp.1 <- aux2[aux2$SampleTimepoint==time & aux2$SampleID %in% supp1$id[!is.na(supp1$RNA)],]
  temp.1$BestResponse <- factor(temp.1$BestResponse,levels=c("PR","SD","PD"))
  colord <- temp.1[order(temp.1$SixMonth),"PatientID"]
  rownames(temp.1) <- temp.1$PatientID
  temp.1 <- temp.1[colord,]
  temp <- temp.1
  temp <- temp[,c("Allograft Rejection","Apical Junction","Bile Acid Metabolism","Cholesterol Homeostasis","Complement","Epithelial Mesenchymal Transition","Fatty Acid Metabolism","Hedgehog Signaling","Il6 Jak Stat3 Signaling","Interferon Alpha Response","Kras Signaling Dn","Myc Targets V2","Myogenesis","Peroxisome","Xenobiotic Metabolism")]
  temp <- t(scale(temp,center=TRUE,scale=TRUE))
  colnames(temp) <- temp.1$PatientID
  column_dend = as.dendrogram(hclust(dist(t(temp), method = "manhattan"), method = "ward.D"))
  column_dend = dendextend::rotate(column_dend, temp.1$Pat[order(temp.1$BestResponse,temp.1$SixMonth)])
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
  h4=Heatmap(as.matrix(temp),
             show_heatmap_legend = FALSE,
             column_title = paste0(time),
             column_dend_height = unit(3, "cm"),
             column_dend_reorder =FALSE,
             cluster_columns=column_dend,
             cluster_rows=TRUE,
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


pdf("Figure3_DE.pdf")
lapply(c("Baseline","On-Treatment"),function(x) {
  clintemp <- aux2[aux2$SampleTimepoint==x,]
  clintemp$Call <- ifelse(clintemp$CD8HiLo=="High",ifelse(clintemp$HedgehogHiLo=="High","H/H","H/L"),ifelse(clintemp$HedgehogHiLo=="High","L/H","L/L"))
  fit6 = survfit(formula=Surv(as.numeric(PFS),as.numeric(factor(PFSCensor)))~Call,data=clintemp)
  plot.data=fortify(fit6,surv.connect=TRUE)
  p <- ggplot(plot.data, aes(time, surv,color=strata)) +
    geom_step(size=2) +
    geom_point(data = subset(plot.data, n.censor > 0),size=3,shape=3,aes(color=strata)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values=cols2) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    labs(color="Quartile",x="Time (Weeks)",y="PFS (%)") +
    theme(legend.position=c(.77,.77),
          legend.key=element_blank(),
          legend.background=element_blank(),
          aspect.ratio=1.0,
          legend.text=element_text(size=12),
          axis.text.y=element_text(size=15),
          axis.text.x=element_text(size=15),
          axis.title=element_text(size=16)) +
    coord_cartesian(xlim=c(0,105),expand=FALSE)
  q <- ggsurvplot(fit6, data = clintemp, risk.table = TRUE,break.x.by=15,xlim=c(0,105),ylab="PFS",palette = c("purple","blue","red","gray50"))
  print(p)
  print(q)
})
dev.off()




