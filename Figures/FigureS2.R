library(data.table)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(scales)



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
          "PD"="#D9D9D9"
)


aux <- fread("data/SampleSourceData.txt",data.table=FALSE)
supp1 <- fread("TableS1.txt",data.table=FALSE)


supp1$id <- paste0(supp1$`Patient ID`,";",supp1$`Sample Timepoint`)
rownames(supp1) <- supp1$id
aux$`PDL1, % Positive Cells` <- supp1[aux$SampleID,"PDL1, % Positive Cells"]
aux$`CD68, % Positive Cells` <- supp1[aux$SampleID,"CD68, % Positive Cells"]
aux$`FOXP3, % Positive Cells` <- supp1[aux$SampleID,"FOXP3, % Positive Cells"]
aux$`CD8, % Positive Cells` <- supp1[aux$SampleID,"CD8, % Positive Cells"]


lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}




pdf("FigureS2_ABC.pdf",width=10,height=8)
do.call(rbind,lapply(c("PDL1, % Positive Cells","CD68, % Positive Cells","FOXP3, % Positive Cells"),function(x) {
  temp <- aux[aux$SampleTimepoint != "Progression",]
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
  p <- ggplot(temp,aes(x=BestResponse,y=Feature)) +
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



pdf("FigureS2_DEFG.pdf")
do.call(rbind,lapply(c("PDL1, % Positive Cells","CD8, % Positive Cells","CD68, % Positive Cells","FOXP3, % Positive Cells"),function(x) {
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



cdtest1 <- cor.test(aux$`CD8, % Positive Cells`[!is.na(aux$`CD8, % Positive Cells`) & !is.na(aux$`T-cell (CD8)`)],aux$`T-cell (CD8)`[!is.na(aux$`CD8, % Positive Cells`) & !is.na(aux$`T-cell (CD8)`)],method="pearson")
cdtest2 <- cor.test(aux$`CD8, % Positive Cells`[!is.na(aux$`CD8, % Positive Cells`) & !is.na(aux$`T-cell (CD8)`)],aux$`T-cell (CD8)`[!is.na(aux$`CD8, % Positive Cells`) & !is.na(aux$`T-cell (CD8)`)],method="spearman")

mtest1 <- cor.test(aux$`CD68, % Positive Cells`[!is.na(aux$`CD68, % Positive Cells`) & !is.na(aux$`Macrophage/Monocyte`)],aux$`Macrophage/Monocyte`[!is.na(aux$`CD68, % Positive Cells`) & !is.na(aux$`Macrophage/Monocyte`)],method="pearson")
mtest2 <- cor.test(aux$`CD68, % Positive Cells`[!is.na(aux$`CD68, % Positive Cells`) & !is.na(aux$`Macrophage/Monocyte`)],aux$`Macrophage/Monocyte`[!is.na(aux$`CD68, % Positive Cells`) & !is.na(aux$`Macrophage/Monocyte`)],method="spearman")


pdf("FigureS2_HI.pdf",height=5,width=7)
ggplot(aux,aes(x=`CD8, % Positive Cells`,y=`T-cell (CD8)`)) +
  geom_point(size=3) +
  geom_smooth(method="lm",se=FALSE,color="black",linetype="dashed") +
  labs(x="% CD8+ Cells by IHC",y="CD8 T-cells by MCP-counter",title=paste0("Spearman: p-value = ",scientific(cdtest2$p.value,2),", rho = ",round(cdtest2$estimate,2),"\nPearson: p-value = ",scientific(cdtest1$p.value,2),", r = ",round(cdtest1$estimate,2))) +
  theme_bw() +
  theme(text=element_text(size=14))
ggplot(aux,aes(x=`CD68, % Positive Cells`,y=`Macrophage/Monocyte`)) +
  geom_point(size=3) +
  geom_smooth(method="lm",se=FALSE,color="black",linetype="dashed") +
  labs(x="% CD68+ Cells by IHC",y="Macrophage/Monocyte by MCP-counter",title=paste0("Spearman: p-value = ",scientific(mtest2$p.value,2),", rho = ",round(mtest2$estimate,2),"\nPearson: p-value = ",scientific(mtest1$p.value,2),", r = ",round(mtest1$estimate,2))) +
  theme_bw() +
  theme(text=element_text(size=14))
dev.off()






