##Script to Identify Differentially Expressed Genes and Gene Set Enrichment Analysis
##Relevant to Figure 3 and supplemental figure 2
##Author - Allison L Richards


cols <- c("Leiomyosarcoma (LMS)"="#8DD3C7",
          "UPS/MFH/High Grade MFS"="#FFED6F",
          "Osteosarcoma"="#BEBADA",
          "Chondrosarcoma"="#FB8072",
          "Dedifferentiated/Pleomorphic LPS"="#80B1D3",
          "Liposarcoma (LPS)"="#80B1D3",
          "Small Blue Round Cell (SBRC)"='#e6194b',
          "Vascular"='#3cb44b',
          "ASPS"='#9a6324',
          "Other"="gray50",
          "PR"="mediumpurple3",
          "SD"="violet",
          "PD"="#D9D9D9",
          "PDbyNTL"="#D9D9D9",
          "Baseline"="#CCEBC5",
          "C2"="#FDB462",
          "Progression"="#D95F02"
)


library(data.table)
library(dplyr)
library(tidyr)
library(sleuth)
library(fgsea)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(msigdbr)
library(GSVA)
library(scales)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(ggrepel)
library(glmnet)
library(nnet)
library(survminer)
library(lubridate)
library(survival)
library(ggfortify)


clin <- fread("XXX.txt",data.table=FALSE)

kall #List of RNAseq abundance files from dbGap
kall2 <- sapply(strsplit(kall,"/"),function(x) {x[grep("abundance",x)-1]})


#For overall differences between PR vs Not
des <- do.call(rbind,lapply(rownames(clin),function(x) {
  tab <- clin[x,]
  if(tab$SampleID %in% kall2) {
    tab$path <- gsub("/abundance.tsv","",kall[grepl(tab$SampleID,kall)])
    colnames(tab) <- gsub("SampleID","sample",colnames(tab))
    tab$PRvNot <- ifelse(tab$BestResponse=="PR",1,0)
    tab
  }
}))


##For conversion analysis
##Add in the group names and the group names for the sample from the opposite time point
#B1: Cold Group at baseline
temp <- des[!is.na(des$OppGroup) & des$Group=="B1",]
des2 <- des[des$`CMO Patient ID` %in% temp$`CMO Patient ID` & des$`Sample Timepoint`=="Baseline",]


#Use "des" or "des2" to run sleuth via Sleuth_Prep_Run.R
#bsub -n 8 -R rusage[mem=40] -W 1:59 'Rscript Sleuth_Prep_Run.R -d des.txt -t Gene -o Sleuth_Prep_Result.Rdata'
#bsub -n 8 -R rusage[mem=40] -W 1:59 'Rscript Sleuth_Prep_Run.R -d des2.txt -t Gene -o Sleuth_Prep_Result_Conversion.Rdata'



#Read in sleuth_prep results
mods <- list.files("Directory",pattern="Sleuth_Prep_Result",full.names = TRUE)


runthings <- function(x) {
  loaded <- sleuth_load(x)
  if(!grepl("Conversion",x)) {
    loaded <- sleuth_fit(loaded, ~`Sample Timepoint` + `Patient ID#` + Cohort + PRvNot + Purity, 'coh_PRvNot')
  } else {
    loaded <- sleuth_fit(loaded, ~Cohort + OppGroup, 'coh_OppGroup')
  }
  use <- names(models(loaded))
  if(length(use) > 0) {
    testcol <- sapply(strsplit(use,"_"),function(x){x[2]})
    fulltable <- data.frame(name=use,test=testcol,stringsAsFactors = FALSE)
    fulltab <- do.call(rbind,lapply(fulltable$name,function(z) {
      loaded <- sleuth_wt(loaded,fulltable$test[fulltable$name==z],z)
      tab <- sleuth_results(loaded, test=fulltable$test[fulltable$name==z], which_model=z, show_all = FALSE)
      cat(paste0(z," = ",length(tab$target_id[tab$qval < 0.05])," -> ",length(tab$target_id[tab$qval < 0.3])),sep="\n")
      tab$test <- fulltable$test[fulltable$name==z]
      tab$model <- z
      tab
    }))
    fulltab
  }
}

largetab <- do.call(rbind,lapply(mods,runthings))
##Table with all genes, differentially expressed have q-value < 0.05



##Gene Set Enrichment Analysis

#Let's use mutsig to be consistent
all_gene_setsh = as.data.frame(msigdbr(species = "Homo sapiens"))
myHall = all_gene_setsh[all_gene_setsh$gs_cat=="H",]
myHall2 <- myHall[,c("gs_name","gene_symbol")]
colnames(myHall2) <- c("ont","gene")


test <- unique(largetab$model)

#Relevant for Figure S3G
pdf("GSEAPlots.pdf")
ori <- do.call(rbind,lapply(test,function(x) {
  temp <- largetab[largetab$model == x,]
  temp <- temp[order(temp$qval),]
  temp <- temp[!duplicated(temp$target_id),]
  rank <- temp$b[rev(order(temp$b))]
  names(rank) <- temp$target_id[rev(order(temp$b))]
  edo2 <- GSEA(rank, TERM2GENE=myHall2, verbose=FALSE,by='fgsea',minGSSize = 15,maxGSSize = 600,nPerm = 10000,pvalueCutoff = 1)
  saveRDS(edo2, file = paste0(x,".rds"))
  set <- as.data.frame(edo2)
  set <- set[,colnames(set)!="core_enrichment"]
  set <- set[set$p.adjust < 0.05,]
  p <- gseaplot2(edo2, geneSetID = 1:5,title=paste0("Model: ",x,"\nTop 5 Pathways"))
  print(p)
  fgseaRes <- as.data.frame(edo2)
  fgseaRes$model <- x
  fgseaRes
}))
dev.off()



#Used the output of this to choose pathways for Figure 3B and 3C

pdf("Figure3B.pdf",width=5,height=4)
temp <- ori[ori$model=="coh_PRvNot" & ori$p.adjust < 0.05,]
temp <- temp[order(temp$p.adjust),]
use <- c(temp$ID[temp$model=="coh_PRvNot" & temp$NES<0 & temp$p.adjust < 0.05][1:10],temp$ID[temp$model=="coh_PRvNot" & temp$NES>0 & temp$p.adjust < 0.05][1:10]) #Note: There are only 6 significantly down-regulated pathways
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



#To create Figure 3C
full <- sleuth_load("Sleuth_Prep_Result.Rdata")
mat <- as.data.frame(sleuth_to_matrix(full,"obs_norm","tpm"))


mat2 <- if(TRUE) {
  use <- c(ori$ID[ori$model=="coh_PRvNot" & ori$NES<0 & ori$p.adjust < 0.05][1:10],ori$ID[ori$model=="coh_PRvNot" & ori$NES>0 & ori$p.adjust < 0.05][1:10])
  use2 <- ori[ori$ID %in% use & ori$model=="coh_PRvNot",colnames(ori) %in% c("ID","core_enrichment")]
  group <- do.call(rbind,lapply(unique(use2$ID),function(x) {
    data.frame(Pathway=x,Symbol=strsplit(use2$core_enrichment[use2$ID==x],split="/")[[1]],stringsAsFactors = FALSE)
  }))
  group$Pathway <- str_to_title(gsub("_"," ",gsub("HALLMARK_","",group$Pathway)),locale="en")
  group2 <- split(group$Symbol,group$Pathway)
  mattemp <- mat[,colnames(mat) %in% des$sample[des$`Sample Timepoint`!="Progression"]]
  test2all <- gsva(as.matrix(mattemp),group2,method="ssgsea")
  as.data.frame(t(test2all))
}
colnames(mat2) <- str_to_title(gsub("_"," ",gsub("HALLMARK_","",colnames(mat2))),locale="en")
mat2$sample <- rownames(mat2)
temp <- merge(des,mat2,by="sample")


pdf("Figure3C.pdf",width=14.5,height=6.5)
lapply(c("Baseline","C2"),function(time) {
  temp.1 <- temp[temp$`Sample Timepoint`==time,]
  temp.1$BestResponse <- factor(temp.1$BestResponse,levels=c("PR","SD","PDbyNTL","PD"))
  temp.1$Pat <- paste0(temp.1$`CMO Patient ID`)
  colord <- temp.1[order(temp.1$BestResponse),"Pat"]
  rownames(temp.1) <- temp.1$Pat
  temp.1 <- temp.1[colord,]
  temp <- temp.1
  temp <- temp[,as.character(unique(group$Pathway))]
  temp <- t(scale(temp,center=TRUE,scale=TRUE))
  colnames(temp) <- temp.1$Pat
  column_dend = as.dendrogram(hclust(dist(t(temp), method = "manhattan"), method = "ward.D"))
  column_dend = dendextend::rotate(column_dend, temp.1$Pat[order(temp.1$BestResponse)])
  ha2 = HeatmapAnnotation(show_annotation_name = F,
                          show_legend=F,
                          annotation_name_gp = gpar(fontsize = 12),
                          BestResponse=temp.1[,"BestResponse"],
                          Cohort=temp.1[,"Cohort"],
                          col=list(BestResponse=cols,
                                   Cohort=cols
                          ),
                          gp = gpar(fontsize = 18, col="black"),
                          na_col="white")
  h4=Heatmap(as.matrix((temp)),
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
print(hb + hc)
dev.off()



pdf("Figure3A.pdf",width=5,height=6)
lapply(c("coh_PRvNot"),function(x) {
  temp <- largetab[largetab$model ==x,]
  temp2.1 <- myHall$gene_symbol[myHall$gs_name=="HALLMARK_HEDGEHOG_SIGNALING"]
  temp$shh <- ifelse(temp$target_id %in% temp2.1,"Yes","No")
  temp$label <- ifelse(temp$target_id=="PDCD1","PDCD1",ifelse(temp$target_id=="GLI1","GLI1",""))
  temp$logq <- -log(temp$qval,2)
  temp <- temp[(order(temp$shh)),]
  p <- ggplot(temp,aes(x=b,y=logq,fill=shh)) +
    geom_point(size=5,alpha=0.7,shape=21,color="gray50") +
    geom_label_repel(aes(label=label),segment.color="black",label.size = NA,segment.size = 0.5,point.padding = 0.25,max.overlaps=10000,xlim=c(2.5,3),ylim=c(5,7.5)) +
    geom_hline(yintercept= -log(0.05,2),linetype="dashed") +
    scale_fill_manual(values=c("white","red")) +
    labs(title=paste0(x),y="-Log2(q)",x="beta") +
    theme_bw() +
    theme(text=element_text(size=18))
  print(p)
  cat(1)
})
dev.off()




##Modeling of Immune Populations and ssGSEA values for top pathways against PFS
##Calculate ssGSEA scores for all significant pathways (previously only did top 10 up and 6 down)
gsstuff <- do.call(rbind,lapply(c("All"),function(time) {
  mat2 <- if(TRUE) {
    group <- myHall2[myHall2$ont %in% ori$ID[ori$model=="coh_PRvNot" & ori$p.adjust < 0.05],]
    colnames(group) <- c("Pathway","Symbol")
    group$Pathway <- str_to_title(gsub("_"," ",gsub("HALLMARK_","",group$Pathway)),locale="en")
    group2 <- split(group$Symbol,group$Pathway)
    mattemp <- mat[,colnames(mat) %in% des$sample[des$`Sample Timepoint`!="Progression"]]
    test2all <- gsva(as.matrix(mattemp),group2,method="ssgsea")
    as.data.frame(t(test2all))
  }
  colnames(mat2) <- str_to_title(gsub("_"," ",gsub("HALLMARK_","",colnames(mat2))),locale="en")
  mat2$sample <- rownames(mat2)
  mat3 <- gather(mat2,Pathway,gsva,c(1:length(colnames(mat2))-1))
  mat3$time <- time
  mat3
}))



b <- unique(clin$SampleID[clin$`Sample Timepoint`=="Baseline"])
c <- unique(clin$SampleID[clin$`Sample Timepoint`=="C2"])
gsstuff$timepoint <- ifelse(gsstuff$sample %in% b,"Baseline",ifelse(gsstuff$sample %in% c,"C2","UNKNOWN"))


gs2.2 <- do.call(rbind,lapply(unique(gsstuff$timepoint),function(time) {
  paths <- unique(gsstuff$Pathway[gsstuff$timepoint==time & gsstuff$time=="All"])
  temp <- do.call(rbind,lapply(paths,function(x) {
    temp2 <- gsstuff[gsstuff$timepoint==time & gsstuff$time=="All" & gsstuff$Pathway==x,]
    twen <- quantile(temp2$gsva,0.25,na.rm=TRUE)
    fif <- quantile(temp2$gsva,0.5,na.rm=TRUE)
    sev <- quantile(temp2$gsva,0.75,na.rm=TRUE)
    temp2$quart <- ifelse(temp2$gsva < twen,1,ifelse(temp2$gsva < fif,2,ifelse(temp2$gsva < sev,3,4)))
    temp2$min <- min(temp2$gsva)
    temp2$max <- max(temp2$gsva)
    temp2
  }))
  temp
}))



df.2 <- do.call(rbind,lapply(c("Baseline","C2"),function(time) {
  gstemp <- spread(gs2.2[gs2.2$timepoint==time,c("sample","Pathway","gsva")],Pathway,gsva)
  clintemp <- merge(clin,gstemp,by.x="SampleID",by.y="sample")
  temp <- do.call(rbind,lapply(c("All"),function(x) {
    names <- colnames(clintemp)[grepl(paste(c("T.cell.CD8","T.cell","B.cell","Macrophage.Monocyte","Neutrophil","Myeloid.dendritic.cell","NK.cell","Endothelial.cell","Cancer.associated.fibroblast","cytotoxicity.score","Monocyte",colnames(gstemp)),collapse="|"),colnames(clintemp))]
    XP <- data.matrix(clintemp[,names])
    temp2 <- do.call(rbind,lapply(c("PFS"),function(y) {
      clintemp$PFS <- clintemp$`Progression-Free Survival (PFS) in Days` / 7
      YP <- Surv(as.numeric(clintemp$PFS), clintemp$PFSReached)
      lasso <- cv.glmnet(x = XP,y = YP,family = 'cox')
      n <- as.data.frame(as.matrix(coef(lasso, s=lasso$lambda.min)))
      n$Test <- rownames(n)
      n$Test <- ifelse(grepl(" ",n$Test),paste0("`",n$Test,"`"),n$Test)
      if(length(rownames(n)[n$`1`!=0 & rownames(n) != "(Intercept)"]) >0) {
        inmod <- paste(paste0("`",rownames(n)[n$`1`!=0 & rownames(n) != "(Intercept)"],"`"),collapse=" + ")
        clintemp$su <- Surv(as.numeric(clintemp$PFS), clintemp$PFSReached==1)
        mod <- coxph(as.formula(paste0("su ~",inmod)),data=clintemp)
        mod2 <- summary(mod)
        co <- as.data.frame(coef(mod2))
        n$indp <- co[n$Test,"Pr(>|z|)"]
        co["T.cell.CD8.mcpraw","Pr(>|z|)"]
        n$overallp <- mod2$logtest[3]
      } else {
        n$indp <- NA
        n$overallp <- NA
      }
      n$time <- time
      n$options <- x
      n$y <- y
      n
    }))
    temp2
  }))
  temp
}))


#df.2 shows that only CD8 T cells and Shh signaling are predictive of PFS in Baseline and On-Treatment (C2) samples
cols2 <- c("purple","blue","red","gray50")
names(cols2) <- c("H/H","H/L","L/H","L/L")

co <- plyr::count(df.2[df.2$`1`!=0,],vars=c("options","time"))

pdf("Figure3D.pdf")
lapply(rownames(co),function(x) {
  palist <- gsub("`","",df.2$Test[df.2$y=="PFS" & df.2$time==co[x,"time"] & df.2$options==co[x,"options"] & df.2$`1`!=0])
  gstemp <- gs2.2[gs2.2$timepoint == co[x,"time"] & gs2.2$Pathway %in% palist,]
  gstemp$hl <- ifelse(gstemp$quart %in% c(1,2),"Low","High")
  gstemp <- spread(gstemp[,c("sample","Pathway","hl")],Pathway,hl)
  rownames(gstemp) <- gstemp$sample
    clintemp <- merge(clin[,c("SampleID","Progression-Free Survival (PFS) in Days","PFSReached")],gstemp,by.x="SampleID",by.y="sample")
    temp <- clin[clin$SampleID %in% clintemp$SampleID & !is.na(clin$SampleID),c("SampleID",palist[palist %in% colnames(clin)])]
    thing <- colnames(temp)[2]
    colnames(temp) <- c("sample","thing")
    twen <- quantile(temp$thing,0.25,na.rm=TRUE)
    fif <- quantile(temp$thing,0.5,na.rm=TRUE)
    sev <- quantile(temp$thing,0.75,na.rm=TRUE)
    temp$quart <- ifelse(temp$thing < twen,1,ifelse(temp$thing < fif,2,ifelse(temp$thing < sev,3,4)))
    temp$hl <- ifelse(temp$quart %in% c(1,2),"Low","High")
    rownames(temp) <- temp$sample
    clintemp$NEWTHING <- temp[clintemp$SampleID,"hl"]
    colnames(clintemp) <- gsub("NEWTHING",thing,colnames(clintemp))
  clintemp$PFS <- clintemp$`Progression-Free Survival (PFS) in Days` / 7
  colnames(clintemp) <- gsub(palist[1],"One",colnames(clintemp))
  colnames(clintemp) <- gsub(palist[2],"Two",colnames(clintemp))
  clintemp$Call <- ifelse(clintemp$One=="High",ifelse(clintemp$Two=="High","H/H","H/L"),ifelse(clintemp$Two=="High","L/H","L/L"))
  fit6 = survfit(formula=Surv(as.numeric(PFS),as.numeric(factor(PFSReached)))~Call,data=clintemp)
  plot.data=fortify(fit6,surv.connect=TRUE)
  p <- ggplot(plot.data, aes(time, surv,color=strata)) +
    geom_step(size=2) +
    geom_point(data = subset(plot.data, n.censor > 0),size=3,shape=3,aes(color=strata)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values=cols2) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    labs(color="Quartile",x="Time (Months)",y="PFS (%)",title=paste0("PFS with ",co[x,"options"]," in ",co[x,"time"]," (full model p = ",round(unique(df.2$overallp[df.2$y=="PFS" & df.2$time==co[x,"time"] & df.2$options==co[x,"options"]]),4),")\n",palist[1],";",palist[2],"\nMed PFS. First=",round(surv_median(fit6)[[2]][1],1),"(",round(surv_median(fit6)[[3]][1],1),"-",round(surv_median(fit6)[[4]][1],1),"), Second=",round(surv_median(fit6)[[2]][2],1),"(",round(surv_median(fit6)[[3]][2],1),"-",round(surv_median(fit6)[[4]][2],1),")\nThird=",round(surv_median(fit6)[[2]][3],1),"(",round(surv_median(fit6)[[3]][3],1),"-",round(surv_median(fit6)[[4]][3],1),"), Fourth=",round(surv_median(fit6)[[2]][4],1),"(",round(surv_median(fit6)[[3]][4],1),"-",round(surv_median(fit6)[[4]][4],1),")")) +
    theme(legend.position=c(.77,.77),
          legend.key=element_blank(),
          legend.background=element_blank(),
          aspect.ratio=1.0,
          legend.text=element_text(size=12),
          axis.text.y=element_text(size=15),
          axis.text.x=element_text(size=15),
          axis.title=element_text(size=16)) +
    coord_cartesian(xlim=c(0,105),expand=FALSE)
  print(p)
})
dev.off()


































