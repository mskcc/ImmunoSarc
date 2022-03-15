##Script to Organize Data and Get Immune Population Content
##Relevant to all figures with a slightly more annotated version of the heatmap being used in Figure 2
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
library(immunedeconv)
library(ComplexHeatmap)


clin <- fread("XXX.txt",data.table=FALSE)

kall #List of RNAseq abundance files from dbGap
kall2 <- sapply(strsplit(kall,"/"),function(x) {x[grep("abundance",x)-1]})


des <- do.call(rbind,lapply(rownames(clin),function(x) {
  tab <- clin[x,]
  if(tab$SampleID %in% kall2) {
    tab$path <- gsub("/abundance.tsv","",kall[grepl(tab$SampleID,kall)])
    colnames(tab) <- gsub("SampleID","sample",colnames(tab))
    tab
  }
}))

#Use "des" to run sleuth via Sleuth_Prep_Run.R
#bsub -n 8 -R rusage[mem=40] -W 1:59 'Rscript Sleuth_Prep_Run.R -d des.txt -t Gene -o Sleuth_Prep_Result.Rdata'

#Read in sleuth_prep results
full <- sleuth_load("Sleuth_Prep_Result.Rdata")

mat <- as.data.frame(sleuth_to_matrix(full,"obs_norm","tpm"))


##Use immunodeconv to characterize immune cell populations via MCPcounter
res_mcp = deconvolute(mat, "mcp_counter", tumor = TRUE, arrays = FALSE, scale_mrna = FALSE)
# >>> Running mcp_counter
pop1 <- as.data.frame(res_mcp)

rownames(pop1) <- pop1$cell_type
pop1 <- pop1[,-1]
pop1 <- as.data.frame(t(pop1))
colnames(pop1) <- gsub(" ",".",colnames(pop1))
colnames(pop1) <- gsub("[+]","",colnames(pop1))
colnames(pop1) <- gsub("/",".",colnames(pop1))
pop1$SampleID <- rownames(pop1)

clin2 <- merge(clin,pop1,by="SampleID",all=TRUE)


##To identify "Groups" of samples at Baseline and On-treatment time points, we performed hierarchical clustering on these results and broke samples into 2 or 3 groups by manual review of the resulting heatmaps
pdf("IdentifyGroupsOfSamples.pdf",width=12,height=8)
lapply(c("Baseline","C2"),function(type) { 
    temp.1 <- clin2[clin2$`Sample Timepoint`==type & clin2$SampleID %in% kall2,]
    temp <- temp.1
    roword <- c("T.cell.CD8","T.cell","B.cell","Macrophage.Monocyte","Neutrophil","Myeloid.dendritic.cell","NK.cell","Endothelial.cell","Cancer.associated.fibroblast")
    temp$BestResponse <- factor(temp$BestResponse,levels=c("PR","SD","PDbyNTL","PD"))
    temp <- temp[,roword]
    temp <- t(scale(temp,center=TRUE,scale=TRUE))
    colnames(temp) <- temp.1$`CMO Patient ID`
    temp <- temp[roword,]
    rownames(temp.1) <- temp.1$`CMO Patient ID`
    column_dend = as.dendrogram(hclust(dist(t(temp), method = "manhattan"), method = "ward.D"))
    column_dendlog.1 = hclust(dist(t(temp), method = "manhattan"), method = "ward.D")
    column_dend = dendextend::rotate(column_dend, temp.1$`CMO Patient ID`[order(temp.1$BestResponse)])
    ha = HeatmapAnnotation(show_annotation_name = T,
                           show_legend=T,
                           annotation_name_gp = gpar(fontsize = 12),
                           BestResponse=temp.1[,"BestResponse"],
                           Cohort=temp.1[,"Cohort"],
                           Samptype=temp.1[,"Sample Timepoint"],
                           col=list(BestResponse=cols,
                                    Cohort=cols,
                                    Samptype=cols),
                           gp = gpar(fontsize = 18, col="black"),
                           na_col="white")
    h4=Heatmap(as.matrix(temp),
               heatmap_legend_param = list(title = "Zscore"),
               column_dend_height = unit(3, "cm"),
               column_dend_reorder =T,
               cluster_columns=column_dend,
               cluster_rows=FALSE,
               row_dend_reorder = FALSE,
               gap = unit(2, "mm"),
               top_annotation = ha
    )
    print(h4)
})
dev.off()


##After manual inspection, we decided that 3 groups best described samples at baseline and 2 groups best described samples on-treatment (C2)
temp.1 <- clin2[clin2$`Sample Timepoint`=="Baseline" & clin2$SampleID %in% kall2,]
temp <- temp.1
roword <- c("T.cell.CD8","T.cell","B.cell","Macrophage.Monocyte","Neutrophil","Myeloid.dendritic.cell","NK.cell","Endothelial.cell","Cancer.associated.fibroblast")
temp$BestResponse <- factor(temp$BestResponse,levels=c("PR","SD","PDbyNTL","PD"))
temp <- temp[,roword]
temp <- t(scale(temp,center=TRUE,scale=TRUE))
colnames(temp) <- temp.1$`CMO Patient ID`
temp <- temp[roword,]
column_dend = as.dendrogram(hclust(dist(t(temp), method = "manhattan"), method = "ward.D"))
column_dendlog.1 = hclust(dist(t(temp), method = "manhattan"), method = "ward.D")
baselinegroups <- cutree(column_dendlog.1,3) %>% as.data.frame()


temp.1 <- clin2[clin2$`Sample Timepoint`=="C2" & clin2$SampleID %in% kall2,]
temp <- temp.1
roword <- c("T.cell.CD8","T.cell","B.cell","Macrophage.Monocyte","Neutrophil","Myeloid.dendritic.cell","NK.cell","Endothelial.cell","Cancer.associated.fibroblast")
temp$BestResponse <- factor(temp$BestResponse,levels=c("PR","SD","PDbyNTL","PD"))
temp <- temp[,roword]
temp <- t(scale(temp,center=TRUE,scale=TRUE))
colnames(temp) <- temp.1$`CMO Patient ID`
temp <- temp[roword,]
column_dend = as.dendrogram(hclust(dist(t(temp), method = "manhattan"), method = "ward.D"))
column_dendlog.1 = hclust(dist(t(temp), method = "manhattan"), method = "ward.D")
c2groups <- cutree(column_dendlog.1,2) %>% as.data.frame()









