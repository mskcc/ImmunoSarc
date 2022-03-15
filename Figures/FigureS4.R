library(data.table)
library(ggplot2)
library(tidyr)


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

oncoprint_cols = c('HOMDEL'='#3953a4',
                   'AMP'='#ed2224',
                   'missense_inframe_driver'='#0e8040',
                   'missense_inframe'='#67BC45',
                   'truncating_driver'='#010101',
                   'truncating'='#010101',
                   'germline_driver'='sienna4',
                   'germline'='sienna3',
                   'FUSION'='#7F17C1',
                   'METH'='yellow3',
                   'Male'='black',
                   'Yes'='black',
                   'No'='#bfbfbf',
                   "Fusion"="orange",
                   "Missing"="white"
)
oncoprint_cols <- c(oncoprint_cols,cols)



aux <- fread("data/PatientSourceData.txt",data.table=FALSE)
supp1 <- fread("TableS1.txt",data.table=FALSE)



aux2 <- aux[!is.na(aux$BestResponse),colnames(aux) %in% c("Subject","Cohort","BestResponse","Fusion","TP53","CDK4/MDM2","NCOR1/MAP2K4","RB1","CDKN2A","ATRX","TERT","PIK3CA","IDH2","KRAS","IDH1","SMARCB1","EWSR1","VEGFA","NF1","MYH9","GLI1")]


aux3 <- gather(aux2,Symbol,call,2:21)

aux3 <- do.call(rbind,lapply(rownames(aux3),function(x) {
  temp <- aux3[x,]
  if(grepl(",",temp$call)) {
    temp2 <- rbind(temp,temp)
    temp2$call[1] <- strsplit(temp$call,",")[[1]][1]
    temp2$call[2] <- strsplit(temp$call,",")[[1]][2]
  } else {
    temp2 <- temp
  }
  temp2
}))


aux3$Subject <- factor(aux3$Subject,levels=aux$Subject[rev(order(aux$`TargetLesion%atBR`))])


aux3$Symbol <- factor(aux3$Symbol,levels=rev(c("Cohort","BestResponse","Fusion","TP53","CDK4/MDM2","NCOR1/MAP2K4","RB1","CDKN2A","ATRX","TERT","PIK3CA","IDH2","KRAS","IDH1","SMARCB1","EWSR1","VEGFA","NF1","MYH9","GLI1")))

aux3$usefacet <- ifelse(aux3$Symbol %in% c("Cohort","BestResponse","Fusion"),"Facet1",ifelse(aux3$Symbol %in% c("TP53","CDK4/MDM2","NCOR1/MAP2K4","RB1","CDKN2A","ATRX","TERT","PIK3CA","IDH2","KRAS","IDH1","SMARCB1","EWSR1"),"Facet2","Facet3"))


pdf("FigureS4_A.pdf",height=5,width=15)
r <- ggplot(aux3, aes(x = Subject, y = Symbol)) +
  geom_tile(fill = alpha('#bfbfbf',1), color = 'white') +
  geom_tile(data = aux3[aux3$call %in% c("HOMDEL","AMP","Missing") | (aux3$Symbol %in% c("Fusion","SampTimepoint","Cohort","BestResponse") & !(aux3$call %in% c("","none"))),], aes(fill = call, y = Symbol), color = 'white', na.rm = T, size = .15) +
  geom_tile(data = aux3[aux3$call %in% c("missense_inframe"),], aes(fill = call, y = Symbol), color = NA, na.rm = T, height = 6.457/20.372, size = .15) +
  geom_tile(data = aux3[aux3$call %in% c("missense_inframe_driver"),], aes(fill = call, y = Symbol), color = NA, na.rm = T, height = 6.457/20.372, size = .15) +
  geom_tile(data = aux3[aux3$call %in% c("truncating_driver","truncating"),], aes(fill = call, y = Symbol), color = NA, na.rm = T, height = 6.457/20.372, size = .15) +
  geom_segment(x=269.5, xend=270.5, y=1.5, yend=2.5, color="white",size=0.25) +
  geom_segment(x=200.5, xend=201.5, y=8.5, yend=9.5, color="white",size=0.25) +
  # coord_fixed(ratio = 40.372/10.214) +
  scale_fill_manual(values = oncoprint_cols,na.value="gray50") +
  scale_color_manual(values = c("#1F78B4","#1F78B4","#1F78B4"),na.value="gray50") +
  theme(axis.line.x = element_blank(), 
        axis.line.y = element_blank(),
        # axis.text.x = element_text(angle=45,hjust=1),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = 'italic'),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_blank(), legend.position = 'bottom', legend.key.size = unit(0.3,"line"),
        plot.margin = unit(c(0,.25,.25,.25), 'lines')) +
  labs(x = '', y= '') +
  facet_grid(usefacet~ ., space = 'free',scales='free')
print(r)
dev.off()


aux$BvOAll <- rowSums(is.na(aux[,c("Baseline_BvO_AllMut","Shared_BvO_AllMut","OnTreatment_BvO_AllMut")]),na.rm=TRUE)
aux$BvODriv <- rowSums(is.na(aux[,c("Baseline_BvO_Driver","Shared_BvO_Driver","OnTreatment_BvO_Driver")]),na.rm=TRUE)
aux$BvPAll <- rowSums(is.na(aux[,c("Baseline_BvP_AllMut","Shared_BvP_AllMut","Progression_BvP_AllMut")]),na.rm=TRUE)


#4 samples had purity estimated by mutation allele frequency. This plot focuses on FACETS purity
supp1$Purity[supp1$Exome %in% c("s_C_MVMX43_M001_d","s_C_2M22RC_M001_d","s_C_H7RF5A_R002_d","s_C_2M22RC_M002_d")] <- NA
supp1b <- supp1[supp1$`Sample Timepoint`=="Baseline",]
supp1o <- supp1[supp1$`Sample Timepoint`=="On-Treatment",]
supp1p <- supp1[supp1$`Sample Timepoint`=="Progression",]
rownames(supp1b) <- supp1b$`Patient ID`
rownames(supp1o) <- supp1o$`Patient ID`
rownames(supp1p) <- supp1p$`Patient ID`

auxBvOAll <- aux[aux$BvOAll!=3,c("PatientID","Cohort","BestResponse","Baseline_BvO_AllMut","Shared_BvO_AllMut","OnTreatment_BvO_AllMut")]
auxBvOAll$FACETSPurityBaseline <- round(supp1b[auxBvOAll$PatientID,"Purity"],2)
auxBvOAll$FACETSPurityOnTreatment <- round(supp1o[auxBvOAll$PatientID,"Purity"],2)
auxBvOAll2 <- gather(auxBvOAll,Category,NumberMut,4:6)
auxBvOAll2$Xaxis <- paste0(auxBvOAll2$Subject,"\n",auxBvOAll2$FACETSPurityBaseline," > ",auxBvOAll2$FACETSPurityOnTreatment)
auxBvOAll2$Cohort <- factor(auxBvOAll2$Cohort,levels=c("UPS/MFH/High Grade MFS","Liposarcoma","Chondrosarcoma","Leiomyosarcoma","Osteosarcoma","Vascular","Small Blue Round Cell","ASPS","Other"))
auxBvOAll2$Category <- factor(auxBvOAll2$Category,levels=c("Baseline_BvO_AllMut","Shared_BvO_AllMut","OnTreatment_BvO_AllMut"))
auxBvOAll2$BestResponse <- factor(auxBvOAll2$BestResponse,levels=c("PR","SD","PD"))
auxBvOAll2$Xaxis <- factor(auxBvOAll2$Xaxis,levels=unique(auxBvOAll2$Xaxis[order(auxBvOAll2$BestResponse)]))


auxBvODriv <- aux[aux$BvODriv!=3,c("PatientID","Cohort","BestResponse","Baseline_BvO_Driver","Shared_BvO_Driver","OnTreatment_BvO_Driver")]
auxBvODriv$FACETSPurityBaseline <- round(supp1b[auxBvODriv$PatientID,"Purity"],2)
auxBvODriv$FACETSPurityOnTreatment <- round(supp1o[auxBvODriv$PatientID,"Purity"],2)
auxBvODriv2 <- gather(auxBvODriv,Category,NumberMut,4:6)
auxBvODriv2$Xaxis <- paste0(auxBvODriv2$Subject,"\n",auxBvODriv2$FACETSPurityBaseline," > ",auxBvODriv2$FACETSPurityOnTreatment)
auxBvODriv2$Cohort <- factor(auxBvODriv2$Cohort,levels=c("UPS/MFH/High Grade MFS","Liposarcoma","Chondrosarcoma","Leiomyosarcoma","Osteosarcoma","Vascular","Small Blue Round Cell","ASPS","Other"))
auxBvODriv2$Category <- factor(auxBvODriv2$Category,levels=c("Baseline_BvO_Driver","Shared_BvO_Driver","OnTreatment_BvO_Driver"))
auxBvODriv2$BestResponse <- factor(auxBvODriv2$BestResponse,levels=c("PR","SD","PD"))
auxBvODriv2$Xaxis <- factor(auxBvODriv2$Xaxis,levels=unique(auxBvODriv2$Xaxis[order(auxBvODriv2$BestResponse)]))


pdf("FiguresS4_BC.pdf",height=7,width=25)
ggplot(auxBvOAll2,aes(x=Xaxis,y=NumberMut,fill=Category)) +
  geom_segment(aes(y = -0.10,yend=-0.05,x=Xaxis,xend=Xaxis, color = BestResponse),size=10) + 
  geom_bar(stat="identity",position="fill",alpha=0.8) +
  theme_bw() +
  geom_text(color="black",aes(label=NumberMut),position=position_fill(),vjust=1,size=6) +
  theme(axis.text.x = element_text(hjust=1,angle=45),
        panel.grid = element_blank(),
        axis.text.y=element_text(size=16),
        axis.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=14)) +
  geom_abline(slope=0,intercept=0,color="#5F5F5F",size=0.5) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=c("red","chocolate1","gold")) +
  labs(y="Percent of Mutations in Patient",color="Best Response",fill="Sharing of Mutations") +
  facet_grid(. ~ Cohort,scales="free",space="free",labeller=label_wrap_gen())
ggplot(auxBvODriv2,aes(x=Xaxis,y=NumberMut,fill=Category)) +
  geom_segment(aes(y = -0.10,yend=-0.05,x=Xaxis,xend=Xaxis, color = BestResponse),size=10) + 
  geom_bar(stat="identity",position="fill",alpha=0.8) +
  theme_bw() +
  geom_text(color="black",aes(label=NumberMut),position=position_fill(),vjust=1,size=6) +
  theme(axis.text.x = element_text(hjust=1,angle=45),
        panel.grid = element_blank(),
        axis.text.y=element_text(size=16),
        axis.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=14)) +
  geom_abline(slope=0,intercept=0,color="#5F5F5F",size=0.5) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=c("red","chocolate1","gold")) +
  labs(y="Percent of Mutations in Patient",color="Best Response",fill="Sharing of Mutations") +
  facet_grid(. ~ Cohort,scales="free",space="free",labeller=label_wrap_gen())
dev.off()


auxBvPAll <- aux[aux$BvPAll!=3,c("PatientID","Cohort","BestResponse","Baseline_BvP_AllMut","Shared_BvP_AllMut","Progression_BvP_AllMut")]
auxBvPAll$FACETSPurityBaseline <- round(supp1b[auxBvPAll$PatientID,"Purity"],2)
auxBvPAll$FACETSPurityProgression <- round(supp1p[auxBvPAll$PatientID,"Purity"],2)
auxBvPAll2 <- gather(auxBvPAll,Category,NumberMut,4:6)
auxBvPAll2$Xaxis <- paste0(auxBvPAll2$Subject,"\n",auxBvPAll2$FACETSPurityBaseline," > ",auxBvPAll2$FACETSPurityProgression)
auxBvPAll2$Cohort <- factor(auxBvPAll2$Cohort,levels=c("UPS/MFH/High Grade MFS","Liposarcoma","Chondrosarcoma","Leiomyosarcoma","Osteosarcoma","Vascular","Small Blue Round Cell","ASPS","Other"))
auxBvPAll2$Category <- factor(auxBvPAll2$Category,levels=c("Progression_BvP_AllMut","Shared_BvP_AllMut","Baseline_BvP_AllMut"))
auxBvPAll2$BestResponse <- factor(auxBvPAll2$BestResponse,levels=c("PR","SD","PD"))
auxBvPAll2$Xaxis <- factor(auxBvPAll2$Xaxis,levels=unique(auxBvPAll2$Xaxis[order(auxBvPAll2$BestResponse)]))


pdf("FigureS4_D.pdf",height=5,width=10)
ggplot(auxBvPAll2,aes(x=Xaxis,y=NumberMut,fill=Category)) +
  geom_segment(aes(y = -0.10,yend=-0.05,x=Xaxis,xend=Xaxis, color = BestResponse),size=10) + 
  geom_bar(stat="identity",position="fill",alpha=0.8) +
  theme_bw() +
  geom_text(color="black",aes(label=NumberMut),position=position_fill(),vjust=1,size=6) +
  theme(axis.text.x = element_text(hjust=1,angle=45),
        panel.grid = element_blank(),
        axis.text.y=element_text(size=16),
        axis.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=14)) +
  geom_abline(slope=0,intercept=0,color="#5F5F5F",size=0.5) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=c("dodgerblue3","darkorchid3","red")) +
  labs(y="Percent of Mutations in Patient",color="Best Response",fill="Sharing of Mutations")
dev.off()



