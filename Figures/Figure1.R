library(data.table)
library(ggplot2)
library(ggpubr)




cols <- c("Leiomyosarcoma"="#8DD3C7",
          "UPS/MFH/High Grade MFS"="#FFED6F",
          "Osteosarcoma"="#BEBADA",
          "Chondrosarcoma"="#FB8072",
          "Liposarcoma"="#80B1D3",
          "Other"="gray50",
          "Small Blue Round Cell"='#e6194b',
          "Vascular"='#3cb44b',
          "ASPS"='#9a6324'
)


aux <- fread("data/PatientSourceData.txt",data.table = FALSE)

aux <- aux[rev(order(aux$`TargetLesion%atBR`)),]
aux$Subject <- factor(aux$Subject,levels=aux$Subject)
aux$PFS <- aux$PFS.days / 30


pdf("Figure1BC.pdf",width=15,height=6)
p <- ggplot(aux,aes(x=Subject,y=`TargetLesion%atBR`,fill=Cohort)) +
  geom_bar(stat="identity") +
  geom_point(shape=8,aes(x=Subject,y=NonEvaluable)) +
  geom_point(shape=19,aes(x=Subject,y=PDbyNTL)) +
  geom_text(aes(x=Subject,y=90,label=Label1),color="black",angle=90,size=3) +
  geom_hline(yintercept=-30,linetype="dashed") +
  scale_fill_manual(values=cols) +
  theme_bw() +
  coord_cartesian(y=c(-100,100)) +
  labs(y="% Best Change") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=8),
        axis.ticks.x = element_blank(),
        legend.position = "none")
q <- ggplot(aux,aes(x=Subject,y=PFS)) +
  geom_bar(stat="identity",fill="black") +
  geom_text(aes(x=Subject,y=Label2Position,label=Label2,color=as.factor(Label2Color)),angle=90,size=3) +
  theme_bw() +
  scale_color_manual(values=c("white","black")) +
  coord_cartesian(y=c(0,12)) +
  geom_hline(yintercept=5.6,linetype="dashed") +
  labs(y="Progression-Free Survival (PFS) in Months") +
  guides(color="none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=8))
pq <- ggarrange(p,q,nrow = 2,align="v",heights=c(1,1))
print(pq)
dev.off()

