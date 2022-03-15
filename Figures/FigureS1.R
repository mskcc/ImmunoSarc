library(data.table)
library(survival)
library(ggfortify)
library(survminer)
library(VennDiagram)

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

col2 <- cols
names(col2) <- paste0("Cohort=",names(col2))

aux <- fread("data/PatientSourceData.txt",data.table = FALSE)

aux$PFS <- aux$PFS.days / 30
aux$Alive <- aux$Alive.days / 30


pdf("FigureS1_AB.pdf",width=11,height=10)
fit6 = survfit(formula=Surv(as.numeric(PFS),as.numeric(factor(PFSCensor)))~Cohort,data=aux)
plot.data=fortify(fit6,surv.connect=TRUE)
q <- ggsurvplot(fit6, data = aux, risk.table = TRUE,palette=col2[order(names(col2))],break.x.by=5,xlim=c(0,25),ylab="%Progression-Free",xlab="Months Since Treatment Start",fun="pct",size=2,font.x = c(20), font.y = c(20),font.tickslab = c(20),risk.table.height=0.35,risk.table.font=c(7))
print(q)
fit6 = survfit(formula=Surv(as.numeric(Alive),as.numeric(factor(AliveCensor)))~Cohort,data=aux)
plot.data=fortify(fit6,surv.connect=TRUE)
q <- ggsurvplot(fit6, data = aux, risk.table = TRUE,palette=col2[order(names(col2))],break.x.by=5,xlim=c(0,30),ylab="%Alive",xlab="Months Since Treatment Start",fun="pct",size=2,font.x = c(20), font.y = c(20),font.tickslab = c(20),risk.table.height=0.35,risk.table.font=c(7))
print(q)
dev.off()



pdf("FigureS1_C.pdf")
q <- draw.quad.venn(area1=length(aux$Subject[aux$RNA==1]),
                    area2=length(aux$Subject[aux$TCR==1]),
                    area3=length(aux$Subject[aux$WES==1]),
                    area4=length(aux$Subject[aux$IHC==1]),
                    n12=length(aux$Subject[aux$RNA==1 & aux$TCR==1]),
                    n13=length(aux$Subject[aux$RNA==1 & aux$WES==1]),
                    n14=length(aux$Subject[aux$RNA==1 & aux$IHC==1]),
                    n23=length(aux$Subject[aux$TCR==1 & aux$WES==1]),
                    n24=length(aux$Subject[aux$TCR==1 & aux$IHC==1]),
                    n34=length(aux$Subject[aux$WES==1 & aux$IHC==1]),
                    n123=length(aux$Subject[aux$RNA==1 & aux$TCR==1 & aux$WES==1]),
                    n124=length(aux$Subject[aux$RNA==1 & aux$TCR==1 & aux$IHC==1]),
                    n134=length(aux$Subject[aux$RNA==1 & aux$WES==1 & aux$IHC==1]),
                    n234=length(aux$Subject[aux$TCR==1 & aux$WES==1 & aux$IHC==1]),
                    n1234=length(aux$Subject[aux$RNA==1 & aux$TCR==1 & aux$WES==1 & aux$IHC==1]),
                    fill=c("red","blue","yellow","gray50"),
                    category=c("RNA","TCR","WES","IHC"),
                    cex=rep(3,15),
                    cat.cex=rep(3,4))
print(q)
dev.off()

