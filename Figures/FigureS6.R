library(data.table)
library(ggplot2)
library(ggfortify)
library(survival)
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


aux <- fread("data/SampleSourceData.txt",data.table=FALSE)
aux2 <- fread("data/PatientSourceData.txt",data.table=FALSE)


lm_eqn <- function(df){
  colnames(df) <- c("x","y")
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}




pdf("FigureS6_AB.pdf",width=6,height=5)
temptab2 <- do.call(rbind,lapply(c("Baseline","On-Treatment"),function(time) {
  temp <- aux[aux$SampleTimepoint == time & !(aux$Subject %in% aux2[!is.na(aux2$PDbyNTL),"Subject"]) & !is.na(aux$ExpressedNeoantigens),]
  temptab3 <- do.call(rbind,lapply(c("PR","SD","PD"),function(x) {
    temp$BR <- ifelse(temp$BestResponse==x,1,0)
    t <- t.test(temp$ExpressedNeoantigens[temp$BR==0],temp$ExpressedNeoantigens[temp$BR==1])
    data.frame(Response=x,Timepoint=time,NumTrue=length(temp$ExpressedNeoantigens[temp$BR==1]),NumFalse=length(temp$ExpressedNeoantigens[temp$BR==0]),MeanTrue=t$estimate[2],MeanFalse=t$estimate[1],p=t$p.value,stringsAsFactors = FALSE)
  }))
  temptab3$label <- paste0(temptab3$Response,": ",round(temptab3$p,3))
  title <- paste0(paste(temptab3$label,collapse="\n"))
  colnames(temp) <- gsub("Cohort","color",colnames(temp)) 
  temp$BestResponse <- factor(temp$BestResponse,levels=c("PR","SD","PD"))
  p <- ggplot(temp,aes(x=BestResponse,y=ExpressedNeoantigens)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.25,size=5,alpha=0.7,aes(color=color)) +
    labs(x="Best Response",y="# of Expressed\nStrong Binding Neoantigens",title=title) +
    theme_bw() +
    scale_color_manual(values=cols,na.value="gray50") +
    theme(text = element_text(size=14),
          title = element_text(size=12),
          axis.title = element_text(size=14),
          axis.text.x = element_text(angle=45,hjust=1))
  print(p)
  temptab3
}))
dev.off()



pdf("FigureS6CD.pdf")
ggplot(aux[aux$SampleTimepoint=="Baseline" & !is.na(aux$Neoantigens),],aes(x=TMB,y=Neoantigens)) +
  geom_point(size=3) +
  theme_bw() +
  labs(title="All Strong Binders at Baseline",y="# of Strong Binding Neoantigens") +
  geom_smooth(method = "lm", se=FALSE, color="black", linetype="dashed",formula = y ~ x) +
  geom_text(x = 2, y = 100, label = lm_eqn(aux[aux$SampleTimepoint=="Baseline" & !is.na(aux$Neoantigens),c("TMB","Neoantigens")]), parse = TRUE)
ggplot(aux[aux$SampleTimepoint=="Baseline" & !is.na(aux$Neoantigens),],aes(x=TMB,y=ExpressedNeoantigens)) +
  geom_point(size=3) +
  theme_bw() +
  labs(title="All Expressed Strong Binders at Baseline",y="# of Expressed Strong Binding Neoantigens") +
  geom_smooth(method = "lm", se=FALSE, color="black", linetype="dashed",formula = y ~ x) +
  geom_text(x = 2, y = 25, label = lm_eqn(aux[aux$SampleTimepoint=="Baseline" & !is.na(aux$Neoantigens),c("TMB","ExpressedNeoantigens")]), parse = TRUE)
dev.off()


pdf("FigureS6_E.pdf")
temp <- aux2
temp$PFS <- temp$PFS.days / 7
temp <- temp[!is.na(temp$hlaA1101),]
fit6 = survfit(formula=Surv(as.numeric(PFS),as.numeric(factor(PFSCensor)))~hlaA1101,data=temp)
sd <- survdiff(Surv(as.numeric(PFS), PFSCensor==1)~hlaA1101,data=temp)
p.val1 <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
plot.data=fortify(fit6,surv.connect=TRUE)
p <- ggplot(plot.data, aes(time, surv,color=strata)) +
  geom_step(size=2) +
  geom_point(data = subset(plot.data, n.censor > 0),size=3,shape=3,aes(color=strata)) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values=c("blue","gray80")) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(color="HLA-A-11-01",x="Time (Weeks)",y="PFS (%)",title=paste0("PFS with A Allele\np-val=",round(p.val1,3),"\nN: HLA-A-11-01=",length(temp$hlaA1101[temp$hlaA1101=="1"]),", Not=",length(temp$hlaA1101[temp$hlaA1101=="0"]),"\nMed PFS. Mutation=",round(surv_median(fit6)[[2]][1],1),"(",round(surv_median(fit6)[[3]][1],1),"-",round(surv_median(fit6)[[4]][1],1),"), WT=",round(surv_median(fit6)[[2]][2],1),"(",round(surv_median(fit6)[[3]][2],1),"-",round(surv_median(fit6)[[4]][2],1),")")) +
  theme(legend.position=c(.77,.77),
        legend.key=element_blank(),
        legend.background=element_blank(),
        aspect.ratio=1.0,
        legend.text=element_text(size=12),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_text(size=16)) +
  coord_cartesian(xlim=c(0,100),expand=FALSE)
q <- ggsurvplot(fit6, data = temp, risk.table = TRUE,palette=c("blue","gray80"),break.x.by=6,xlim=c(0,100),ylab="PFS")
print(p)
print(q)
data.frame(HLA="A",Allele="HLA-A-11-01",NumAll=length(temp$hlaA1101[temp$hlaA1101=="1"]),NumNot=length(temp$hlaA1101[temp$hlaA1101=="0"]),MedOSAll=surv_median(fit6)[[2]][1],UpOSAll=surv_median(fit6)[[3]][1],DownOSAll=surv_median(fit6)[[4]][1], MedOSNot=surv_median(fit6)[[2]][2],UpOSNot=surv_median(fit6)[[3]][2],DownOSNot=surv_median(fit6)[[4]][2],p=p.val1,stringsAsFactors = FALSE)
dev.off()




pdf("FigureS6_FtoM.pdf",width=8,height=8)
temptab <- do.call(rbind,lapply(colnames(aux)[grepl("diversity",colnames(aux))],function(y) {
  temp <- aux[!(aux$Subject %in% aux2[!is.na(aux2$PDbyNTL),"Subject"]),]
  colnames(temp) <- gsub(y,"Test",colnames(temp))
  temptab2 <- do.call(rbind,lapply(c("Baseline","On-Treatment"),function(time) {
    temp <- temp[temp$SampleTimepoint == time & !is.na(temp$Test),]
    temptab3 <- do.call(rbind,lapply(c("PR","SD","PD"),function(x) {
      temp$BR <- ifelse(temp$BestResponse==x,1,0)
      t <- t.test(temp$Test[temp$BR==0],temp$Test[temp$BR==1])
      data.frame(TCR=y,Response=x,Timepoint=time,NumTrue=length(temp$Test[temp$BR==1]),NumFalse=length(temp$Test[temp$BR==0]),MeanTrue=t$estimate[2],MeanFalse=t$estimate[1],p=t$p.value,stringsAsFactors = FALSE)
    }))
    temptab3$label <- paste0(temptab3$Response,": ",round(temptab3$p,3))
    title <- paste0(temptab3$Timepoint,"\n",paste(temptab3$label,collapse="\n"))
    colnames(temp) <- gsub("Cohort","color",colnames(temp)) 
    temp$BestResponse <- factor(temp$BestResponse,levels=c("PR","SD","PD"))
    p <- ggplot(temp,aes(x=BestResponse,y=Test)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width=0.25,size=5,alpha=0.7,aes(color=color)) +
      labs(x="BestResponse",y=y,title=title) +
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





