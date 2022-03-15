##Script for Kaplan-Meier Plots
##Relevant to Figure 2, and supplemental figures 1, 5 and 6
##Author - Allison L Richards


library(data.table)
library(dplyr)
library(tidyr)
library(survminer)
library(lubridate)
library(survival)
library(ggfortify)
library(RColorBrewer)



clin <- fread("XXX.txt",data.table=FALSE)


#This script will work for all above figures but as an example, we'll compare quartiles of TMB which is derived directly from Tempo

clin <- clin[!is.na(clin$`Sample Timepoint`) & clin$`Sample Timepoint`=="Baseline",] #TMB was only characterized at baseline
clin$TMBquart <- ifelse(clin$TMB < quantile(clin$TMB,0.25,na.rm=TRUE),1,
                         ifelse(clin$TMB < quantile(clin$TMB,0.5,na.rm=TRUE),2,
                                ifelse(clin$TMB < quantile(clin$TMB,0.75,na.rm=TRUE),3,4)))


pdf("PFS_TMBquart.pdf")
aux <- do.call(rbind,lapply(c("TMBquart"),function(x) {
  temp <- clin
  temp$PFS <- temp$`Progression-Free Survival (PFS) in Days` / 7
  colnames(temp) <- gsub(x,"Test",colnames(temp))
  fit6 = survfit(formula=Surv(as.numeric(PFS),as.numeric(factor(PFSReached)))~Test,data=temp)
  ttab <- do.call(rbind,lapply(c("1,2","1,3","1,4","2,3","2,4","3,4"),function(pair) {
    first <- strsplit(pair,",")[[1]][1]
    second <- strsplit(pair,",")[[1]][2]
    temp2 <- temp[temp$Test %in% c(first,second),]
    sd <- survdiff(Surv(as.numeric(PFS), PFSReached==1)~Test,data=temp2)
    p.val1 <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
    data.frame(first=first,second=second,p=p.val1,stringsAsFactors = FALSE)
  }))
  ttab$label <- paste0(ttab$first," vs ",ttab$second,": ",round(ttab$p,3))
  linklabel <- paste(ttab$label,collapse="\n")
  plot.data=fortify(fit6,surv.connect=TRUE)
  p <- ggplot(plot.data, aes(time, surv,color=strata)) +
    geom_step(size=2) +
    geom_point(data = subset(plot.data, n.censor > 0),size=3,shape=3,aes(color=strata)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values=brewer.pal(4,"Reds")) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    labs(color="Quartile",x="Time (Months)",y="PFS (%)",title=paste0("PFS with ",x,"\nMed PFS. First=",round(surv_median(fit6)[[2]][1],1),"(",round(surv_median(fit6)[[3]][1],1),"-",round(surv_median(fit6)[[4]][1],1),"), Second=",round(surv_median(fit6)[[2]][2],1),"(",round(surv_median(fit6)[[3]][2],1),"-",round(surv_median(fit6)[[4]][2],1),")\nThird=",round(surv_median(fit6)[[2]][3],1),"(",round(surv_median(fit6)[[3]][3],1),"-",round(surv_median(fit6)[[4]][3],1),"), Fourth=",round(surv_median(fit6)[[2]][4],1),"(",round(surv_median(fit6)[[3]][4],1),"-",round(surv_median(fit6)[[4]][4],1),")\n",linklabel)) +
    theme(legend.position=c(.77,.77),
          legend.key=element_blank(),
          legend.background=element_blank(),
          aspect.ratio=1.0,
          legend.text=element_text(size=12),
          axis.text.y=element_text(size=15),
          axis.text.x=element_text(size=15),
          axis.title=element_text(size=16)) +
    coord_cartesian(xlim=c(0,105),expand=FALSE)
  q <- ggsurvplot(fit6, data = temp, risk.table = TRUE,palette=brewer.pal(4,"Reds"),break.x.by=6,xlim=c(0,105),ylab="PFS")
  print(p)
  print(q)
  data.frame(Test=x,NumFirst=length(temp$`CMO Patient ID`[temp$Test=="1"]),NumSecond=length(temp$`CMO Patient ID`[temp$Test=="2"]),NumThird=length(temp$`CMO Patient ID`[temp$Test=="3"]),NumFourth=length(temp$`CMO Patient ID`[temp$Test=="4"]),MedFirst=surv_median(fit6)[[2]][1],MedFirstDown=surv_median(fit6)[[3]][1],MedFirstUp=surv_median(fit6)[[4]][1],MedSecond=surv_median(fit6)[[2]][2],MedSecondDown=surv_median(fit6)[[3]][2],MedSecondUp=surv_median(fit6)[[4]][2],MedThird=surv_median(fit6)[[2]][3],MedThirdDown=surv_median(fit6)[[3]][3],MedThirdUp=surv_median(fit6)[[4]][3],MedFourth=surv_median(fit6)[[2]][4],MedFourthDown=surv_median(fit6)[[3]][4],MedFourthUp=surv_median(fit6)[[4]][4],p12=ttab$p[ttab$first==1&ttab$second==2],p13=ttab$p[ttab$first==1&ttab$second==3],p14=ttab$p[ttab$first==1&ttab$second==4],p23=ttab$p[ttab$first==2&ttab$second==3],p24=ttab$p[ttab$first==2&ttab$second==4],p34=ttab$p[ttab$first==3&ttab$second==4],stringsAsFactors = FALSE)
}))
dev.off()











