library(data.table)
library(ggplot2)
library(survival)
library(ggfortify)
library(RColorBrewer)
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



aux <- fread("data/PatientSourceData.txt",data.table=FALSE)



pdf("FigureS5_AD.pdf",width=8,height=8)
temptab <- do.call(rbind,lapply(c("BaselineTMB","BaselineFGA"),function(y) {
  temp <- aux[is.na(aux$PDbyNTL),]
  colnames(temp) <- gsub(y,"Test",colnames(temp))
  temp <- temp[!is.na(temp$Test),]
  temptab3 <- do.call(rbind,lapply(c("PR","SD","PD"),function(x) {
    temp$BR <- ifelse(temp$BestResponse==x,1,0)
    t <- t.test(temp$Test[temp$BR==0],temp$Test[temp$BR==1])
    data.frame(Response=x,Test=y,NumTrue=length(temp$Test[temp$BR==1]),NumFalse=length(temp$Test[temp$BR==0]),MeanTrue=t$estimate[2],MeanFalse=t$estimate[1],p=t$p.value,stringsAsFactors = FALSE)
  }))
  temptab3$label <- paste0(temptab3$Response,": ",round(temptab3$p,3))
  title <- paste0(paste(temptab3$label,collapse="\n"))
  colnames(temp) <- gsub("Cohort","color",colnames(temp)) #whatever the test is will already be changed to Group
  temp$BestResponse <- factor(temp$BestResponse,levels=c("PR","SD","PD"))
  p <- ggplot(temp,aes(x=BestResponse,y=Test)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.25,size=5,alpha=0.7,aes(color=color)) +
    labs(x="Cohort",y=y,title=title) +
    theme_bw() +
    scale_color_manual(values=cols,na.value="gray50") +
    theme(text = element_text(size=18),
          title = element_text(size=12),
          axis.title = element_text(size=18),
          axis.text.x = element_text(angle=45,hjust=1))
  print(p)
  temptab3
}))
dev.off()



pdf("FigureS5_BE.pdf",width=8,height=8)
temptab2 <- do.call(rbind,lapply(c("BaselineTMB","BaselineFGA"),function(y) {
  temp <- aux
  colnames(temp) <- gsub(y,"Test",colnames(temp))
  temp <- temp[!is.na(temp$Test),]
  tests <- unique(temp$Cohort[!(temp$Cohort %in% c("Other"))])
  counttest <- plyr::count(temp$Cohort)
  tests <- tests[tests %in% counttest$x[counttest$freq >= 3]]
  temptab3 <- do.call(rbind,lapply(tests,function(x) {
    temp$Co <- ifelse(temp$Cohort==x,1,0)
    t <- t.test(temp$Test[temp$Co==0],temp$Test[temp$Co==1])
    data.frame(Cohort=x,Test=y,NumTrue=length(temp$Test[temp$Co==1]),NumFalse=length(temp$Test[temp$Co==0]),MeanTrue=t$estimate[2],MeanFalse=t$estimate[1],p=t$p.value,stringsAsFactors = FALSE)
  }))
  temptab3$label <- paste0(temptab3$Cohort,": ",round(temptab3$p,3))
  title <- paste0(paste(temptab3$label,collapse="\n"))
  colnames(temp) <- gsub("BestResponse","color",colnames(temp)) #whatever the test is will already be changed to Group
  temp$Cohort <- factor(temp$Cohort,levels=c("Liposarcoma","Chondrosarcoma","Osteosarcoma","UPS/MFH/High Grade MFS","Leiomyosarcoma","Vascular","Small Blue Round Cell","ASPS","Other"))
  p <- ggplot(temp,aes(x=Cohort,y=Test)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.25,size=5,alpha=0.7,aes(color=color)) +
    labs(x="Cohort",y=y,title=title) +
    theme_bw() +
    scale_color_manual(values=cols,na.value="gray50") +
    theme(text = element_text(size=18),
          title = element_text(size=12),
          axis.title = element_text(size=18),
          axis.text.x = element_text(angle=45,hjust=1))
  print(p)
  temptab3
}))
dev.off()



pdf("FigureS5_CF.pdf")
aux <- do.call(rbind,lapply(c("BaselineTMBQuart","BaselineFGAQuart"),function(x) {
  temp <- aux
  temp$PFS <- temp$PFS.days / 7
  colnames(temp) <- gsub(x,"Test",colnames(temp))
  temp <- temp[!is.na(temp$Test),]
  temp$PFSReached <- ifelse(grepl("[+]",temp$Label2),0,1)
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
    labs(color="Quartile",x="Time (Weeks)",y="PFS (%)",title=paste0("PFS with ",x,"\nMed PFS. First=",round(surv_median(fit6)[[2]][1],1),"(",round(surv_median(fit6)[[3]][1],1),"-",round(surv_median(fit6)[[4]][1],1),"), Second=",round(surv_median(fit6)[[2]][2],1),"(",round(surv_median(fit6)[[3]][2],1),"-",round(surv_median(fit6)[[4]][2],1),")\nThird=",round(surv_median(fit6)[[2]][3],1),"(",round(surv_median(fit6)[[3]][3],1),"-",round(surv_median(fit6)[[4]][3],1),"), Fourth=",round(surv_median(fit6)[[2]][4],1),"(",round(surv_median(fit6)[[3]][4],1),"-",round(surv_median(fit6)[[4]][4],1),")\n",linklabel)) +
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

