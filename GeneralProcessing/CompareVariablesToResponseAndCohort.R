##Script describing comparisons of various metrics across responses and cohorts
#Relevant to Figure 2, and supplemental figures 2, 3, 5, and 6 and supplemental table 1
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
library(ggplot2)


clin <- fread("XXX.txt",data.table=FALSE)


##Throughout the manuscript, we compare many variables to sarcoma cohort and best response
#Diversity of TCR receptors - can be retrieved from following TCRProcessing.R (S6F-S6M)
#IHC markers and change in IHC markers - can be derived from data deposited to dbGap (2A-C, S2A-S2G)
#CD8 Tcells - can be derived using SetUpData.R (S3B-S3C)
#Hedgehog Signaling - can be derived using GSEACalculations.R (S3D-S3E)
#Gene expression - can be derived from IdentifyingDEG.R (S3A, S3F)
#TMB and FGA - results from Tempo (S5A, S5B, S5D, S5E)
#Expressed Neoantigens - can be found taking results from Tempo and using snp-pileup to identify neoantigens with at least 1 supporting read in the corresponding RNAseq bam (S6A, S6B)


##This test case will be shown with data for IHC markers but can be extrapolated to the above variables
test <- colnames(clin)[grepl("%",colnames(clin))]



pdf("BoxPlots_AcrossCohort.pdf",width=8,height=8)
cohtab <- do.call(rbind,lapply(c("Cohort"),function(x) {
  temptab <- do.call(rbind,lapply(test,function(y) {
    temp <- clin
    colnames(temp) <- gsub(x,"Group",gsub(y,"Test",colnames(temp)))
    temptab2 <- do.call(rbind,lapply(c("Baseline","C2"),function(time) {
      temp <- temp[temp$`Sample Timepoint` == time & !is.na(temp$`Sample Timepoint`),]
      temp <- temp[!is.na(temp$Group) & !is.na(temp$Test) & temp$Group != "PDbyNTL",]
      tests <- unique(temp$Group[!(temp$Group %in% c("Other","PDbyNTL"))])
      temp2 <- plyr::count(temp$Group)
      tests <- tests[tests %in% temp2$x[temp2$freq >= 3]]
      temptab3 <- do.call(rbind,lapply(tests,function(piece) {
        temp$Piece <- ifelse(temp$Group==piece,1,0)
        t <- t.test(temp$Test[temp$Piece==0],temp$Test[temp$Piece==1])
        data.frame(Piece=piece,Test=y,Time=time,NumTrue=length(temp$Test[temp$Piece==1]),NumFalse=length(temp$Test[temp$Piece==0]),MeanTrue=t$estimate[2],MeanFalse=t$estimate[1],p=t$p.value,stringsAsFactors = FALSE)
      }))
      temptab3$label <- paste0(temptab3$Piece,": ",round(temptab3$p,3))
      title <- paste0("Time: ",time,"\n",paste(temptab3$label,collapse="\n"))
      colnames(temp) <- gsub("BestResponse|Cohort","color",colnames(temp)) #whatever the test is will already be changed to Group
      temp$Group <- factor(temp$Group,levels=c("PR","SD","PDbyNTL","PD","Liposarcoma (LPS)","Chondrosarcoma","Osteosarcoma","UPS/MFH/High Grade MFS","Leiomyosarcoma (LMS)","Vascular","Small Blue Round Cell (SBRC)","ASPS","Other"))
      p <- ggplot(temp,aes(x=as.factor(Group),y=Test)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width=0.25,size=5,alpha=0.7,aes(color=color)) +
        labs(x=x,y=y,title=title) +
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
  temptab
}))
dev.off()


pdf("BoxPlots_AcrossResponses.pdf",width=8,height=8)
besttab <- do.call(rbind,lapply(c("BestResponse"),function(x) {
  temptab <- do.call(rbind,lapply(test,function(y) {
    temp <- clin
    colnames(temp) <- gsub(x,"Group",gsub(y,"Test",colnames(temp)))
    temptab2 <- do.call(rbind,lapply(c("Baseline","C2"),function(time) {
      temp <- temp[temp$`Sample Timepoint` == time & !is.na(temp$`Sample Timepoint`),]
      temp <- temp[!is.na(temp$Group) & !is.na(temp$Test) & temp$Group != "PDbyNTL",]
      tests <- unique(temp$Group[!(temp$Group %in% c("Other","PDbyNTL"))])
      temp2 <- plyr::count(temp$Group)
      tests <- tests[tests %in% temp2$x[temp2$freq >= 3]]
      temptab3 <- do.call(rbind,lapply(tests,function(piece) {
        temp$Piece <- ifelse(temp$Group==piece,1,0)
        t <- t.test(temp$Test[temp$Piece==0],temp$Test[temp$Piece==1])
        data.frame(Piece=piece,Test=y,Time=time,NumTrue=length(temp$Test[temp$Piece==1]),NumFalse=length(temp$Test[temp$Piece==0]),MeanTrue=t$estimate[2],MeanFalse=t$estimate[1],p=t$p.value,stringsAsFactors = FALSE)
      }))
      temptab3$label <- paste0(temptab3$Piece,": ",round(temptab3$p,3))
      title <- paste0("Time: ",time,"\n",paste(temptab3$label,collapse="\n"))
      colnames(temp) <- gsub("BestResponse|Cohort","color",colnames(temp)) #whatever the test is will already be changed to Group
      temp$Group <- factor(temp$Group,levels=c("PR","SD","PDbyNTL","PD","Liposarcoma (LPS)","Chondrosarcoma","Osteosarcoma","UPS/MFH/High Grade MFS","Leiomyosarcoma (LMS)","Vascular","Small Blue Round Cell (SBRC)","ASPS","Other"))
      p <- ggplot(temp,aes(x=as.factor(Group),y=Test)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width=0.25,size=5,alpha=0.7,aes(color=color)) +
        labs(x=x,y=y,title=title) +
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
  temptab
}))
dev.off()





#Figures with both Baseline and On-treatment samples together
plot <- function(x,test) {
  temp <- clin
  colnames(temp) <- gsub(x,"Feature",colnames(temp))
  colnames(temp) <- gsub(test,"Test",colnames(temp))
  temp <- temp[!is.na(temp$`Sample Timepoint`) & temp$`Sample Timepoint`!="Progression" & !is.na(temp$Feature) & !is.na(temp$Test),]
  testname <- "Best Response"
  p <- ggplot(temp,aes(x=Test,y=Feature)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=3,aes(color=Cohort),width = 0.2) +
    scale_color_manual(values=cols) +
    theme_bw() +
    labs(x=testname,y=x) +
    theme(text = element_text(size=20)) +
    facet_grid(. ~ `Sample Timepoint`,scales = "free",space="free")
  print(p)
}

#For these plots, we combine PD and PDbyNTL as they are all progression samples
clin$BR <- clin$BestResponse
clin$BR <- ifelse(clin$BR=="PDbyNTL","PD",clin$BR)
clin$BR <- factor(clin$BR,levels=c("PR","SD","PD"))

pdf("IHCPlots_CompareResponse.pdf",width=10,height=5)
lapply(colnames(clin)[grepl("% Positive",colnames(clin))],plot,test="BR")
dev.off()




##Figures to compare the change in IHC markers
plot2 <- function(x,test) {
  temp <- clin
  colnames(temp) <- gsub(x,"Feature",colnames(temp))
  colnames(temp) <- gsub(test,"Test",colnames(temp))
  temp <- temp[!is.na(temp$`Sample Timepoint`) & temp$`Sample Timepoint`!="Progression" & !is.na(temp$Feature) & !is.na(temp$Test),]
  testname <- "Best Response"
  temp2 <- spread(temp[,c("CMO Patient ID","Test","Cohort","Sample Timepoint","Feature")],`Sample Timepoint`,Feature)
  temp2$num <- rowSums(!is.na(temp2[,c("Baseline","C2")]))
  temp2 <- temp2[temp2$num == 2,]
  temp2$Change <- temp2$C2 - temp2$Baseline
  p <- ggplot(temp2,aes(x=Test,y=Change)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=3,aes(color=Cohort),width = 0.2) +
    scale_color_manual(values=cols) +
    theme_bw() +
    labs(x=testname,y=paste0("Change in\n",x)) +
    theme(text = element_text(size=20))
  print(p)
}


pdf("IHCPlots_CompareChangeInResponse.pdf",width=8,height=5)
lapply(colnames(clin)[grepl("% Positive",colnames(clin))],plot2,test="BR")
dev.off()




##In order to create supplemental table 1, IHC markers were included in a linear model to determine effect on response
##Overall fit pvalue from https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

clin$PRvNot <- ifelse(clin$BestResponse=="PR",1,0)

besttab3 <- do.call(rbind,lapply(c("PRvNot"),function(x) {
  temptab <- do.call(rbind,lapply(c("PDL1, % Positive Cells","CD8, % Positive Cells","CD68, % Positive Cells","FOXP3, % Positive Cells","PD-1, % Positive Cells"),function(y) {
    temp <- clin
    colnames(temp) <- gsub(x,"Group",gsub(y,"Test",colnames(temp)))
    temptab2 <- do.call(rbind,lapply(c("Baseline","C2"),function(time) {
      temp <- temp[temp$`Sample Timepoint` == time & !is.na(temp$`Sample Timepoint`),]
      temp <- temp[!is.na(temp$Group) & !is.na(temp$Test),] 
      tests <- unique(temp$Group[!(temp$Group %in% c("Other","PDbyNTL"))])
      aux <- plyr::count(temp$Group)
      tests <- tests[tests %in% aux$x[aux$freq >= 3]]
      modsum <- summary(lm(Test ~ Cohort + Group,data=temp))
      modtemp <- (lm(Test ~ Cohort + Group,data=temp))
      p <- modsum$coefficients[length(rownames(modsum$coefficients)),4]
      totalp <- lmp(modtemp)
      data.frame(Group=x,Test=y,Time=time,modelp=totalp,Testp=p,stringsAsFactors = FALSE)
    }))
    temptab2
  }))
  temptab
}))




besttab4 <- do.call(rbind,lapply(c("PRvNot"),function(x) {
  temptab <- do.call(rbind,lapply(c("PDL1, % Positive Cells","CD8, % Positive Cells","CD68, % Positive Cells","FOXP3, % Positive Cells","PD-1, % Positive Cells"),function(y) {
    temp <- clin
    colnames(temp) <- gsub(x,"Group",gsub(y,"Test",colnames(temp)))
    temp <- temp[!is.na(temp$Group) & !is.na(temp$Test),]
    temp2 <- spread(temp[temp$`Sample Timepoint` != "Progression" & !is.na(temp$`Sample Timepoint`),c("CMO Patient ID","Test","Cohort","Sample Timepoint","Group")],`Sample Timepoint`,Test)
    temp2$num <- rowSums(!is.na(temp2[,c("Baseline","C2")]))
    temp2 <- temp2[temp2$num == 2,]
    temp2$Change <- temp2$C2 - temp2$Baseline
    modsum <- summary(lm(Change ~ Cohort + Group,data=temp2))
    modtemp <- (lm(Change ~ Cohort + Group,data=temp2))
    p <- modsum$coefficients[length(rownames(modsum$coefficients)),4]
    totalp <- lmp(modtemp)
    data.frame(Group=x,Test=y,Time="Subtract",modelp=totalp,Testp=p,stringsAsFactors = FALSE)
  }))
  temptab
}))





