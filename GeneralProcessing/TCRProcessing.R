##Script to Create Bash Scripts for Processing TCR-seq
##Relevant to supplemental figure 2
##Author - Allison L Richards



##This script is based off Mixcr and VDJtools manuals


#Load in fastqs from dbGap
fast <- list.files("XXX",full.names=TRUE)

bashrun <- do.call(rbind,lapply(fast[grepl("R1",fast)],function(r1) {
  samp <- strsplit(r1,"/|_")[[1]]
  r2 <- fast[grepl(samp,fast) & grepl("R2",fast)]
  one <- paste0("mkdir ",samp)
  two <- paste0("bsub -n 2 -R rusage[mem=60] -W 10:00 -J mixcr",samp," 'java -jar mixcr-3.0.11/mixcr.jar align -OsaveOriginalReads=true -OvjAlignmentOrder=JThenV -s hs -r ",samp,"/mixcr.log -t 1 ",r1," ",r2," ",samp,"/",samp,"_mixcr_ABDG.vdjca'")
  rbind(one,two)
}))
write.table(bashrun,"BashAlignMixcr.sh",quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")


##That was quick, now run assembly
files2 <- list.files("XXX",pattern="ABDG.vdjca",recursive = TRUE,full.names = TRUE)

bashass <- do.call(rbind,lapply(files2,function(x) {
  loc <- strsplit(x,"/")[[1]]
  samp <- strsplit(loc[grep("vdjca",loc)],"_")[[1]][1]
  loc <- paste(loc[1:(length(loc)-1)],collapse="/")
  two <- paste0("bsub -n 2 -R rusage[mem=60] -W 10:00 -J mixcr2",samp," 'java -jar mixcr-3.0.11/mixcr.jar assemble -OcloneClusteringParameters.searchParameters=twoMismatches -f -r ",samp,"/mixcr.log -t 1 ",x," ",loc,"/",samp,"_mixcr.clns'")
  rbind(two)
}))
write.table(bashass,"BashAssembleMixcr.sh",quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")


##Next step is exporting the clones
files3 <- list.files("XXX",pattern="mixcr.clns",recursive = TRUE,full.names = TRUE)

bashexp <- do.call(rbind,lapply(files3,function(x) {
  loc <- strsplit(x,"/")[[1]]
  samp <- strsplit(loc[grep("clns",loc)],"_")[[1]][1]
  loc <- paste(loc[1:(length(loc)-1)],collapse="/")
  two <- paste0("bsub -n 1 -R rusage[mem=10] -W 1:00 -J mixcr3",samp," 'java -jar mixcr-3.0.11/mixcr.jar exportClones ",x," ",loc,"/",samp,"_mixcr_clones.txt'")
  rbind(two)
}))
write.table(bashexp,"BashExportMixcr.sh",quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")


##OK, now that this is all done, we can move onto VDJ tools which needs all the files in the same location and a metadata file
files4 <- list.files("XXX",pattern="clones.txt",recursive = TRUE,full.names = TRUE)

bashmove <- do.call(rbind,lapply(files4,function(x) {
  start <- x
  loc <- strsplit(x,"/")[[1]]
  samp <- strsplit(loc[grep("clones",loc)],"_")[[1]][1]
  end <- paste0("XXY",samp,"_",batch,"_mixcr_clones.txt")
  two <- paste0("scp ",start," ",end)
  rbind(two)
}))
write.table(bashmove,"BashReorganizeMixcr.sh",quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")



#Pull together metadata file
meta <- fread("PhenotypeTableFordbGap_052421.txt",data.table=FALSE)
#Pieces needed: file_name, sample_id, batch_tcr
files5 <- list.files("XXX",pattern="mixcr_clones")

colnames(meta) <- gsub(" |[?]|[(]|[)]|#|[.]|[/]|%|,|&|\\^|","",colnames(meta))
meta2 <- do.call(rbind,lapply(rownames(meta),function(x) {
  temp <- meta[x,]
  sample <- temp$SampleID
  thing <- files5[grepl(sample,files5)]
  if(length(thing) > 0 & !is.na(sample)) {
    id <- gsub("_mixcr_clones[.]txt","",thing)
    tcr <- strsplit(id,"_")[[1]]
    tcr <- paste0(tcr[length(tcr) - 1],"_",tcr[length(tcr)])
    temp2 <- cbind(data.frame(file_name=thing,sample_id=id,batchtcr=tcr,stringsAsFactors = FALSE),temp)
    temp2
  }
}))
write.table(meta2,"MetaFileForConvert.txt",row.names=FALSE,sep="\t",quote=FALSE)



##Followed by following scripts in terminal
# java -jar vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S MiXcr -m MetaFileForConvert.txt Converted_ >Convert.o 2>Convert.e
# java -jar vdjtools-1.2.1/vdjtools-1.2.1.jar CalcBasicStats -m metadata.txt Output >CalcBasStat.o 2>CalcBasStat.e


