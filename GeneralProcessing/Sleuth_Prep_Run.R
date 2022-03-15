##Script to Run Sleuth Prep for RNAseq Analysis
##Relevant to Figure 3 and supplemental figure 2
##Author - Allison L Richards


library(data.table)
library(reshape2)
library(sleuth)
library(argparse)


parser=ArgumentParser()
parser$add_argument('-d','--designmatrix',type='character',help='Design Matrix for Sleuth')
parser$add_argument('-t','--type',type='character',help="Either 'Gene' or 'Transcript'")
parser$add_argument('-o','--outputfile',type='character',help="Full file name for Rdata output")

args=parser$parse_args()

designpath = args$designmatrix
type = args$type
outpath = args$outputfile

df2 <- fread(designpath,data.table=FALSE)


anno <- fread("Hugo_ENST_ensembl75.txt",data.table=FALSE) #Annotation file pairing ENST, ENSG and Hugo Symbol of Ensembl version 75

if (type == "Transcript" ) {
  sc <- sleuth_prep(df2,extra_bootstrap_summary=TRUE,read_bootstrap_tpm=TRUE,target_mapping=anno)
  sleuth_save(sc,outpath)
} else if (type == "Gene") {
  scg <- sleuth_prep(df2,extra_bootstrap_summary=TRUE,read_bootstrap_tpm=TRUE,target_mapping=anno,aggregation_column="hugo",gene_mode=TRUE)
  sleuth_save(scg,outpath)
} else {
  cat('You\'ve chosen a type that isn\'t available')
}
