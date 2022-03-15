# Executable -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(stringr)
  library(dplyr)
})

if (!interactive()) {
  
  parser=ArgumentParser()
  parser$add_argument('-a', '--arriba', type='character', help='file path to arriba output')
  parser$add_argument('-fc', '--fusioncatcher', type='character', help='file path to fusioncatcher output')
  parser$add_argument('-o', '--output', type='character', help='output path')
  args=parser$parse_args()
  
  arriba.filename = args$arriba
  fc.filename = args$fusioncatcher
  outputdir = args$output
  

  # Renaming to ensure consistency ------------------------------------------
  arriba.data <- fread(arriba.filename) %>% mutate(SR = split_reads1 + split_reads2) %>%
    # standardized columns
    setnames(c('#gene1','gene2','breakpoint1','breakpoint2','SR','discordant_mates','reading_frame','fusion_transcript'),
             c('gene1','gene2','breakpoint1','breakpoint2','SR','PR','reading_frame','fusion_transcript')) %>%
    select(-c(direction1,direction2,`strand1(gene/fusion)`,`strand1(gene/fusion)`,`strand2(gene/fusion)`,
              split_reads1,split_reads2,coverage1,coverage2,closest_genomic_breakpoint1,
              closest_genomic_breakpoint2,filters,read_identifiers)) %>% data.table()
  fc.data <- fread(fc.filename) %>% filter(!grepl('no_protein',Fusion_description)) %>%
    mutate(`Fusion_point_for_gene_1(5end_fusion_partner)` = gsub(':\\+|:-','',`Fusion_point_for_gene_1(5end_fusion_partner)`),
           `Fusion_point_for_gene_2(3end_fusion_partner)` = gsub(':\\+|:-','',`Fusion_point_for_gene_2(3end_fusion_partner)`)) %>%
    # standardized columns
    setnames(c('Gene_1_symbol(5end_fusion_partner)','Gene_2_symbol(3end_fusion_partner)',
               'Fusion_point_for_gene_1(5end_fusion_partner)','Fusion_point_for_gene_2(3end_fusion_partner)',
               'Spanning_unique_reads','Spanning_pairs','Predicted_effect','Fusion_sequence'),
             c('gene1','gene2','breakpoint1','breakpoint2','SR','PR','reading_frame','fusion_transcript')) %>%
    select(gene1,gene2,Fusion_description,breakpoint1,breakpoint2,SR,PR,reading_frame,fusion_transcript) %>% data.table()
  
  # merging -----------------------------------------------------------------
  merge(arriba.data,fc.data,by = c('gene1','gene2','breakpoint1','breakpoint2'),suffixes = c('.arriba','.fc')) %>%
    setnames('confidence','arriba.confidence') %>% 
    # call confidence and pizzly
    mutate(call.confidence = case_when(
      is.na(reading_frame.fc) & !is.na(reading_frame.arriba) ~ 'Arriba',
      !is.na(reading_frame.fc) & is.na(reading_frame.arriba) ~ 'Fusion catcher',
      !is.na(reading_frame.fc) & !is.na(reading_frame.arriba) ~ 'Arriba;Fusion catcher',
      is.na(reading_frame.fc) & is.na(reading_frame.arriba) ~ 'None Called'
    ))  %>%
    select(gene1,gene2,type,call.confidence,breakpoint1,breakpoint2,site1,site2,
           PR.arriba,SR.arriba,PR.fc,SR.fc,arriba.confidence,reading_frame.arriba,reading_frame.fc,
           fusion_transcript.arriba,fusion_transcript.fc,peptide_sequence,Fusion_description) %>%
    data.table()  -> arriba.fc.merged
  if(any(grepl('None Called',unlist(arriba.fc.merged)))){
    stop('Some of the calls are weird :/')
  }
  
  print(paste0('Writing to: ',outputdir))
  write.table(arriba.fc.merged,outputdir,quote = F,sep = '\t',row.names = F)
}
