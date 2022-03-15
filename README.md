# ImmunoSarc



## Figures
Figures directory includes source data and code necessary to recreate all figures and supplemental figures.

## GeneralProcessing
GeneralProcessing directory includes scripts used for processing and overall analysis of the data.

Scripts to process RNA-seq are found in GeneralProcessing and include:

```
RunSTAR.sh
RunKallisto.sh
RunFusionCatcher.sh
RunArriba.sh
CombineFusionCalls.R
```

Further information on references and versions can be found in the Methods section of the manuscript.

Other GeneralProcessing scripts require output from Tempo (for WES), Kallisto and above scripts (for RNA-seq) and MiXCR and VDJtools (for TCR-seq), as well as patient and sample data that are provided in the Figures directory or in supplemental material.

SetUpData.R - Script to organize data and identify Immune groups. This script also includes a heatmap similar to Figure 2D and 2E except missing a handful of tracks.

CompareVariablesToResponseAndCohort.R - Script to recreate plots and p-values for figures comparing variables (e.g. IHC values, TCR diversity, TMB, etc) to sample cohort or patient response. Relevant to Figure 2, and supplemental figures 2, 3, 5, and 6.

KMPlots.R - Script on Kaplan-Meier plots and how we compare survival. Relevant to Figure 2, and supplemental figures 1, 5 and 6.

Sleuth_Prep_Run.R - Script for running sleuth_prep.

DifferentialExpressionAnalysis.R - Script describing the differential expression analysis, gene set enrichment analysis and modeling included in Figure 3 and supplemental figure 3.

## Contact
E-mail any questions to richara4@mskcc.org.
