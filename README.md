# NLR_Classes_MaravilhaR
The NLR_Classes_MaravilhaR R script takes 3 files: an excel file labelled “NLR Classes Legumes w RefPlantNLR.xlsx”, "RefPlantNLR.xlsx" and "Fabales_LOCUS_09012023.tsv".
 
The R script was heavily based on Jiorgos Kourelis' scripts. It takes: 
1) NLR Classes Legumes w RefPlantNLR.xlsx - the NLRtracker output aligned to the RefPlantNLR Supplemental datasets Supplemental_dataset_01_RefPlantNLR_v.20210712_481_AA, Supplemental_dataset_09_RefPlantNLR_vv.20210712_481_AA_NBARC_deduplictated, Supplemental_dataset_13_RefPlantNLR_v.20210712_481_Other-NBARC and TN_OTHER
2) RefPlantNLR.xlsx - The RefPlantNLR datasets alone for removal purposes
3) Fabales_LOCUS_09012023.tsv - All the loci of the employed Fabaceae species with the respective sequence names, for deduplication purposes
And generates graphs for comparative NLRome analysis.

# DESEQ_Single_Pipeline_RitaMaravilha_2024
Takes count data from Stringtie, performs differential expression analysis and contains data visualization sections

# Trim_BLASTp_Results
Takes BLASTp results, filters them and contains data visualization sections
