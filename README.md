More information can be found in: 
https://www.biorxiv.org/content/10.64898/2026.01.25.701577v1
Marques, R. M., Santos, C., Pai, H., Patto, M. C. V., Kamoun, S., & Kourelis, J. (2026). NLR immune receptors can exhibit tissue-specific expression patterns across legume species. bioRxiv, 2026-01.
Maravilha Marques, R., Santos, C., Pai, H., Vaz Patto, M. C., Kamoun, S., & Kourelis, J. (2026). Legume NLR immune receptors exhibit tissue-specific expression patterns across species - scripts and data. Zenodo. https://doi.org/10.5281/zenodo.18173999

# NLR_Classes_MaravilhaR_012026.R
The NLR_Classes_MaravilhaR R script takes 3 files: an excel file labelled “NLR Classes Legumes w RefPlantNLR.xlsx”, "RefPlantNLR.xlsx" and "Fabales_LOCUS_09012023.tsv".
 
The R script was heavily based on Jiorgos Kourelis' scripts. It takes: 
1) NLR Classes Legumes w RefPlantNLR.xlsx - the NLRtracker output aligned to the RefPlantNLR Supplemental datasets Supplemental_dataset_01_RefPlantNLR_v.20210712_481_AA, Supplemental_dataset_09_RefPlantNLR_vv.20210712_481_AA_NBARC_deduplictated, Supplemental_dataset_13_RefPlantNLR_v.20210712_481_Other-NBARC and TN_OTHER
2) RefPlantNLR.xlsx - The RefPlantNLR datasets alone for removal purposes
3) Fabales_LOCUS_09012023.tsv - All the loci of the employed Fabaceae species with the respective sequence names, for deduplication purposes
And generates graphs for comparative NLRome analysis.

# DESEQ_Single_Pipeline_RitaMaravilha_012026.R
Takes count data from Stringtie, performs differential expression analysis and contains data visualization sections

# Trim_BLASTp_Results
Takes BLASTp results, filters them and contains data visualization sections

Kourelis, J., Sakai, T., Adachi, H., & Kamoun, S. (2021). RefPlantNLR is a comprehensive collection of experimentally validated plant disease resistance proteins from the NLR family. PLoS Biology, 19(10), e3001124.
