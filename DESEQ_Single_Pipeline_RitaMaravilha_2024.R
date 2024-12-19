pkgs <- rownames(installed.packages())
if(!"dplyr" %in% pkgs) install.packages("dplyr")
if(!"ggplot2" %in% pkgs) install.packages("ggplot2")
if(!"stringr" %in% pkgs) install.packages("stringr")
if(!"tidyverse" %in% pkgs) install.packages("tidyverse")
if(!"readxl" %in% pkgs) install.packages("readxl")
if(!"VennDiagram" %in% pkgs) install.packages("VennDiagram")
if(!"pheatmap" %in% pkgs) install.packages("pheatmap")
if(!"RColorBrewer" %in% pkgs) install.packages("RColorBrewer")
if(!"PoiClaClu" %in% pkgs) install.packages("PoiClaClu")
# Require necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
#BiocManager::install("IsoformSwitchAnalyzeR")
library(DESeq2)
library(edgeR)
library(apeglm)
#library(IsoformSwitchAnalyzeR)
#library(dplyr)
library(ggplot2)
#library(stringr)
library(tidyverse)
library(readxl)
library(VennDiagram)
library(pheatmap)
library("RColorBrewer")
library("PoiClaClu")
setwd("D:/Dropbox (KamounLab)/Kamounity folder/Rita/Tissue_Spec/Counts_2024")
setwd("C:/Users/Rita/Dropbox (KamounLab)/Kamounity folder/Rita/Tissue_Spec/Counts_2024")
setwd("D:/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/Counts_2024")
setwd("C:/Users/RitaML/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/Counts_2024")

##### Creating BUSCO id counts #####
setwd("D:/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/BUSCO PCA")

# before reading in the tabular BUSCO data, remove the # from the third line
# and replace Busco id by Busco_id (the space isn't a good character)
BUSCO_Car <- read.table("BUSCO_Car_table.tabular", 
                        sep="\t", # separator is a tab
                        header=TRUE, 
                        fill=TRUE, # Write NA where no value is provided
                        quote="") # ignores " and '
BUSCO_Mtr <- read.table("BUSCO_Mtr_table.tabular", sep="\t", header=TRUE, 
                        fill=TRUE, quote="")
BUSCO_Lsa <- read.table("BUSCO_Lsa_table.tabular", sep="\t", header=TRUE, 
                        fill=TRUE, quote="")
BUSCO_Psa <- read.table("BUSCO_Psa_table.tabular", sep="\t", header=TRUE, 
                        fill=TRUE, quote="")
BUSCO_Lcu <- read.table("BUSCO_Lcu_table.tabular", sep="\t", header=TRUE, 
                        fill=TRUE, quote="")
BUSCO_Pvu <- read.table("BUSCO_Pvu_table.tabular", sep="\t", header=TRUE, 
                        fill=TRUE, quote="")
BUSCO_Gma <- read.table("BUSCO_Gma_table.tabular", sep="\t", header=TRUE, 
                        fill=TRUE, quote="")

# the differences in BUSCO number are due to duplications, so let's only keep 
#the max score. The total number of BUSCOs should be 5366
BUSCO_Car <- BUSCO_Car %>%
  mutate(Species = "Cicer arietinum") %>% # add the Species name
  group_by(Busco_id) %>% # group by the BUSCO id
  arrange(desc(Score)) %>% # keep the one with the highest score
  slice(1) #slice the duplicates
BUSCO_Gma <- BUSCO_Gma %>% mutate(Species = "Glycine max") %>% 
  group_by(Busco_id) %>% arrange(desc(Score)) %>% slice(1)
BUSCO_Lcu <- BUSCO_Lcu %>% mutate(Species = "Lens culinaris") %>% 
  group_by(Busco_id) %>% arrange(desc(Score)) %>% slice(1)
BUSCO_Lsa <- BUSCO_Lsa %>% mutate(Species = "Lathyrus sativus") %>% 
  group_by(Busco_id) %>% arrange(desc(Score)) %>% slice(1)
BUSCO_Mtr <- BUSCO_Mtr %>% mutate(Species = "Medicago truncatula") %>% 
  group_by(Busco_id) %>% arrange(desc(Score)) %>% slice(1)
BUSCO_Psa <- BUSCO_Psa %>% mutate(Species = "Pisum sativum") %>% 
  group_by(Busco_id) %>% arrange(desc(Score)) %>% slice(1)
BUSCO_Pvu <- BUSCO_Pvu %>% mutate(Species = "Phaseolus vulgaris") %>% 
  group_by(Busco_id) %>% arrange(desc(Score)) %>% slice(1)

  #[, -c("Status", "")]
#BUSCO_all <- list(BUSCO_Car, BUSCO_Gma, BUSCO_Lcu, BUSCO_Lsa, BUSCO_Mtr, 
#                  BUSCO_Psa, BUSCO_Pvu) %>% # remove unimportant columns
#  reduce(left_join, by = c("Busco_id", "Description", "OrthoDB.url"))
BUSCO_all <- rbind(BUSCO_Car, BUSCO_Gma, BUSCO_Lcu, BUSCO_Lsa, BUSCO_Mtr, 
BUSCO_Psa, BUSCO_Pvu)

# have a gene id corresponded to the transcript id for deduplicating
BUSCO_all <- BUSCO_all %>%
  ungroup() %>%
  mutate(gene_id = case_when(
    Species == "Lathyrus sativus" ~ str_sub(BUSCO_all$Sequence, end = -4),
    Species == "Glycine max" ~ str_sub(BUSCO_all$Sequence, end = -5),
    Species != "Medicago truncatula" ~ str_sub(BUSCO_all$Sequence, end = -3),
    TRUE ~ Sequence
  ))

write.csv(BUSCO_all, "BUSCO_ids_Tissue_Spec.csv")
setwd("D:/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/Counts_2024")
write.csv(BUSCO_all, "BUSCO_ids_Tissue_Spec.csv")
BUSCO_all <- read.csv("BUSCO_ids_Tissue_Spec.csv", row.names = 1)


##### make a file with all expression data from all species #####
# First, let's keep only the transcripts found with Stringtie
# Car
counts_Car <- read.csv(file="gene_count_matrix_Car.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Car$Gene_id <- sapply(strsplit(counts_Car$gene_id, "\\|"), function(x) tail(x, 1))
write.csv(counts_Car, "gene_count_matrix_Car_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
counts_Car <- counts_Car[grep("^cicar.CDCFrontier.gnm1.ann1.Ca_", counts_Car$Gene_id),]
counts_Car <- counts_Car %>%
  mutate(Species = "Cicer arietinum")
write.csv(counts_Car, "gene_count_matrix_Car_without_MSTRG.csv")

# Gma
counts_Gma <- read.csv(file="gene_count_matrix_Gma.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Gma$Gene_id <- sapply(strsplit(counts_Gma$gene_id, "\\|"), function(x) tail(x, 1))
# Replace gene-GLYMA_ by Glyma. and remove the v4 at the end
counts_Gma <- counts_Gma %>%
  mutate(Gene_id = str_replace_all(Gene_id, "gene-GLYMA_", "Glyma.")) 
write.csv(counts_Gma, "gene_count_matrix_Gma_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
counts_Gma <- counts_Gma[grep("^Glyma.", counts_Gma$Gene_id),]
counts_Gma <- counts_Gma %>%
  mutate(Gene_id = str_sub(Gene_id, end = -3)) %>% #remove the last characters to match the seqname
  mutate(Species = "Glycine max")
write.csv(counts_Gma, "gene_count_matrix_Gma_without_MSTRG.csv")

# Lcu
counts_Lcu <- read.csv(file="gene_count_matrix_Lcu.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Lcu$Gene_id <- sapply(strsplit(counts_Lcu$gene_id, "\\|"), function(x) tail(x, 1))
write.csv(counts_Lcu, "gene_count_matrix_Lcu_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
counts_Lcu <- counts_Lcu[grep("^Lcu.2RBY.", counts_Lcu$Gene_id),]
counts_Lcu <- counts_Lcu %>%
  mutate(Species = "Lens culinaris")
write.csv(counts_Lcu, "gene_count_matrix_Lcu_without_MSTRG.csv")

# Lsa
counts_Lsa <- read.csv(file="gene_count_matrix_Lsa.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Lsa$Gene_id <- sapply(strsplit(counts_Lsa$gene_id, "\\|"), function(x) tail(x, 1))
write.csv(counts_Lsa, "gene_count_matrix_Lsa_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
counts_Lsa <- counts_Lsa[grep("^g", counts_Lsa$Gene_id),]
counts_Lsa <- counts_Lsa %>%
  mutate(Species = "Lathyrus sativus")
write.csv(counts_Lsa, "gene_count_matrix_Lsa_without_MSTRG.csv")

# Mtr
counts_Mtr <- read.csv(file="gene_count_matrix_Mtr.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Mtr$Gene_id <- sapply(strsplit(counts_Mtr$gene_id, "\\|"), function(x) tail(x, 1))
# make the sequence name match the protein ID
counts_Mtr <- counts_Mtr %>%
  mutate(Gene_id = str_replace_all(Gene_id, "MtrunA17Chr0c0", "MtrunA17_Chr")) %>%
  mutate(Gene_id = str_replace_all(Gene_id, "R", "g")) %>%
  mutate(Gene_id = str_replace_all(Gene_id, "MtrunA17Chr", "MtrunA17_Chr"))
write.csv(counts_Mtr, "gene_count_matrix_Mtr_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
counts_Mtr <- counts_Mtr[grep("^MtrunA17", counts_Mtr$Gene_id),]
counts_Mtr <- counts_Mtr %>%
  mutate(Species = "Medicago truncatula")
write.csv(counts_Mtr, "gene_count_matrix_Mtr_without_MSTRG.csv")

# Psa
counts_Psa <- read.csv(file="gene_count_matrix_Psa.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Psa$Gene_id <- sapply(strsplit(counts_Psa$gene_id, "\\|"), function(x) tail(x, 1))
write.csv(counts_Psa, "gene_count_matrix_Psa_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
counts_Psa <- counts_Psa[grep("^LOC", counts_Psa$Gene_id),]
counts_Psa <- counts_Psa %>%
  mutate(Species = "Pisum sativum")
write.csv(counts_Psa, "gene_count_matrix_Psa_without_MSTRG.csv")
# this is a bit annoying but I'll need to match the locus to the seqname

# Pvu
counts_Pvu <- read.csv(file="gene_count_matrix_Pvu.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Pvu$Gene_id <- sapply(strsplit(counts_Pvu$gene_id, "\\|"), function(x) tail(x, 1))
# Replace gene-PHAVU_ by phavu.G19833.gnm2.ann1.Phvul. and remove the g at the end
counts_Pvu <- counts_Pvu %>%
  mutate(Gene_id = str_replace_all(Gene_id, "gene-PHAVU_", "phavu.G19833.gnm2.ann1.Phvul."))
write.csv(counts_Pvu, "gene_count_matrix_Pvu_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
counts_Pvu <- counts_Pvu[grep("^phavu.G19833.gnm2.ann1.Phvul.", counts_Pvu$Gene_id),]
counts_Pvu <- counts_Pvu %>%
  mutate(Gene_id = str_sub(Gene_id, end = -2)) %>% #remove the last characters to match the seqname
  mutate(Species = "Phaseolus vulgaris")
write.csv(counts_Pvu, "gene_count_matrix_Pvu_without_MSTRG.csv")

# Make a counts file with all species
counts_all <- counts_Car
# change column names to Leaf and Root instead of sample names
colnames(counts_all)[2:7] <- c("Leaf 1", "Leaf 2", "Leaf 3", 
                               "Root 1", "Root 2", "Root 3") 
counts_all <- counts_all %>%
  full_join(counts_Gma, by=c("Species","Gene_id","Leaf 1"="GmaL1", "Leaf 2"="GmaL2", "Leaf 3"="GmaL3",
                             "Root 1"="GmaR1", "Root 2"="GmaR2", "Root 3"="GmaR3", "gene_id")) %>%
  full_join(counts_Lcu, by=c("Species","Gene_id","Leaf 1"="LcuL1", "Leaf 2"="LcuL2", "Leaf 3"="LcuL3",
                             "Root 1"="LcuR1", "Root 2"="LcuR2", "Root 3"="LcuR3", "gene_id")) %>%
  full_join(counts_Lsa, by=c("Species","Gene_id","Leaf 1"="LsaL1", "Leaf 2"="LsaL2", "Leaf 3"="LsaL3",
                             "Root 1"="LsaR1", "Root 2"="LsaR2", "Root 3"="LsaR3", "gene_id")) %>%
  full_join(counts_Mtr, by=c("Species","Gene_id","Leaf 1"="MtrL1", "Leaf 2"="MtrL2", "Leaf 3"="MtrL3",
                             "Root 1"="MtrR1", "Root 2"="MtrR2", "Root 3"="MtrR3", "gene_id")) %>%
  full_join(counts_Psa, by=c("Species","Gene_id","Leaf 1"="PsaL1", "Leaf 2"="PsaL2", "Leaf 3"="PsaL3",
                             "Root 1"="PsaR1", "Root 2"="PsaR2", "Root 3"="PsaR3", "gene_id")) %>%
  full_join(counts_Pvu, by=c("Species","Gene_id","Leaf 1"="PvuL1", "Leaf 2"="PvuL2", "Leaf 3"="PvuL3",
                             "Root 1"="PvuR1", "Root 2"="PvuR2", "Root 3"="PvuR3", "gene_id"))
# 28270 + 52820 + 55845 + 31721 + 293348 + 62857 + 28131 = 552992
write.csv(counts_all, "All_species_counts_without_MSTRG_not_filtered_or_normalized.csv")

#### BUSCO PCA ######
BUSCO_all <- read.csv("BUSCO_ids_Tissue_Spec.csv", row.names = 1)
counts_all <- read.csv("All_species_counts_without_MSTRG_not_filtered_or_normalized.csv", row.names = 1)
LOCUS <- read.tsv("Fabales_LOCUS_29022024.tsv")
LOCUS_tissue_spec <- LOCUS %>%
  mutate(Genome = case_when(
    Organism == "Pisum sativum ZW6" ~ "CAAS_Psat_ZW6_1.0",
    TRUE ~ Genome
  )) %>%
  filter(Genome %in% c("cicar.CDCFrontier.gnm1.ann1.nRhs.protein", 
                       "Glycine max Williams 82", 
                       "Lens culinaris Redberry v2.0",
                         "gp_hifiasm_assembly.v1.0", "MtrunA17r5.0-ANR-EGN-r1.9",
                         "CAAS_Psat_ZW6_1.0", "phavu.G19833.gnm2.ann1.PB8d.protein"))
LOCUS_tissue_spec <- LOCUS_tissue_spec %>%
  mutate(Locus_id = case_when(
    Organism == "Pisum sativum ZW6" ~ Locus,
    Organism != "Medicago truncatula" ~ gsub("[.][^.]+$", "", seqname), # remove everything after the last .
    TRUE ~ seqname
  ))
LOCUS_tissue_spec <- LOCUS_tissue_spec %>%
  # copy the info in coded_by when the phaseolus vulgaris or cercis canadensis
  #locus is incomplete
  mutate(Locus = case_when(
    Locus == "Phvul" ~ str_sub(LOCUS_tissue_spec$coded_by, end = -3),
    TRUE ~ Locus))

counts_all_LOCUS <- counts_all %>%
  left_join(LOCUS_tissue_spec, by= c("Gene_id" = "Locus_id")) 

# Pisum sativum has LOC instead of XP in Gene_id, so let's correct that before
#joining dataframes
counts_all_LOCUS <- counts_all_LOCUS %>%
  drop_na(Organism) %>%
  mutate(Gene_id = case_when(
    Species == "Pisum sativum" ~ seqname,
    TRUE ~ Gene_id
  ))

BUSCO_counts_all_2 <- left_join(BUSCO_all, counts_all_LOCUS, by=c("Species", 
                                                          "gene_id"="Gene_id"
                                                          ))

setwd("D:/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/Counts_2024")
write.csv(BUSCO_counts_all_2, "BUSCO_counts_all_duplicated.csv")
BUSCO_counts_all_2 <- read.csv("BUSCO_counts_all_duplicated.csv", row.names=1)
BUSCO_counts_all_2 <- BUSCO_counts_all_2 %>%
  group_by(Busco_id, Genome) %>%
  arrange(desc(seqname)) %>%
  arrange(desc(Length.x)) %>% # keep the larger
  slice(1) #37170/7 = 5310 BUSCOs

count(BUSCO_counts_all_2$Status=="Missing") # 1189

setwd("D:/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/Counts_2024")
write.csv(BUSCO_counts_all_2, "BUSCO_counts_all.csv")

# Let's try to make a PCA plot
# Load metadata
metadata <- read_xlsx("RNASeq_Tissue_metadata.xlsx", 
                      #header = TRUE, sep = ",", 
)
metadata$Tissue <- relevel(factor(metadata$Tissue), ref = 'Root') #make root as ref
metadata$Species <- relevel(factor(metadata$Species), ref = 'Medicago truncatula')

# Prepare the counts file
BUSCO_counts_all_counts <- BUSCO_counts_all_2 %>%
  ungroup() %>%
  select(Busco_id, Leaf.1, Leaf.2, Leaf.3, Root.1, Root.2, Root.3, Species)
# Make a counts file per species and then merge the 7
Psa_BUSCO_counts <- BUSCO_counts_all_counts %>%
  filter(Species=="Pisum sativum")
colnames(Psa_BUSCO_counts)=c("Busco_id","PsaL1","PsaL2","PsaL3","PsaR1","PsaR2",
                             "PsaR3","Species")
Psa_BUSCO_counts$Species <- NULL
Gma_BUSCO_counts <- BUSCO_counts_all_counts %>%
  filter(Species=="Glycine max")
colnames(Gma_BUSCO_counts)=c("Busco_id","GmaL1","GmaL2","GmaL3","GmaR1","GmaR2",
                             "GmaR3","Species")
Gma_BUSCO_counts$Species <- NULL
Lcu_BUSCO_counts <- BUSCO_counts_all_counts %>%
  filter(Species=="Lens culinaris")
colnames(Lcu_BUSCO_counts)=c("Busco_id","LcuL1","LcuL2","LcuL3","LcuR1","LcuR2",
                             "LcuR3","Species")
Lcu_BUSCO_counts$Species <- NULL
Mtr_BUSCO_counts <- BUSCO_counts_all_counts %>%
  filter(Species=="Medicago truncatula")
colnames(Mtr_BUSCO_counts)=c("Busco_id","MtrL1","MtrL2","MtrL3","MtrR1","MtrR2",
                             "MtrR3","Species")
Mtr_BUSCO_counts$Species <- NULL
Car_BUSCO_counts <- BUSCO_counts_all_counts %>%
  filter(Species=="Cicer arietinum")
colnames(Car_BUSCO_counts)=c("Busco_id","CarL1","CarL2","CarL3","CarR1","CarR2",
                             "CarR3","Species")
Car_BUSCO_counts$Species <- NULL
Lsa_BUSCO_counts <- BUSCO_counts_all_counts %>%
  filter(Species=="Lathyrus sativus")
colnames(Lsa_BUSCO_counts)=c("Busco_id","LsaL1","LsaL2","LsaL3","LsaR1","LsaR2",
                             "LsaR3","Species")
Lsa_BUSCO_counts$Species <- NULL
Pvu_BUSCO_counts <- BUSCO_counts_all_counts %>%
  filter(Species=="Phaseolus vulgaris")
colnames(Pvu_BUSCO_counts)=c("Busco_id","PvuL1","PvuL2","PvuL3","PvuR1","PvuR2",
                             "PvuR3","Species")
Pvu_BUSCO_counts$Species <- NULL

require(plyr)
BUSCO_counts_all_counts_pivot <- join_all(
  list(Psa_BUSCO_counts, Gma_BUSCO_counts, Lcu_BUSCO_counts, Mtr_BUSCO_counts,
       Car_BUSCO_counts, Lsa_BUSCO_counts, Pvu_BUSCO_counts), 
  by = "Busco_id")
write.csv(BUSCO_counts_all_counts_pivot, "BUSCO_counts_all_counts_pivot.csv")

# Filter rowmean > 5
counts_BUSCO <- read.csv("BUSCO_counts_all_counts_pivot.csv", row.names=2)
counts_BUSCO$X <- NULL
# I won't filter anything for now to maximize the information from BUSCOs
#counts$rowmean <- rowMeans(counts)
#counts_filtered <- subset(counts, rowmean > 5)
#counts_filtered <- counts_filtered[,c(1:6)]

# TMM Normalization
# make the DGEList:
# I need to remove the NAs first for the PCA to work. 
# This means I will only analyse the BUSCOs with orthologs in all species: 
#3782/5366. This is a limitation of this methodology.
counts_BUSCO <- counts_BUSCO %>%
  drop_na()
TMM <- DGEList(counts_BUSCO)

#/ calculate TMM normalization factors:
TMM <- calcNormFactors(TMM)

#/ get the normalized counts:
counts_filtered_TMM <- cpm(TMM, log=FALSE)
counts_filtered_TMM_log2 <- cpm(counts_BUSCO, log=T)
write.csv(as.data.frame(counts_filtered_TMM), file=paste("DESeq/Raw/deseq2_BUSCO_TMM.csv", sep = "_"))
counts_filtered_TMM <- read.csv(paste("DESeq/Raw/deseq2_BUSCO_TMM.csv", sep="_"), row.names = 1)
dds <- DESeqDataSetFromMatrix(countData=round(counts_filtered_TMM), 
                              colData = metadata,
                              design = ~ Species,
                              tidy = F)

DEdds <- DESeq(dds)
vsd <- vst(DEdds) #normalization
png(filename = paste("DESeq/Output/Leaf_vs_Root_BUSCO_PCA_Species.png", sep="_"), 
    width = 1250,height = 750, units = "px")
PCA <- plotPCA(vsd, intgroup = #c(
                 "Species" #"Species")
)
PCA + scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
  coord_fixed(ratio =2, xlim = NULL, ylim = NULL) + theme_minimal() + 
  geom_label(aes(label = name))
dev.off()

png(filename = paste("DESeq/Output/Leaf_vs_Root_BUSCO_PCA_Tissue.png", sep="_"), 
    width = 1250,height = 750, units = "px")
PCA <- plotPCA(vsd, intgroup = #c(
                 "Tissue"#"Species")
)
PCA + scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
  coord_fixed(ratio =2, xlim = NULL, ylim = NULL) + theme_minimal() + 
  geom_label(aes(label = name))
dev.off()
#ggsave("PCA 3782 BUSCO orthologs TMM.png", width=1250, height=720, units='px', dpi = 300)

##### BUSCO DEGs #####
dds <- DESeqDataSetFromMatrix(countData=round(counts_filtered_TMM), 
                              colData = metadata,
                              design = ~ Tissue,
                              tidy = F)

DEdds <- DESeq(dds)

resultsNames(DEdds) 
model.matrix(~ Tissue, data = metadata) # 0 = Root Sample, 1 = Leaf Sample

# apeglm shrinkage, apply normalization to the data
resNorm_Tissue <- lfcShrink(DEdds, coef="Tissue_Leaf_vs_Root", type="apeglm") # tissue_Leaf_vs_Root
# Visualize the data as a data frame
resNorm_Tissue_df <- as.data.frame(resNorm_Tissue)
# Get the rownames into a column named ID
resNorm_Tissue_df$ID = rownames(resNorm_Tissue_df)
rownames(resNorm_Tissue_df) <- NULL
# Provide a summary: how many genes are there and how many are up or down regulated?
# Add info to the excel file "DESeq2_Summary.xlsx" in the "Analysis" folder
summary(resNorm_Tissue)
write.csv(as.data.frame(resNorm_Tissue_df), file="DESeq/Sum/deseq2_BUSCO_All_genes.csv")
resNorm_Tissue_df <- read.csv("DESeq/Sum/deseq2_BUSCO_All_genes.csv")

# Subset for adjusted p value < 0.05
res <- subset(resNorm_Tissue, padj < 0.05) 
res_df <- as.data.frame(res)
# Subset for log2FC higher than 1 (absolute value)
res_deseq <- res[abs(res$log2FoldChange) > 1,]   
res_deseq_df <- as.data.frame(res_deseq)
# Save results; substitute species name by using Ctrl + Shift + Alt + M
write.csv(as.data.frame(res_deseq), file="DESeq/Sum/deseq2_BUSCO_DEGs.csv")
res_deseq <- read.csv("DESeq/Sum/deseq2_BUSCO_DEGs.csv")

##### BUSCO Heatmap #####
setwd("D:/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/Counts_2024")

library(heatmaply)
#counts_BUSCO <- read.csv("BUSCO_counts_all_counts_pivot.csv", row.names=2)
#counts_BUSCO$X <- NULL
counts_BUSCO <- read.csv(paste("DESeq/Raw/deseq2_BUSCO_TMM.csv", sep="_"), row.names = 1)
# Filter rowmean > 5
counts_BUSCO_filtered <- counts_BUSCO
#counts_BUSCO_filtered$rowmean <- rowMeans(counts_BUSCO_filtered, na.rm=TRUE)
#counts_BUSCO_filtered <- subset(counts_BUSCO_filtered, rowmean > 5000)
#counts_BUSCO_filtered <- counts_BUSCO_filtered[,c(1:42)] #remove the rowmean column

counts_BUSCO_filtered[counts_BUSCO_filtered==0] <- NA #to perform a log transformation, I need to remove 0s
counts_BUSCO_filtered <- counts_BUSCO_filtered %>%
  drop_na()
counts_BUSCO_filtered_more <- counts_BUSCO_filtered
counts_BUSCO_filtered_more$rowmean <- rowMeans(counts_BUSCO_filtered_more, na.rm=TRUE)
counts_BUSCO_filtered_more <- subset(counts_BUSCO_filtered_more, rowmean > 250)
counts_BUSCO_filtered_more <- counts_BUSCO_filtered_more[,c(1:42)] #remove the rowmean column

rc <- colorspace::rainbow_hcl(ncol(counts_BUSCO_filtered))
rc <- c(rep(c("#003b84","#ff7400","#00b67e","#488200","#FFBB00","#3f0051","#FF1100"),
            times=1, each=6))
heatmaply(log(counts_BUSCO_filtered),xlab = "samples", 
          fontsize_row = 8, fontsize_col = 7,
          colors = viridis(256, option = "mako"), grid_gap = 0.1, Rowv = FALSE,
          label_names = c("Tissue/Accession", "Gene", "Expression Group"),
          key.title = "log(Expression of BUSCOs)",
          col_side_colors = rc, file = "DESeq/Output/heatmaply_BUSCO_clustering_250_TMM_2.html"
          )
# 10203at72025 is a BUSCO with high leaf expression and low root expression in all
counts_BUSCO[counts_BUSCO==0] <- NA #to perform a log transformation, I need to remove 0s

counts_BUSCO <- counts_BUSCO %>%
  drop_na()
  
heatmaply(log(counts_BUSCO),xlab = "samples", 
          fontsize_row = 8, fontsize_col = 7,
          colors = viridis(256, option = "mako"), grid_gap = 0.1, Rowv = FALSE,
          label_names = c("Tissue/Accession", "Gene", "Expression Group"),
          key.title = "log(Expression of BUSCOs)",
          col_side_colors = rc
)

png(filename = "DESeq/Output/BUSCOs_QC_Heatmap_250_TMM.png", 
    width = 1920,height = 1080, units = "px")
pheatmap(log(counts_BUSCO_filtered_more), #clustering_distance_rows = poisd$dd, 
         #clustering_distance_cols = poisd$dd, 
         col = viridis(256, option = "mako"))
dev.off()
########### DESeq2 Pipeline ###########
# set working directory
setwd("D:/Dropbox (KamounLab)/Kamounity folder/Rita/Tissue_Spec/Counts_2024")
setwd("C:/Users/Rita/Dropbox (KamounLab)/Kamounity folder/Rita/Tissue_Spec/Counts_2024")
setwd("D:/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/Counts_2024")
setwd("C:/Users/RitaML/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/Counts_2024")

# Load metadata
metadata <- read_xlsx("RNASeq_Tissue_metadata.xlsx", 
                      #header = TRUE, sep = ",", 
)
metadata$Tissue <- relevel(factor(metadata$Tissue), ref = 'Root') #make root as ref

choice <- menu(c("Car", "Gma", "Lcu", "Lsa", "Mtr", "Psa", "Pvu"
                 #, "all"
                 ), graphics = TRUE , title="Which species do you want to analyse?")

if (choice == 1) {
  species <- "Car"
  species_full <- "Cicer_arietinum"
  col_data <- metadata[c(13:18),]
  cat("The species is Cicer arietinum")
} else if (choice == 2) {
  species <- "Gma"
  species_full <- "Glycine_max"
  col_data <- metadata[c(37:42),]
  cat("The species is Glycine max")
} else if (choice == 3) {
  species <- "Lcu"
  species_full <- "Lens_culinaris"
  col_data <- metadata[c(25:30),]
  cat("The species is Lens culinaris")
} else if (choice == 4) {
  species <- "Lsa"
  species_full <- "Lathyrus_sativus"
  col_data <- metadata[c(1:6),]
  cat("The species is Lathyrus sativus")
} else if (choice == 5) {
  species <- "Mtr"
  species_full <- "Medicago_truncatula"
  col_data <- metadata[c(19:24),]
  cat("The species is Medicago truncatula")
} else if (choice == 6) {
  species <- "Psa"
  species_full <- "Pisum_sativum"
  col_data <- metadata[c(7:12),]
  cat("The species is Pisum sativum")
} else if (choice == 7) {
  species <- "Pvu"
  species_full <- "Phaseolus_vulgaris"
  col_data <- metadata[c(31:36),]
  cat("The species is Phaseolus vulgaris")
#} else if (choice == 8) {
#  species <- "all"
#  species_full <- "all"
#  col_data <- metadata
#  cat ("You chose to check all species")
} else {
  cat("That species does not exist within this dataset")
}

# the counts file is the output of "make a file with all expression data from 
#all species" from StringTie
counts <- read.csv(file=paste("gene_count_matrix", species, "without_MSTRG.csv", 
                              sep="_"), 
                   header=TRUE, sep=",",row.names = 1)
# sum the reads with the same gene ID
counts <- counts[,2:8] #remove all the columns without reads
counts <- counts %>%
  group_by(Gene_id) %>% 
  summarise_all(.funs = sum,na.rm=T)
# make Gene_id the row names
genes <- counts$Gene_id
counts$Gene_id <- NULL #make the counts file only have read counts
counts <- as.data.frame(counts)
rownames(counts) <- genes


# Filter rowmean > 5
counts$rowmean <- rowMeans(counts)
counts_filtered <- subset(counts, rowmean > 5)
counts_filtered <- counts_filtered[,c(1:6)] #remove the rowmean column

# TMM Normalization
# make the DGEList:
TMM <- DGEList(counts_filtered)

#/ calculate TMM normalization factors:
TMM <- calcNormFactors(TMM)

#/ get the normalized counts:
counts_filtered_TMM <- cpm(TMM, log=FALSE)
#counts_filtered_TMM_log2 <- cpm(counts_filtered, log=T)
write.csv(as.data.frame(counts_filtered_TMM), file=paste("DESeq/Raw/deseq2", species, "TMM.csv", sep = "_"))
counts_filtered_TMM <- read.csv(paste("DESeq/Raw/deseq2", species, "TMM.csv", sep="_"), row.names = 1)

print(col_data$ID)
dds <- DESeqDataSetFromMatrix(countData=round(counts_filtered_TMM), 
                              colData = col_data,
                              design = ~ Tissue,
                              tidy = F)

DEdds <- DESeq(dds)
resultsNames(DEdds) 
model.matrix(~ Tissue, data = col_data) # 0 = Root Sample, 1 = Leaf Sample

# apeglm shrinkage, apply normalization to the data
resNorm_Tissue <- lfcShrink(DEdds, coef="Tissue_Leaf_vs_Root", type="apeglm") # tissue_Leaf_vs_Root
# Visualize the data as a data frame
resNorm_Tissue_df <- as.data.frame(resNorm_Tissue)
# Get the rownames into a column named ID
resNorm_Tissue_df$ID = rownames(resNorm_Tissue_df)
rownames(resNorm_Tissue_df) <- NULL
# Provide a summary: how many genes are there and how many are up or down regulated?
# Add info to the excel file "DESeq2_Summary.xlsx" in the "Analysis" folder
summary(resNorm_Tissue)
write.csv(as.data.frame(resNorm_Tissue_df), file=paste("DESeq/Sum/deseq2", species, "All_genes.csv", sep = "_"))
resNorm_Tissue_df <- read.csv(paste("DESeq/Sum/deseq2", species, "All_genes.csv", sep = "_"))

# Subset for adjusted p value < 0.05
res <- subset(resNorm_Tissue, padj < 0.05) 
res_df <- as.data.frame(res)
# Subset for log2FC higher than 1 (absolute value)
res_deseq <- res[abs(res$log2FoldChange) > 1,]   
res_deseq_df <- as.data.frame(res_deseq)
# Save results; substitute species name by using Ctrl + Shift + Alt + M
write.csv(as.data.frame(res_deseq), file=paste("DESeq/Sum/deseq2", species, "DEGs.csv", sep = "_"))
res_deseq <- read.csv(paste("DESeq/Sum/deseq2", species, "DEGs.csv", sep = "_"))

vsd <- vst(DEdds) #normalization
res_tibble <- (as_tibble(res_deseq, rownames = "gene"))
rownames(res_tibble) <- res_tibble$gene
# Check genes ordered by the log2FC
top_genes <- res_tibble %>% arrange(log2FoldChange)

# subset vsd matrix based on gene list
topVarGenes <- head(order(-rowVars(assay(vsd))),100)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

##### Volcano plot with all reads (not DEG only) #####
table <- resNorm_Tissue_df %>%
  mutate(diffexpressed = NULL) %>%
  mutate(
    diffexpressed = case_when(
      # if log2Foldchange > 1 and Pvalue < 0.05, set as "UP" 
      pvalue < 0.05 & log2FoldChange < -1 ~ "Root", 
      #  table$diffexpressed[table$logFC > 1 & table$Pvalue < 0.05] <- "UP"
      # if log2Foldchange < -1 and Pvalue < 0.05, set as "DOWN"
      pvalue < 0.05 & log2FoldChange > 1 ~ "Leaf", 
      #  table$diffexpressed[table$logFC < -1 & table$Pvalue < 0.05] <- "DOWN"
      TRUE ~ "NO"
    ))

mycolors <- c("olivedrab3", "brown", "black")
names(mycolors) <- c("Leaf", "Root", "NO")
png(filename = paste("DESeq/Output/Volcano", species_full,"TMM.png", sep="_"),
    width = 10,height = 6, units = "in", res = 600)
ggplot(data=table, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, size=baseMean)) + 
  geom_point() + 
  ggtitle(paste(species_full, "root vs leaf Volcano plot")) +
  theme_minimal() +
  scale_color_manual(values = mycolors) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off()

# breaks are the interval for the histogram
png(filename = paste("DESeq/Output/Histogram", species_full,"TMM_pvalue.png", sep="_"),
    width = 10,height = 6, units = "in", res = 600)
hist(table$pvalue[table$baseMean > 1], breaks = 20, col = "grey50", border = "white")
dev.off()

######## LogFC histogram #######
png(filename = paste("DESeq/Output/Histogram", species_full,"logFC_2_TMM.png", sep="_"),
   width = 10,height = 6, units = "in", res = 600) 
# Calculate histogram, but do not draw it
my_hist=hist(table$log2FoldChange , breaks=250  , plot=F)
# Color vector according to log2FC
my_color= ifelse(my_hist$breaks < -1, "brown" , ifelse (my_hist$breaks >=1, "olivedrab3", rgb(0.2,0.2,0.2,0.2) ))
# Final plot
plot(my_hist, col=my_color , border=F , main=paste("Histogram of", species_full, 
                                                   "log2 Fold Change Expression", 
                                                   sep=" "), 
     xlab="log(Fold Change)", xlim=c(-15,15) )
dev.off()

########################## Venn Diagram Leaf vs Root ######################
cpm_deseq <- read.csv(paste("DESeq/Raw/deseq2", species,"TMM.csv", sep = "_"), row.names = 1)

# Change to cpm files, 0.5 cut off <- <- <- <-
# Our library is 10M reads/sample; count 10 -> cpm 1
# Filter rowmean > 5
cpm_deseq$rowmean <- rowMeans(cpm_deseq)
cpm_deseq$ID <- rownames(cpm_deseq)
row.names(cpm_deseq) <- NULL

# Edit sample names according to Species
new_counts <- cpm_deseq
new_counts$averageLeaf <- rowMeans(cpm_deseq[, 1:3])
new_counts$averageRoot <- rowMeans(cpm_deseq[, 4:6])
new_counts$devLeaf <- apply(cpm_deseq[, 1:3], MARGIN =1, FUN = sd)
new_counts$devRoot <- apply(cpm_deseq[, 4:6], MARGIN =1, FUN = sd)
new_counts <- subset(new_counts, rowmean > 0.5)

# Get which genes are being expressed in leaves and roots (cpm_deseq >0.5)
new_counts <- new_counts %>%
  mutate(Leaf = case_when(
    averageLeaf>0.5 ~ "Yes", 
    averageLeaf<0.5 ~ "No",
    TRUE ~ "No"
    )) 

new_counts <- new_counts %>%
  mutate(Root = case_when(
    averageRoot>0.5 ~ "Yes", 
    averageRoot<0.5 ~ "No",
    TRUE ~ "No"
  )) 
genes_Leaf <- new_counts %>%
  filter(Leaf=="Yes") 
genes_Root <- new_counts %>%
  filter(Root=="Yes")
genes_Leaf <- genes_Leaf$ID
genes_Root <- genes_Root$ID

v <- venn.diagram(list(Leaf=genes_Leaf, Root=genes_Root),
                  fill = c("olivedrab3", "brown"),
                  alpha = c(0.5, 0.5),
                  filename=NULL)
png(filename = paste("DESeq/Output/Leaf_vs_Root", species, "Venn.png", sep="_"), width = 5,height = 5, units = "in", res = 600)
grid.newpage()
grid.draw(v)
dev.off()

################# Poisson Distance between samples and QC ######################
# Quality Control
library("pheatmap")

poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$name)
colnames(samplePoisDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
png(filename = paste("DESeq/Output/Leaf_vs_Root", species, "QC_Heatmap.png", sep="_"), 
    width = 5,height = 5, units = "in", res = 600)
pheatmap(samplePoisDistMatrix, clustering_distance_rows = poisd$dd, 
         clustering_distance_cols = poisd$dd, col = colors)
dev.off()

png(filename = paste("DESeq/Output/Leaf_vs_Root", species, "QC_PCA.png", sep="_"), 
    width = 1250,height = 750, units = "px")
PCA <- plotPCA(vsd, intgroup = #c(
  "Tissue"#, "ID")
  )
PCA + scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
  coord_fixed(ratio =2, xlim = NULL, ylim = NULL) + theme_minimal() + 
  geom_label(aes(label = name))
dev.off()

####### Histogram Daniel Luedke #####
my_color= ifelse(my_hist$breaks < -1, "brown" , ifelse (my_hist$breaks >=1, "olivedrab3", rgb(0.2,0.2,0.2,0.2) ))

#Fold change distribution. Shows a histogram with the frequency of the fold change distribution for a subset of genes e.g. all NLRs as definde in the subsetting step
png(filename = paste("DESeq/Output/Histogram", species,"logFC.png", sep="_"),
    width = 8,height = 10, units = "in", res = 600)
hist(resNorm_Tissue_df$log2FoldChange, breaks=250, #ylim=c(0,1500), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border=F,#"darkgrey", 
     col = my_color#, show.value=TRUE
     )
abline(v = c(-1,0,1), 
       col=c("brown","black","olivedrab3"), 
       lwd=1, lty=2)
dev.off()

############# single gene plot counts Daniel Luedke ##########
#generates a plot for the representation normalized expression for root and shoot tissue
plotCounts(dds, gene="Glyma.01G000100", intgroup = "Tissue", returnData = TRUE) %>% 
  ggplot() + aes(Tissue, count, color=Tissue) + 
  geom_jitter(width = 0.2, size = 4, shape = 4, stroke = 3) + ggtitle("Gene of Interest") +
  theme_bw()+scale_y_continuous(limits = c(0, 10), oob = scales::squish)

################ Total NLRs ########################
# Read file with all NLRs for each species. 
#The seqname does not have alternative splicing
TotalNLRs <- read.csv("Locus and Type Legume NLRs 06032024.csv", sep=",")

TSSpecies <- c("Lathyrus sativus", "Pisum sativum", "Lens culinaris", 
               "Cicer arietinum", "Medicago truncatula", "Phaseolus vulgaris",
               "Glycine max")

# filter to have only the NLRs from the species we want
TotalNLRs <- TotalNLRs %>%
  filter(Species %in% TSSpecies) 

#TotalNLRs <- TotalNLRs %>%
#  mutate(seqname = case_when(
#    Species == "Glycine max" ~ str_sub(TotalNLRs$seqname, end = -3), #remove the splicing
#    Species == "Phaseolus vulgaris" ~ str_sub(TotalNLRs$seqname, end = -3),
#    TRUE ~ seqname
#  ))

Lathyrus <- read.csv("DESeq/Sum/deseq2_Lsa_All_genes.csv", row.names = 1)
Pisum <- read.csv("DESeq/Sum/deseq2_Psa_All_genes.csv", row.names = 1)
Lens <- read.csv("DESeq/Sum/deseq2_Lcu_All_genes.csv", row.names = 1)
Cicer <- read.csv("DESeq/Sum/deseq2_Car_All_genes.csv", row.names = 1)
Medicago <- read.csv("DESeq/Sum/deseq2_Mtr_All_genes.csv", row.names = 1)
Phaseolus <- read.csv("DESeq/Sum/deseq2_Pvu_All_genes.csv", row.names = 1)
Glycine <- read.csv("DESeq/Sum/deseq2_Gma_All_genes.csv", row.names = 1)

#next: trim out the NLRs
Lathyrus <- Lathyrus %>% #revise ID
  filter(ID %in% TotalNLRs$seqname) %>%
  mutate(Species = "Lathyrus sativus")
Pisum <- Pisum %>%
  filter(ID %in% TotalNLRs$Locus) %>%
  mutate(Species = "Pisum sativum")
Lens <- Lens %>%
  filter(ID %in% TotalNLRs$seqname) %>%
  mutate(Species = "Lens culinaris")
Cicer <- Cicer %>%
  filter(ID %in% TotalNLRs$seqname) %>%
  mutate(Species = "Cicer arietinum")
Medicago <- Medicago %>% #revise ID
  filter(ID %in% TotalNLRs$Locus) %>%
  mutate(Species = "Medicago truncatula")
Phaseolus <- Phaseolus %>% #revise ID
  filter(ID %in% TotalNLRs$seqname) %>%
  mutate(Species = "Phaseolus vulgaris")
Glycine <- Glycine %>% #revise ID
  filter(ID %in% TotalNLRs$seqname) %>%
  mutate(Species = "Glycine max")


########################## NLRs logFC dotplot ###########################

NLRs <- rbind(Lathyrus, Pisum, Lens, Cicer, Medicago, Phaseolus, Glycine)
write.csv(NLRs, "DESeq/All_NLRs.csv")
NLRs <- read.csv("DESeq/All_NLRs.csv", row.names = 1)

NLRs <- NLRs %>%
  mutate(cond = if_else(log2FoldChange > 1, "olivedrab3", 
                        if_else(log2FoldChange < -1, "brown4", "black")))
#NLRs <- NLRs %>%
#  #drop_na(log2FoldChange) %>%
#  mutate(cond = case_when( #there's something messing up the tidyverse/dplyr packages... only load what I need
#    log2FoldChange>1 ~ 'olivedrab3',
#    log2FoldChange<-1 ~ 'brown4',
#    TRUE ~ 'black'   #anything that does not meet the criteria above
#  ))

ggplot(NLRs, aes(Species, log2FoldChange)) +
  geom_jitter(colour=NLRs$cond) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), 
             color = "black", linewidth=0.5)
ggsave("DESeq/Output/NLRs/Total_NLRs_logFC_per_Species_jitter3.png", width=1250, height=720, units='px', dpi = 125)

#Fold change distribution. Shows a histogram with the frequency of the fold change distribution for a subset of genes e.g. all NLRs as definde in the subsetting step
png(filename = "DESeq/Output/NLRs/Histogram_All_NLRs_logFC.png", 
    width = 1250,height = 720, units = "px")
hist(NLRs$log2FoldChange, breaks=250, #ylim=c(0,1500), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border="darkgrey",#,#"darkgrey", 
     col = "black"#, show.value=TRUE
)
abline(v = c(-1,0,1), 
       col=c("brown","black","olivedrab3"), 
       lwd=1, lty=2)
dev.off()

choice <- menu(c("Car", "Gma", "Lcu", "Lsa", "Mtr", "Psa", "Pvu"
), graphics = TRUE , title="Which species do you want to analyse?")

if (choice == 1) {species <- "Car"
  species_full <- "Cicer arietinum"
  cat("The species is Cicer arietinum")
} else if (choice == 2) {species <- "Gma"
  species_full <- "Glycine max"
  cat("The species is Glycine max")
} else if (choice == 3) {species <- "Lcu"
  species_full <- "Lens culinaris"
  cat("The species is Lens culinaris")
} else if (choice == 4) {species <- "Lsa"
  species_full <- "Lathyrus sativus"
  cat("The species is Lathyrus sativus")
} else if (choice == 5) {species <- "Mtr"
  species_full <- "Medicago truncatula"
  cat("The species is Medicago truncatula")
} else if (choice == 6) {species <- "Psa"
  species_full <- "Pisum sativum"
  cat("The species is Pisum sativum")
} else if (choice == 7) {species <- "Pvu"
  species_full <- "Phaseolus vulgaris"
  cat("The species is Phaseolus vulgaris")
} else {cat("That species does not exist within this dataset")}

NLRs_species <- NLRs[NLRs$Species==species_full,]
my_hist=hist(NLRs_species$log2FoldChange, breaks=100, plot = F, #ylim=c(0,1500), 
             #xlim=c(-15,15), 
             #lwd=2, 
             #lend=100#, border="darkgrey",#"darkgrey", 
             #col = "black"#, show.value=TRUE
)
# Color vector according to log2FC
my_color= ifelse(my_hist$breaks < -1, "brown" , ifelse (my_hist$breaks >=1, "olivedrab3", rgb(0.2,0.2,0.2,0.2) ))

#png(filename = paste("DESeq/Output/Histogram", species,"logFC.png", sep="_"),
#    width = 8,height = 10, units = "in", res = 600)
#hist(NLRs_species$log2FoldChange, breaks=250, #ylim=c(0,1500), 
#     xlim=c(-15,15), lwd=2, 
#     lend=100, border=F,#"darkgrey", 
#     col = my_color#, show.value=TRUE
#)
#abline(v = c(-1,0,1), 
#       col=c("brown","black","olivedrab3"), 
#       lwd=1, lty=2)
#dev.off()

#NLRs_species <- NLRs[NLRs$Species==species_full,]
#Fold change distribution. Shows a histogram with the frequency of the fold change distribution for a subset of genes e.g. all NLRs as definde in the subsetting step
png(filename = paste("DESeq/Output/NLRs/", species,"/Histogram_", species, 
                     "_NLRs_logFC_100_breaks_50_max.png", sep=""), 
    width = 1250,height = 720, units = "px")
hist(NLRs_species$log2FoldChange, breaks=100, ylim=c(0,50), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border="darkgrey",#"darkgrey", 
     col = my_color#, show.value=TRUE
)
abline(v = c(-1,0,1), 
       col=c("brown","black","olivedrab3"), 
       lwd=1, lty=2)
dev.off()


################ NLR Expression ##########################
deseq_All <- read.csv("DESeq/Raw/All_species_counts_without_MSTRG_not_filtered_or_normalized.csv", 
                      row.names =1)
# sum the reads with the same gene ID
deseq_All <- deseq_All[,2:9] #remove all the columns without reads
deseq_All <- deseq_All %>%
  group_by(Species, Gene_id) %>% 
  summarise_all(.funs = sum,na.rm=T) #552992 to 550640
# make Gene_id the row names
genes <- deseq_All$Gene_id
Species <- deseq_All$Species
deseq_All$Gene_id <- NULL #make the deseq_All file only have read counts
deseq_All$Species <- NULL
Species_Gene_id <- data.frame(Gene_id = genes, Species = Species) 
deseq_All <- as.data.frame(deseq_All)
rownames(deseq_All) <- genes

# Filter rowmean > 5
deseq_All$rowmean <- rowMeans(deseq_All)
deseq_All_filtered <- subset(deseq_All, rowmean > 5) # from 550640 to 211280
deseq_All_filtered <- deseq_All_filtered[,c(1:6)] #remove the rowmean column

# TMM Normalization
# make the DGEList:
TMM <- DGEList(deseq_All_filtered)
#/ calculate TMM normalization factors:
TMM <- calcNormFactors(TMM)
#/ get the normalized counts:
deseq_All_filtered_TMM <- cpm(TMM, log=FALSE)
deseq_All_filtered_TMM_df <- as.data.frame(deseq_All_filtered_TMM)
deseq_All_filtered_TMM_df$Gene_id <- row.names(deseq_All_filtered_TMM_df)
deseq_All_filtered_TMM_df <- left_join(deseq_All_filtered_TMM_df, Species_Gene_id, by="Gene_id")
#counts_filtered_TMM_log2 <- cpm(counts_filtered, log=T)
write.csv(as.data.frame(deseq_All_filtered_TMM_df), file="DESeq/Raw/All_species_counts_without_MSTRG_filted_and_normalized.csv")
deseq_All_filtered_TMM_df <- read.csv(file="DESeq/Raw/All_species_counts_without_MSTRG_filted_and_normalized.csv", row.names = 1)

NLRs_expression <- deseq_All_filtered_TMM_df %>%
  filter(Gene_id %in% NLRs$ID)

write.csv(NLRs_expression, file="DESeq/Raw/All_species_NLR_expression.csv") 

NLRs <- NLRs %>%
  mutate(cond = case_when(
    log2FoldChange>1 ~ 'olivedrab3',
    log2FoldChange<-1 ~ 'brown4',
    TRUE ~ 'black'   #anything that does not meet the criteria above
  ))

ggplot(NLRs, aes(Species, log2FoldChange)) +
  geom_jitter(colour=NLRs$cond) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), color = "black", linewidth=0.5)
#geom_point(data = ds, aes(y = mean), colour = 'red', size = 3)
ggsave("DESeq/Output/NLRs/Total_NLRs_logFC_jitter.png", width=1250, height=720, units='px', dpi = 125)

All_info_NLRs <- left_join(NLRs_expression, NLRs, by = c("Species", "Gene_id" = "ID"))
write.csv(All_info_NLRs, file = "DESeq/Raw/All_species_NLR_expression_logFC.csv")

##### NLR graphs per species #####

# In Geneious go to Workflow > Extract Sequences by Name > 
#Copy Sequence Names and Residues to excel

# choose the species
choice <- menu(c("Car", "Gma", "Lcu", "Lsa", "Mtr", "Psa", "Pvu"
                 #, "all"
), graphics = TRUE , title="Which species do you want to analyse?")

if (choice == 1) {
  species <- "Car"
  species_full <- "Cicer arietinum"
  cat("The species is Cicer arietinum")
} else if (choice == 2) {
  species <- "Gma"
  species_full <- "Glycine max"
  cat("The species is Glycine max")
} else if (choice == 3) {
  species <- "Lcu"
  species_full <- "Lens culinaris"
  cat("The species is Lens culinaris")
} else if (choice == 4) {
  species <- "Lsa"
  species_full <- "Lathyrus sativus"
  cat("The species is Lathyrus sativus")
} else if (choice == 5) {
  species <- "Mtr"
  species_full <- "Medicago truncatula"
  cat("The species is Medicago truncatula")
} else if (choice == 6) {
  species <- "Psa"
  species_full <- "Pisum sativum"
  cat("The species is Pisum sativum")
} else if (choice == 7) {
  species <- "Pvu"
  species_full <- "Phaseolus vulgaris"
  cat("The species is Phaseolus vulgaris")
  #} else if (choice == 8) {
  #  species <- "all"
  #  species_full <- "all"
  #  col_data <- metadata
  #  cat ("You chose to check all species")
} else {
  cat("That species does not exist within this dataset")
}

data <- NLRs_expression %>%
  filter(Species == species_full) %>%
  pivot_longer(
    cols = !c(Gene_id, Species),
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  mutate(Value=as.numeric(gsub(",", ".", Value)))

write.csv(data, paste("NLR_expression", species, "pivoted.csv", sep="_"))

# Heatmap 
# All_NLRs_Cicer_arietinum_CPM_DESeq_Heatmap.png
ggplot(data, aes(Sample, Gene_id, fill= Value)) + 
  geom_tile()  +
  scale_fill_gradient2(low="white", mid= "olivedrab3", high="brown4", 
                       midpoint = max(data$Value)/2) + 
  theme(axis.text.y = element_text(size = 5)) #+
#theme_bw()
ggsave(paste("DESeq/Output/NLRs/All_NLRs", species_full, "CPM_Heatmap.png", sep="_"), width=1250, 
       height=720, units='px', dpi = 125)

# All_NLRs_Cicer_arietinum_CPM_DESeq_Heatmap_bw.png
ggplot(data, aes(Sample, Gene_id, fill= Value)) + 
  geom_tile()  +
  scale_fill_gradient2(low="white", mid= "darkblue", high="black", 
                       midpoint = max(data$Value)/2) + 
  theme(axis.text.y = element_text(size = 5))
ggsave(paste("DESeq/Output/NLRs/All_NLRs", species_full, "CPM_Heatmap_bw.png", sep="_"), width=1250, 
       height=720, units='px', dpi = 125)

######### New counts ##############
# Edit sample names according to Species
new_counts <- NLRs_expression
new_counts$averageLeaf <- rowMeans(NLRs_expression[, 1:3])
new_counts$averageRoot <- rowMeans(NLRs_expression[, 4:6])
new_counts$devLeaf <- apply(NLRs_expression[, 1:3], MARGIN =1, FUN = sd)
new_counts$devRoot <- apply(NLRs_expression[, 4:6], MARGIN =1, FUN = sd)

new_counts <- new_counts %>%
  mutate(Leaf = case_when(
    averageLeaf>0.5 ~ "Yes", 
    averageLeaf<0.5 ~ "No",
    TRUE ~ "No"
  )) 
 
new_counts <- new_counts %>%
  mutate(Root = case_when(
    averageRoot>0.5 ~ "Yes", 
    averageRoot<0.5 ~ "No",
    TRUE ~ "No"
  )) 

write.csv(new_counts, file="DESeq/Raw/All_species_NLR_expression_Venn_table.csv") 

###### specific NLR Venn ######
new_counts <- read.csv("DESeq/Raw/All_species_NLR_expression_Venn_table.csv", row.names=1)

specific_counts <- new_counts %>%
  filter(Species == species_full)

genes_Leaf <- specific_counts %>%
  filter(Leaf=="Yes")
genes_Root <- specific_counts %>%
  filter(Root=="Yes")
genes_Leaf <- genes_Leaf$Gene_id
genes_Root <- genes_Root$Gene_id

v <- venn.diagram(list(Leaf=genes_Leaf, Root=genes_Root),
                  fill = c("olivedrab3", "brown"),
                  alpha = c(0.5, 0.5),
                  main = paste("Venn Diagram of", species_full,"NLRs expressed in Leaves and Roots"),
                  filename=NULL)
png(filename = paste("DESeq/Output/NLRs/", species, "Leaf_vs_Root_All_NLRs_Venn.png", 
                     sep="_"), width = 5,height = 5, units = "in", res = 600)
grid.newpage()
grid.draw(v)
dev.off()
############################# NLR Class Expression ###############
# Read file with all NLRs for each species. The seqname does not have alternative splicing
TotalNLRs <- read.csv("Locus and Type Legume NLRs 06032024.csv", sep=",")

TSSpecies <- c("Lathyrus sativus", "Pisum sativum", "Lens culinaris", 
               "Cicer arietinum", "Medicago truncatula", "Phaseolus vulgaris",
               "Glycine max")

# filter to have only the NLRs from the species we want
TotalNLRs <- TotalNLRs %>%
  filter(Species %in% TSSpecies) 

new_counts <- read.csv(file="DESeq/Raw/All_species_NLR_expression_Venn_table.csv", row.names=1) 

TotalClassNLRs_1 <- inner_join(All_info_NLRs, TotalNLRs, by= c("Gene_id"="seqname", "Species"))
TotalClassNLRs_2 <- inner_join(All_info_NLRs, TotalNLRs, by = c("Gene_id"="Locus", "Species"))
TotalClassNLRs <- full_join(TotalClassNLRs_1,TotalClassNLRs_2)

TotalClassNLRs <- TotalClassNLRs %>%
  mutate(Locus = case_when(
    is.na(Locus) ~ coded_by,
    TRUE ~ Locus
  )) %>%
  group_by(Locus) %>%
  arrange(desc(Locus)) %>%
  slice(1) #Slice out alternative splicing
TotalClassNLRs <- TotalClassNLRs %>%
  mutate(seqname=case_when(
    is.na(seqname) ~ Gene_id,
    TRUE ~ seqname
  ))

TotalClassNLRs <- left_join(TotalClassNLRs, new_counts)

write.csv(TotalClassNLRs, "DESeq/Output/NLRs/All_NLRs_Expression_logFC_Locus_Class.csv")

##### Specific Pie chart #####
LeafNLRs <- TotalClassNLRs %>%
  filter(Species == species_full) %>%
  group_by(Class) %>%
  count(Leaf)
OnlyLeafNLRs <- LeafNLRs %>%
  filter(Leaf=="Yes")
RootNLRs <- TotalClassNLRs %>%
  filter(Species == species_full) %>%
  group_by(Class) %>%
  count(Root)
OnlyRootNLRs <- RootNLRs %>%
  filter(Root=="Yes")

# For leaves
# Compute the position of labels
data <- OnlyLeafNLRs %>% 
  arrange(order(Class)) %>%
  mutate(Sum=sum(OnlyLeafNLRs$n)) %>%
  mutate(prop = (n / Sum) *100) %>%
  mutate(ypos = cumsum(prop)+ 1.1*prop )

colors = c("turquoise","lightgoldenrod2", "thistle3",  "salmon", "dodgerblue")
png(filename = paste("DESeq/Output/NLRs/Leaf_All_NLRs", species, "Class_Pie.png", 
                     sep="_"), width = 5,height = 5, units = "in", res = 600)
pie(data$n, labels = paste(sep=", ", data$Class, data$n), 
    main=paste(species_full, "NLR Expression per Class in Leaves", sep = " "),
    col= colors)
dev.off()

# For Roots
# Compute the position of labels
data <- OnlyRootNLRs %>% 
  arrange(order(Class)) %>%
  mutate(Sum=sum(OnlyRootNLRs$n)) %>%
  mutate(prop = (n / Sum) *100) %>%
  mutate(ypos = cumsum(prop)+1*prop )

colors = c("turquoise","lightgoldenrod2", "thistle3", 
           #"green", #Phaseolus
           "salmon", "dodgerblue")
png(filename = paste("DESeq/Output/NLRs/Root_All_NLRs", species, "Class_Pie.png", 
                     sep="_"), width = 5,height = 5, units = "in", res = 600)
pie(data$n, labels = paste(sep=", ", data$Class, data$n), 
    main=paste(species_full, "NLR Expression per Class in Roots", sep = " "),
    col= colors)
dev.off()

############################ Venn Diagram NLRs ##########################
TotalClassNLRs <- read.csv("DESeq/Output/NLRs/All_NLRs_Expression_logFC_Locus_Class.csv", row.names=1)

genes_Leaf <- TotalClassNLRs %>%
  filter(Leaf=="Yes")
genes_Root <- TotalClassNLRs %>%
  filter(Root=="Yes")
genes_Leaf <- genes_Leaf$Gene_id
genes_Root <- genes_Root$Gene_id

v <- venn.diagram(list(Leaf=genes_Leaf, Root=genes_Root),
                  fill = c("olivedrab3", "brown"),
                  alpha = c(0.5, 0.5),
                  filename=NULL)
png(filename = "DESeq/Output/NLRs/Leaf_vs_Root_All_NLRs_Venn.png", width = 5,height = 5, units = "in", res = 600)
grid.newpage()
grid.draw(v)
dev.off()

############### Expression of Conserved NLRs ###############
Conserved_nodes <- read_xlsx("NB_ARC_NLR_nodes_filtered.xlsx")

# Separating the aminoacid sequence from the IDs
Conserved_nodes <- Conserved_nodes %>%
  separate(NLR_sequence, c("seqname","sequence"), sep=": ")

# Counting how many aminoacids the protein has into a new column called 'Length'
Conserved_nodes$sequence <- str_replace_all(Conserved_nodes$sequence," ","") # remove extra spaces
Conserved_nodes$sequence <- str_replace_all(Conserved_nodes$sequence,"-","") # remove gaps -
Conserved_nodes <- Conserved_nodes %>%
  mutate("Length" = nchar(Conserved_nodes$sequence))

# make seqname not have transcript info
Conserved_nodes <- Conserved_nodes %>%
  mutate(seqname = case_when(
    Species %in% c("Medicago truncatula", "Pisum sativum") 
    ~ seqname,
    # remove everything after the last .
    TRUE ~ gsub("[.][^.]+$", "", seqname)
  )) %>%
  mutate(seqname=case_when(
    # find _() in the seqname. These will not have a match in Gene_id,
    # we need to correct for that before joining the conserved NLRs with 
    # expression information
    grepl("_(", seqname, fixed = TRUE) ~ str_sub(Conserved_nodes$seqname, end = -5),
    TRUE ~ seqname
  ))

TotalClassNLRs <- TotalClassNLRs %>%
  mutate(Gene_id = case_when(
    Species == "Pisum sativum" ~ seqname,
    TRUE ~ Gene_id
  ))
Conserved_nodes_expression <- full_join(Conserved_nodes, TotalClassNLRs, 
                                        by = c("seqname"="Gene_id"))
# Tidying the table
Conserved_nodes_expression <- Conserved_nodes_expression %>%
  mutate(Leaf.1=case_when(
    is.na(Leaf.1) ~ 0,
    TRUE ~ Leaf.1)) %>% 
  mutate(Leaf.2=case_when(is.na(Leaf.2) ~ 0,TRUE ~ Leaf.2)) %>%
  mutate(Leaf.3=case_when(is.na(Leaf.3) ~ 0,TRUE ~ Leaf.3)) %>%
  mutate(Root.1=case_when(is.na(Root.1) ~ 0,TRUE ~ Root.1)) %>%
  mutate(Root.2=case_when(is.na(Root.2) ~ 0,TRUE ~ Root.2)) %>%
  mutate(Root.3=case_when(is.na(Root.3) ~ 0,TRUE ~ Root.3)) %>%
  mutate(Leaf=case_when(is.na(Leaf) ~ "No",TRUE ~ Leaf)) %>%
  mutate(Root=case_when(is.na(Root) ~ "No",TRUE ~ Root)) %>%
  filter(!is.na(NLR_name)) %>%
  mutate(Species.x = case_when(
    startsWith(as.character(NLR_name), "XP_") ~ "Pisum sativum",
    startsWith(as.character(NLR_name), "Mtrun") ~ "Medicago truncatula",
    startsWith(as.character(NLR_name), "phavu") ~ "Phaseolus vulgaris",
    startsWith(as.character(NLR_name), "Lcu") ~ "Lens culinaris",
    startsWith(as.character(NLR_name), "Glyma") ~ "Glycine max",
    startsWith(as.character(NLR_name), "g") ~ "Lathyrus sativus",
    startsWith(as.character(NLR_name), "cicar") ~ "Cicer arietinum",
    TRUE ~ Species.x
  ))
Conserved_nodes_expression <- subset(Conserved_nodes_expression, select=
                                       -c(Class.y, Species.y, seqname.y.y, seqname.y,
                                          Type, RefSeq, Name, gene, X1, ...13, 
                                          ...14))
columns <- c("Species", "Class", "Node", "NLR_variant", "seqname", "sequence",
             "Comment", "variant_length", "Leaf_1", "Leaf_2", "Leaf_3",
             "Root_1","Root_2","Root_3","Mean_expression","log2FC","lfcSE",
             "p-value","FDR","color","original","Organism","Genome",
             "original_sequence","original_length","Locus","Coded_by",
             "Leaf_mean_expression", "Root_mean_expression", "Leaf_stdev_expression",
             "Root_stdev_expression", "Leaf_expressed", "Root_expressed")
colnames(Conserved_nodes_expression) <- columns

# Create a column with Tissue assignment (NE, L, R, Both)
Conserved_nodes_expression <- Conserved_nodes_expression %>%
  mutate(Tissue = case_when(
    Leaf_expressed == "Yes" & Root_expressed == "No" ~ "Leaf",
    Leaf_expressed == "No" & Root_expressed == "Yes" ~ "Root",
    Leaf_expressed == "Yes" & Root_expressed == "Yes" ~ "Both",
    Leaf_expressed == "No" & Root_expressed == "No" ~ "Not Expressed",
    TRUE ~ "Invalid"
  ))

write.csv(Conserved_nodes_expression, "Conserved_nodes_expression.csv")
### Bar chart with presence/absence of NLRs with Class info ####
# Count how many NLRs are in each class, for each Species
TotalClassNLRs <- TotalClassNLRs %>%
  mutate(Tissue = case_when(
    Leaf == "Yes" & Root == "No" ~ "Leaf",
    Leaf == "No" & Root == "Yes" ~ "Root",
    Leaf == "Yes" & Root == "Yes" ~ "Both",
    Leaf == "No" & Root == "No" ~ "Not Expressed",
    TRUE ~ "Invalid"
  )) %>%
  mutate(logFC = case_when(
    log2FoldChange < -1 ~ "More expressed in Roots",
    log2FoldChange > 1 ~ "More expressed in Leaves",
    TRUE ~ "Same expression"
  )) %>%
  mutate(Tissue_logFC = case_when(
    Tissue == "Leaf" ~ "Leaf-specific",
    Tissue == "Root" ~ "Root-specific",
    Tissue == "Not Expressed" ~ "Not Expressed",
    Tissue == "Both" & logFC == "More expressed in Leaves" ~ "Leaf tendency",
    Tissue == "Both" & logFC == "More expressed in Roots" ~ "Root tendency",
    Tissue == "Both" & logFC == "Same expression" ~ "Same expression",
    TRUE ~ "Invalid"
  ))

write.csv(TotalClassNLRs, "Total_Class_NLRs_logFC.csv")
TotalClassNLRs <- read.csv("Total_Class_NLRs_logFC.csv", row.names = 1)

TotalClassNLRs <- TotalClassNLRs %>%
  mutate(cond_logFC = case_when(
    Tissue_logFC == "Not Expressed" ~ "#000000",
    Tissue_logFC == "Same expression" ~ "blue",
    Tissue_logFC == "Root-specific" ~ "#834c3b",
    Tissue_logFC == "Leaf-specific" ~ "#488200",
    Tissue_logFC == "Leaf tendency" ~ "#b6c630",
    Tissue_logFC == "Root tendency" ~ "#FFAA00",
    TRUE ~ "red"
  ))

TotalClassNLRs$Species = factor(TotalClassNLRs$Species,
                                levels =c("Lathyrus sativus", "Pisum sativum",
                                          "Lens culinaris", "Medicago truncatula",
                                          "Cicer arietinum","Glycine max", 
                                          "Phaseolus vulgaris"), ordered =TRUE) 


ggplot(TotalClassNLRs, aes(Species, log2FoldChange)) +
  geom_jitter(colour=TotalClassNLRs$cond_logFC, size = 1.5) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), 
             color = "black", linewidth=0.5)
ggsave("DESeq/Output/NLRs/Total_Class_NLRs_logFC_per_Species_jitter5.png", width=1250, height=720, units='px', dpi = 125)


TotalClassNLRs_Number <- TotalClassNLRs %>%
  group_by(Species, Tissue) %>%
  count(Class)
# OR
TotalClassNLRs_Number <- TotalClassNLRs %>%
  group_by(Species, Tissue_logFC) %>%
  count(Class)

TotalClassNLRs_Number <- TotalClassNLRs_Number %>%
  group_by(Species) %>%
  mutate(ClassPercentage = n/(sum(n))*100) %>%#Percentage of NLRs with that class
  mutate(Total = sum(n))

TotalClassNLRs_Number$Species = factor(TotalClassNLRs_Number$Species, 
                                               levels =c("Lathyrus sativus", "Pisum sativum", 
                                                         "Lens culinaris", "Medicago truncatula",
                                                         "Cicer arietinum",
                                                         "Glycine max", "Phaseolus vulgaris"), ordered =TRUE) 

TotalClassNLRs_Number$Class = factor(TotalClassNLRs_Number$Class, 
                                           levels =c("CC-NLR", "TIR-NLR", 
                                                     "CCR-NLR", 
                                                     "CCG10-NLR",
                                                     "TNP", "Other"), ordered =TRUE) 
TotalClassNLRs_Number$Tissue = factor(TotalClassNLRs_Number$Tissue, 
                                            levels =c("Not Expressed", "Both", "Leaf", "Root"), ordered =TRUE) 
TotalClassNLRs_Number$Tissue_logFC = factor(TotalClassNLRs_Number$Tissue_logFC, 
                                      levels =c("Not Expressed", "Leaf-specific",
                                                "Leaf tendency", "Same expression", 
                                                 "Root tendency", "Root-specific"), ordered =TRUE) 

ggplot() +
  geom_bar(data=TotalClassNLRs_Number, aes(y = n, x = Tissue, fill = Class)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species) + theme_bw() + 
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_number_tendency.png", width = 1250, height = 720, units="px", dpi = 120)


ggplot() +
  geom_bar(data=TotalClassNLRs_Number, aes(y = ClassPercentage, x = Tissue, fill = Class)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species) + theme_bw() + scale_y_continuous(minor_breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_percentage_tendency.png", width = 1250, height = 720, units="px", dpi = 120)

TotalClassNLRs_Number <- TotalClassNLRs_Number %>%
  group_by(Species, Tissue_logFC) %>%
  mutate(Percentage_Tissue = sum(ClassPercentage)) %>%
  mutate(Number_Tissue_logFC = sum(n))

write.csv(TotalClassNLRs_Number, "Total_Class_NLRs_logFC_tendency_Summary.csv")
TotalClassNLRs_Number <- read.csv("Total_Class_NLRs_logFC_tendency_Summary.csv")


my.lines = data.frame(x=1, y=seq(from=7, to=798, by = 7))
Conserved_NLRs_pivot[,c(12:13)] <- NULL

ggplot(TotalClassNLRs_Number,                                # Draw heatmap-like plot
       aes(Sample, n, fill = CPM)) +  geom_tile() + 
  scale_y_discrete(breaks = TotalClassNLRs_Number$n 
                   # , labels = Conserved_NLRs$seqname
  ) + #geom_hline(yintercept = 0.5 + 0:25908, colour = "black", linewidth = 0.5) +
  
  # scale_z_log10() +
  scale_fill_gradient2(low="white", mid = "#06F9C6",high="#7F0318", trans = "log" , na.value="white") + 
  ggtitle("Expression of conserved NLRs", subtitle = "Leaves and Roots") +
  ylab("Medicago truncatula homologue") + #geom_segment(data=my.lines, aes(x=1, y, xend=6.5, yend=y), linewidth=0.01, inherit.aes=F)
  theme_minimal() +
  # this line allows to divide the graph in NLR homologue groups
  geom_hline(yintercept = seq(from=0.5, to=798, by = 7), linewidth=0.00001)

ggplot() +
  geom_bar(data=TotalClassNLRs_Number, aes(y = n, x = Tissue, fill = Class)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species) + theme_bw() + 
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
# The NA is probably a truncated NLR, not considered for the NLR study

ggplot() +
  geom_bar(data=Conserved_NLRs_bar, aes(y = ClassPercentage, x = Tissue, fill = Class)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species.x) + theme_bw() + scale_y_continuous(minor_breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))

ggplot() +
  geom_bar(data=Conserved_NLRs_bar, aes(y = n, x = Tissue, fill = Class)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species.x) + theme_bw() + 
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
# The NA is probably a truncated NLR, not considered for the NLR study

ggplot() +
  geom_bar(data=Conserved_NLRs_bar, aes(y = ClassPercentage, x = Tissue, fill = Class)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species.x) + theme_bw() + scale_y_continuous(minor_breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))


##### Create line expression plot where which line is a species #####
# Overlaying two plots
# Review later
Conserved_nodes_expression$Clade <- paste(Conserved_nodes_expression$Class, 
                                         Conserved_nodes_expression$Node, sep="_")
Conserved_nodes_expression <- Conserved_nodes_expression[-which(Conserved_nodes_expression$Clade == "CCG10-NLR_12"), ]
write.csv(Conserved_nodes_expression, "Conserved_nodes_expression_no_problematic_CCG10.csv")
Conserved_nodes_expression <- read.csv("Conserved_nodes_expression_no_problematic_CCG10.csv", row.names=1)

Conserved_nodes_expression$Clade2 <- gsub("-", "_", Conserved_nodes_expression$Clade)
Conserved_nodes_expression_2 <- Conserved_nodes_expression %>%
  drop_na(log2FC) %>%
  mutate(Species_Color = case_when(
    Species == "Phaseolus vulgaris" ~ "#FF1100",
    Species == "Glycine max" ~ "#ff7400",
    Species == "Cicer arietinum" ~ "#FFBB00",
    Species == "Medicago truncatula" ~ "#488200",
    Species == "Lens culinaris" ~ "#00b67e",
    Species == "Pisum sativum" ~ "#003b84",
    Species == "Lathyrus sativus" ~ "#3f0051",
    TRUE ~ "black"
  ))

library("multcompView")

anova_1 <- aov(log2FC ~ Clade2, data = Conserved_nodes_expression_2)
summary(anova_1)
tukey_1 <- TukeyHSD(anova_1)
print(tukey_1)
cld_1 <- multcompLetters4(anova_1, tukey_1) # compact letter display
print(cld_1)
# table with factors and 3rd quantile
Tk_1 <- group_by(Conserved_nodes_expression_2, Clade2) %>%
  summarise(mean=mean(log2FC), quant = quantile(log2FC, probs = 0.75)) %>%
  arrange(desc(mean))
# extracting the compact letter display and adding to the Tk table
cld_1 <- as.data.frame.list(cld_1$Clade2)
Tk_1$cld <- cld_1$Letters
print(Tk_1)
#Conserved_nodes_expression_2$Clade2 <- factor(Conserved_nodes_expression_2$Clade2,
#                                   levels=c("CC", "TIR",
#                                            "CCG10", "CCR", "TNP", "Other"))
library("scales")
ggplot(Conserved_nodes_expression_2, aes(x=Clade2, y=log2FC#, fill=Clade2
                                         )) + 
  geom_boxplot() +
  labs(x="NLR Clade", y="log2FC") +
  theme_bw() + 
  #scale_y_continuous(labels = scales::percent) +
  geom_jitter(color=Conserved_nodes_expression_2$Species_Color, size=1, alpha=0.9,
              width = 0.2, height = 0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  labs(title = "logFC of the expression of Conserved NLRs", subtitle = "In leaf and root tissue") +
  geom_text(data = Tk_1, aes(x = Clade2, y = quant, label = cld), size = 3, vjust=-1, hjust =-1) #+
  #scale_fill_brewer(c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
  #scale_fill_manual(values=c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
# saving the final figure
ggsave("DESeq/Output/NLRs/Conserved_NLR_nodes_logFC_jitterplot_2.png", width = 1250, height = 720, units="px", dpi = 120)


ggplot(subset(Conserved_nodes_expression, Species=="Phaseolus vulgaris"), 
       aes(Clade, log2FC)) + 
  geom_jitter(color = "#FF1100", size=2, width = 0.2, height = 0) + 
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Glycine max"), 
            color = "#ff7400", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Cicer arietinum"), 
            color = "#FFBB00", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Medicago truncatula"), 
            color = "#488200", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Lens culinaris"), 
            color = "#00b67e", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Pisum sativum"), 
            color = "#003b84", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Lathyrus sativus"), 
            color = "#3f0051", size=2, width = 0.2, height = 0) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  labs(title = "logFC of the expression of Conserved NLRs", subtitle = "In leaf and root tissue")
  
ggsave("DESeq/Output/NLRs/Conserved_NLR_nodes_logFC_jitterplot.png", width=1250, 
       height=720, units='px', dpi = 100)

ggplot(subset(Conserved_nodes_expression, Species=="Phaseolus vulgaris"), 
       aes(Clade, Leaf_mean_expression)) + 
  geom_jitter(color = "#FF1100", size=2, width = 0.2, height = 0) + 
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Glycine max"), 
            color = "#ff7400", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Cicer arietinum"), 
            color = "#FFBB00", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Medicago truncatula"), 
            color = "#488200", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Lens culinaris"), 
            color = "#00b67e", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Pisum sativum"), 
            color = "#003b84", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Lathyrus sativus"), 
            color = "#3f0051", size=2, width = 0.2, height = 0) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  labs(title = "Expression of Conserved NLRs", subtitle = "In leaf tissue") +
  ylim(c(0,55))


ggsave("DESeq/Output/NLRs/Conserved_NLR_nodes_TMM_Leaf_jitterplot.png", width=1250, 
       height=720, units='px', dpi = 100)

ggplot(subset(Conserved_nodes_expression, Species=="Phaseolus vulgaris"), 
       aes(Clade, Root_mean_expression)) + 
  geom_jitter(color = "#FF1100", size=2, width = 0.2, height = 0) + 
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Glycine max"), 
            color = "#ff7400", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Cicer arietinum"), 
            color = "#FFBB00", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Medicago truncatula"), 
            color = "#488200", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Lens culinaris"), 
            color = "#00b67e", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Pisum sativum"), 
            color = "#003b84", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(Conserved_nodes_expression, Species=="Lathyrus sativus"), 
            color = "#3f0051", size=2, width = 0.2, height = 0) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  labs(title = "Expression of Conserved NLRs", subtitle = "In root tissue") +
  ylim(c(0,55))

ggsave("DESeq/Output/NLRs/Conserved_NLR_nodes_TMM_Root_jitterplot.png", width=1250, 
       height=720, units='px', dpi = 100)
