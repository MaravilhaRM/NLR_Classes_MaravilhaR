# Author: Rita Maravilha Marques
# Script to perform expression analysis and generate figures

# Install required packages if needed
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

# Load required libraries
library(DESeq2) # for differential expression analysis
library(edgeR) # for normalizing RNA Seq data
library(apeglm) # for normalizing logFC values
#library(IsoformSwitchAnalyzeR) #would enable statistical identification of isoform switches
library(dplyr) # for data wrangling
library(ggplot2) # visualization
library(stringr) # for common string operations in R
library(tidyverse) # for data wrangling
library(readxl) # read Excel files
library(VennDiagram) # generate Venn Diagrams
library(pheatmap) # generate heatmaps
library("RColorBrewer") # Colors
library("PoiClaClu") # for performing classification and clustering of RNA Seq samples based on Poisson model

# Set relative working directory according to the script's location
setwd("../Tissue_Spec/Counts_2024/112024")

##### Creating BUSCO id counts #####
setwd("../Tissue_Spec/BUSCO PCA")

# before reading in the tabular BUSCO data, remove the # from the third line
# and replace Busco id by Busco_id (the space isn't a good character)
BUSCO_Car <- read.table("BUSCO_Car_table.tabular", 
                        sep="\t", # separator is a tab
                        header=TRUE, # has a header
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
#the max score. The total number of BUSCOs for the fabales dataset should be 5366
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

# join all BUSCO ids on a single table
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

# write the results
write.csv(BUSCO_all, "BUSCO_ids_Tissue_Spec_112024.csv")
setwd("../Tissue_Spec/Counts_2024/112024")

# read if needed
BUSCO_all <- read.csv("BUSCO_ids_Tissue_Spec_112024.csv", row.names = 1)

##### make a file with all expression data from all species #####
# First, let's keep only the transcripts found with Stringtie
# Car - Cicer arietinum
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


# Gma - Glycine max
counts_Gma <- read.csv(file="gene_count_matrix_Gma.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Gma$Gene_id <- sapply(strsplit(counts_Gma$gene_id, "\\|"), function(x) tail(x, 1))
# Replace gene-GLYMA_ by Glyma.
counts_Gma <- counts_Gma %>%
  mutate(Gene_id = str_replace_all(Gene_id, "gene-GLYMA_", "Glyma.")) 
write.csv(counts_Gma, "gene_count_matrix_Gma_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
counts_Gma <- counts_Gma[grep("^Glyma.", counts_Gma$Gene_id),]
# remove the v4 at the end
counts_Gma <- counts_Gma %>%
  mutate(Gene_id = str_sub(Gene_id, end = -3)) %>% #remove the last characters to match the seqname
  mutate(Species = "Glycine max")
write.csv(counts_Gma, "gene_count_matrix_Gma_without_MSTRG.csv")


# Lcu - Lens culinaris
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


# Lsa - Lathyrus sativus
counts_Lsa <- read.csv(file="gene_count_matrix_Lsa.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Lsa$Gene_id <- sapply(strsplit(counts_Lsa$gene_id, "\\|"), function(x) tail(x, 1))
# Replace gene- by nothing
counts_Lsa <- counts_Lsa %>%
  mutate(Gene_id = str_replace_all(Gene_id, "gene-", "")) 
write.csv(counts_Lsa, "gene_count_matrix_Lsa_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
#counts_Lsa <- counts_Lsa[grep("^g", counts_Lsa$Gene_id),]
counts_Lsa <- counts_Lsa[grep("^LATHSAT", counts_Lsa$Gene_id),]
counts_Lsa <- counts_Lsa %>%
  mutate(Species = "Lathyrus sativus")
write.csv(counts_Lsa, "gene_count_matrix_Lsa_without_MSTRG.csv")


# Mtr - Medicago truncatula
counts_Mtr <- read.csv(file="gene_count_matrix_Mtr.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Mtr$Gene_id <- sapply(strsplit(counts_Mtr$gene_id, "\\|"), function(x) tail(x, 1))
# make the sequence name match the protein ID
counts_Mtr <- counts_Mtr %>%
  #mutate(Gene_id = str_replace_all(Gene_id, "MtrunA17Chr0c0", "MtrunA17_Chr")) %>%
  #mutate(Gene_id = str_replace_all(Gene_id, "R", "g")) %>% #Rs are repeat regions
  mutate(Gene_id = str_replace_all(Gene_id, "MtrunA17Chr", "MtrunA17_Chr"))
write.csv(counts_Mtr, "gene_count_matrix_Mtr_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
# this annotation has some genes beginning with Mt, others with MtrunA17
counts_Mtr <- counts_Mtr[grep("^Mt", counts_Mtr$Gene_id),]
counts_Mtr <- counts_Mtr %>%
  mutate(Species = "Medicago truncatula")
write.csv(counts_Mtr, "gene_count_matrix_Mtr_without_MSTRG.csv")


# Psa - Pisum sativum
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


# Pvu - Phaseolus vulgaris
counts_Pvu <- read.csv(file="gene_count_matrix_Pvu.csv", header=TRUE, 
                       sep=",")
# I need to separate the last part of the string by |
# Extracting the substring after "|"
counts_Pvu$Gene_id <- sapply(strsplit(counts_Pvu$gene_id, "\\|"), function(x) tail(x, 1))
# Replace gene-PHAVU_ by phavu.G19833.gnm2.ann1.Phvul.
counts_Pvu <- counts_Pvu %>%
  mutate(Gene_id = str_replace_all(Gene_id, "gene-PHAVU_", "phavu.G19833.gnm2.ann1.Phvul."))
write.csv(counts_Pvu, "gene_count_matrix_Pvu_with_MSTRG.csv")
# Filter out elements that don't contain the desired substring
counts_Pvu <- counts_Pvu[grep("^phavu.G19833.gnm2.ann1.Phvul.", counts_Pvu$Gene_id),]
#  remove the g at the end
counts_Pvu <- counts_Pvu %>%
  mutate(Gene_id = str_sub(Gene_id, end = -2)) %>% #remove the last characters to match the seqname
  mutate(Species = "Phaseolus vulgaris")
write.csv(counts_Pvu, "gene_count_matrix_Pvu_without_MSTRG.csv")

##### Example on how to add New Lsa genome JIC 2.1.1 gff #####
# Here is an example on how to add a new gff file to the previous Locus file

# Load new file
gff_Lsa <- read_xlsx("JIC_Lsat_v2.1.1_genomic_edited.gff.xlsx", skip=6)
# select only the rows pertaining to gene or mRNA
gff_Lsa <- gff_Lsa %>%
  filter(...3 %in% "mRNA")
gff_Lsa <- gff_Lsa[, 1:13]
colnames(gff_Lsa) <- c("seqname", "source", "feature", "start", "end", "score", 
                       "strand", "frame", "Locus", "Locus2", 
                       "Gene_id", "gene_biotype", "Locus_tag")
# tidy
gff_Lsa <- gff_Lsa %>%
  mutate(Locus = str_replace_all(Locus, "ID=rna-", "")) %>%
  mutate(Locus2 = str_replace_all(Locus2, "Parent=gene-", "")) %>%
  mutate(Gene_id = str_replace_all(Gene_id, "Note=ID:", "")) %>%
  mutate(Gene_id = str_replace_all(Gene_id, "%3B~source:AUGUSTUS", "")) %>%
  mutate(Gene_id = str_replace_all(Gene_id, "%3B~source:GeneMark.hmm3", "")) %>%
  mutate(Gene_id = str_replace_all(Gene_id, "%3B~source:gmst", "")) %>%
  mutate(gene_biotype = str_replace_all(gene_biotype, "gbkey=", "")) %>%
  mutate(Locus_tag = str_replace_all(Locus_tag, "locus_tag=", ""))
write.csv(gff_Lsa, "gff_Lsa.csv")
gff_Lsa <- read.csv("gff_Lsa.csv", row.names=1)

# Load previous Locus file
Locus <- read_csv("Locus and Type Legume NLRs 06032024.csv")

# Replace the previous genome with the newest
Locus <- Locus %>%
  mutate(Genome= str_replace_all(Genome, "gp_hifiasm_assembly.v1.0", "JIC_Lsat_v2.1.1"))

# Join the previous Locus file with the new info
new_Locus <- left_join(Locus, gff_Lsa, by=c("original"="Gene_id"))
new_Locus <- new_Locus %>%
  mutate(Locus.x = case_when(
    Genome == "JIC_Lsat_v2.1.1" ~ Locus2,
    TRUE ~ Locus.x
  )) %>%
  mutate(coded_by = case_when(
    Genome == "JIC_Lsat_v2.1.1" ~ Locus.y,
    TRUE ~ coded_by
  ))
colnames(new_Locus) <- c("Species", "seqname", "Class", "original", "Organism",
                         "Genome", "Type", "sequence", "Length", "Locus", 
                         "coded_by", "RefSeq", "Name", "gene", "X1", "16", "17")
new_Locus <- new_Locus[, 1:17]

# write it up again
write.csv(new_Locus, "Locus and Type Legume NLRs 28102024.csv")
new_Locus <- read.csv("Locus and Type Legume NLRs 28102024.csv", row.names=1)

##### Make a counts file with all species ######
# Start with Cicer arietinum counts and add the others
counts_all <- counts_Car
# change column names to Leaf and Root instead of sample names
colnames(counts_all)[2:7] <- c("Leaf 1", "Leaf 2", "Leaf 3", 
                               "Root 1", "Root 2", "Root 3") 

# join all counts into counts_all
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
# 28269 + 52814 + 55844 + 31721 + 297567 + 62857 + 28131 = 557203

# write the counts file with all species
write.csv(counts_all, "All_species_counts_without_MSTRG_not_filtered_or_normalized_112024.csv")

##### Merge counts with loci info ######
# Read in the BUSCO, counts and Locus files
BUSCO_all <- read.csv("BUSCO_ids_Tissue_Spec_112024.csv", row.names = 1)
counts_all <- read.csv("All_species_counts_without_MSTRG_not_filtered_or_normalized_112024.csv", row.names = 1)
LOCUS <- read_tsv("Fabales_LOCUS_28102024_plus_JIC_Lsa.tsv")

# Rewrite genome to accomodate new annotations
LOCUS_tissue_spec <- LOCUS %>%
  mutate(Genome = case_when(
    Organism == "Pisum sativum ZW6" ~ "CAAS_Psat_ZW6_1.0",
    TRUE ~ Genome
  )) %>%
  filter(Genome %in% c("cicar.CDCFrontier.gnm1.ann1.nRhs.protein", 
                       "Glycine max Williams 82", 
                       "Lens culinaris Redberry v2.0",
                         "JIC_Lsat_v2.1.1", "MtrunA17r5.0-ANR-EGN-r1.9",
                         "CAAS_Psat_ZW6_1.0", "phavu.G19833.gnm2.ann1.PB8d.protein"))

# Accomodate changes in gff annotation files
LOCUS_tissue_spec <- LOCUS_tissue_spec %>%
  mutate(Locus_id = case_when(
    Organism == "Pisum sativum ZW6" ~ Locus,
    Organism == "Lathyrus sativus" ~ coded_by,
    Organism != "Medicago truncatula" ~ gsub("[.][^.]+$", "", seqname), # remove everything after the last .
    TRUE ~ seqname
  ))
LOCUS_tissue_spec <- LOCUS_tissue_spec %>%
  # copy the info in coded_by when the phaseolus vulgaris
  #locus is incomplete
  mutate(Locus = case_when(
    Locus == "Phvul" ~ str_sub(LOCUS_tissue_spec$coded_by, end = -3),
    TRUE ~ Locus))

# Treat alternative splicing variants, keeping only the longest
LOCUS_tissue_spec <- LOCUS_tissue_spec %>%
  # keep the most complete variant (longer length)
  group_by(Locus_id) %>%
  arrange(desc(Length)) %>%
  slice(1)
  #337847 to #280756
  
counts_all <- counts_all %>%
  group_by(Gene_id, Species) %>% 
  #select(-c(gene_id, Gene_id, Species)) %>%
  summarise_at(which(sapply(counts_all, is.numeric)),.funs = sum,na.rm=T)
write.csv(counts_all, "All_species_counts_without_MSTRG_not_filtered_or_normalized_112024_grouped_by_Gene_id.csv", row.names=FALSE )

counts_all_LOCUS <- counts_all %>%
  left_join(LOCUS_tissue_spec, by= c("Gene_id" = "Locus_id")) 

# There's many P. vul missing from the Locus
# Removing the Splicing alternatives to get the Loci 
#(i.e. remove every character after the last dot)
counts_all_LOCUS <- counts_all_LOCUS %>%
  mutate(Organism = case_when(
    Species == "Phaseolus vulgaris" ~ "Phaseolus vulgaris cv. G19833",
    TRUE ~ Organism
  )) %>%
  mutate(Genome = case_when(
    Species == "Phaseolus vulgaris" ~ "phavu.G19833.gnm2.ann1.PB8d.protein",
    TRUE ~ Genome
  )) %>%
  mutate(coded_by= case_when(
    # remove prefix of seqname
    (Species == "Phaseolus vulgaris" & is.na(coded_by)) ~ 
      str_replace_all(seqname, "phavu.G19833.gnm2.ann1.", "") ,
    TRUE ~ coded_by
  )) %>%
  mutate(Locus = case_when(
    # remove characters after the last dot
    (Species == "Phaseolus vulgaris" & is.na(Locus)) ~ sub(".[^.]+$", "", coded_by),
    TRUE ~ Locus
  ))

# Pisum sativum has LOC instead of XP in Gene_id, 
# and Lathyrus sativus has LATHSAT_LOCUS instead of gXXXXX,
#so let's correct that before joining dataframes
counts_all_LOCUS_final <- counts_all_LOCUS %>%
  drop_na(Organism) %>%
  # store transcript info for Lathyrus sativus
  mutate(gene = seqname) %>%
  mutate(Gene_id = case_when(
    Species == "Pisum sativum" ~ seqname,
    # I have to remove the transcript info from the Lathyrus genes
    Species == "Lathyrus sativus" ~ sub(".[^.]+$", "", seqname),
    TRUE ~ Gene_id
  ))

# This file has one transcript per locus
write.csv(counts_all_LOCUS_final, "counts_all_LOCUS_final.csv", row.names=FALSE)

##### BUSCO PCA #####
counts_all_LOCUS_final <- read.csv("counts_all_LOCUS_final.csv")
BUSCO_all <- read.csv("BUSCO_ids_Tissue_Spec_112024.csv", row.names = 1)

BUSCO_counts_all_2 <- left_join(BUSCO_all, counts_all_LOCUS_final, by=c("Species", 
                                                          "gene_id"="Gene_id"
                                                          ))

setwd("../Tissue_Spec/Counts_2024/112024")
write.csv(BUSCO_counts_all_2, "BUSCO_counts_all_duplicated.csv")
BUSCO_counts_all_2 <- read.csv("BUSCO_counts_all_duplicated.csv", row.names=1)
BUSCO_counts_all_2 <- BUSCO_counts_all_2 %>%
  group_by(Busco_id, Species) %>%
  arrange(desc(seqname)) %>%
  arrange(desc(Length.x)) %>% # keep the larger
  slice(1) #37170/7 = 5310 BUSCOs; #37562/7 = 5366 BUSCOs

sum(BUSCO_counts_all_2$Status == "Missing") # 1189

write.csv(BUSCO_counts_all_2, "BUSCO_counts_all.csv")

# Let's make a PCA plot
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
#3503/5366. This is a limitation of this methodology.
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
setwd("../Tissue_Spec/Counts_2024/112024")

library(heatmaply)
#counts_BUSCO <- read.csv("BUSCO_counts_all_counts_pivot.csv", row.names=2)
#counts_BUSCO$X <- NULL
counts_BUSCO <- read.csv(paste("DESeq/Raw/deseq2_BUSCO_TMM.csv", sep="_"), row.names = 1)
# Filter rowmean > 5
counts_BUSCO_filtered <- counts_BUSCO
#counts_BUSCO_filtered$rowmean <- rowMeans(counts_BUSCO_filtered, na.rm=TRUE)
#counts_BUSCO_filtered <- subset(counts_BUSCO_filtered, rowmean > 5000)
#counts_BUSCO_filtered <- counts_BUSCO_filtered[,c(1:42)] #remove the rowmean column

#to perform a log transformation, I need to remove 0s; therefore, I will add one to all counts
counts_BUSCO_filtered <- counts_BUSCO_filtered + 1
counts_BUSCO_filtered <- counts_BUSCO_filtered %>%
  drop_na()
counts_BUSCO_filtered_more <- counts_BUSCO_filtered
counts_BUSCO_filtered_more$rowmean <- rowMeans(counts_BUSCO_filtered_more, na.rm=TRUE)
counts_BUSCO_filtered_more <- subset(counts_BUSCO_filtered_more, rowmean > 250)
counts_BUSCO_filtered_more <- counts_BUSCO_filtered_more[,c(1:42)] #remove the rowmean column

rc <- colorspace::rainbow_hcl(ncol(counts_BUSCO_filtered))
# repeat six times because we have six samples per species
rc <- c(rep(c("#003b84","#ff7400","#00b67e","#488200","#FFBB00","#3f0051","#FF1100"),
            times=1, each=6))
heatmaply(log(counts_BUSCO_filtered_more),xlab = "samples", 
          fontsize_row = 6, fontsize_col = 9,
          colors = viridis(256, option = "mako"), grid_gap = 0.1, Rowv = TRUE,
          label_names = c("Tissue/Accession", "Gene", "Expression Group"),
          key.title = "log(Expression of BUSCOs)",
          col_side_colors = rc, file = "DESeq/Output/heatmaply_BUSCO_clustering_250_TMM_3.pdf"
          )
# 10203at72025 is a BUSCO with high leaf expression and low root expression in all
#counts_BUSCO[counts_BUSCO==0] <- NA #to perform a log transformation, I need to remove 0s
counts_BUSCO_new <- counts_BUSCO + 1
#counts_BUSCO <- counts_BUSCO %>%
#  drop_na()
  
heatmaply(log(counts_BUSCO_new),xlab = "samples", 
          fontsize_row = 6, fontsize_col = 9,
          colors = viridis(256, option = "mako"), grid_gap = 0.1, Rowv = TRUE,
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
######### DESeq2 Pipeline - one species at a time ########### 
# set working directory
setwd("../Tissue_Spec/Counts_2024/112024")

# Load metadata
metadata <- read_xlsx("RNASeq_Tissue_metadata.xlsx", 
                      #header = TRUE, sep = ",", 
)
metadata$Tissue <- relevel(factor(metadata$Tissue), ref = 'Root') #make root as ref

#choose species
choice <- menu(c("Car", "Gma", "Lcu", "Lsa", "Mtr", "Psa", "Pvu"
                 , "all"
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
} else if (choice == 8) {
  species <- "all"
  species_full <- "all"
  col_data <- data.frame(ID=c("L1","L2","L3","R1","R2","R3"), 
                         Tissue=c("Leaf","Leaf","Leaf","Root","Root","Root"))
  col_data$Tissue <- relevel(factor(col_data$Tissue), ref = 'Root') #make root as ref
  cat ("You chose to check all species")
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

print(col_data$ID) #check samples
#perform differential expression analysis
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
resNorm_Tissue_df <- read.csv(paste("DESeq/Sum/deseq2", species, "All_genes.csv", sep = "_"), row.names = 1)

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

##### Volcano plot with all reads (not DEG only) - one species at a time #####
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

##### LogFC histogram - one species at a time #######
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
     xlab="log(Fold Change)", xlim=c(-8,8) )
dev.off()

table <- table %>%
  mutate(Species = case_when(
    startsWith(as.character(ID), "LOC") ~ "Pisum sativum",  
    startsWith(as.character(ID), "Mt") ~ "Medicago truncatula",  
    startsWith(as.character(ID), "phavu") ~ "Phaseolus vulgaris",  
    startsWith(as.character(ID), "Lcu") ~ "Lens culinaris",  
    startsWith(as.character(ID), "Glyma") ~ "Glycine max",  
    startsWith(as.character(ID), "LATH") ~ "Lathyrus sativus",  
    startsWith(as.character(ID), "cicar") ~ "Cicer arietinum",  
    TRUE ~ NA))
table <- table %>%
  mutate(NLR = case_when(
    ID %in% TotalClassNLRs$seqname ~ "NLR",
    ID %in% TotalClassNLRs$original ~ "NLR",
    ID %in% TotalClassNLRs$Locus ~ "NLR",
    TRUE ~ "Other gene"
      ))

table$Species = factor(table$Species,
                       levels =c("Phaseolus vulgaris","Glycine max",
                                 "Cicer arietinum","Lens culinaris",
                                 "Pisum sativum","Lathyrus sativus", 
                                 "Medicago truncatula"
                                 ), ordered =TRUE) 

ggplot(table, aes(x=log2FoldChange, fill=diffexpressed)) +
  geom_histogram(position="identity", #colour="grey40", 
                 #alpha=0.2, 
                 bins = 100) + xlim(c(-5,5)) + #ylim(c(0,1000)) +
  theme_bw() + coord_flip() +
  scale_fill_manual(values=c(NO="grey", Root="brown4", Leaf = "olivedrab3")) + 
  #geom_abline(slope = 1,intercept = c(-1,0,1),color=c("brown","black","olivedrab3"),linetype="dashed") +
  facet_grid(Species ~ .,
             scales = "free_x", # different scales per graph
             space = "free") + # graphs can occupy different spaces
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "none")
ggsave("DESeq/Output/Histograms_All_genes_logFC_per_Species_2.png", width=500, height=720, units='px', dpi = 125)
ggsave("DESeq/Output/Histograms_All_genes_logFC_per_Species_2.svg", width=500, height=720, units='px', dpi = 125)


ggplot(table_NLRs, aes(x=log2FoldChange, fill=diffexpressed)) +
  geom_histogram(position="identity", #colour="grey40", 
                 #alpha=0.2, 
                 bins = 100) + xlim(c(-5,5)) + #ylim(c(0,1000)) +
  theme_bw() + coord_flip() +
  scale_fill_manual(values=c(NO="grey", Root="brown4", Leaf = "olivedrab3")) + 
  #geom_abline(slope = 1,intercept = c(-1,0,1),color=c("brown","black","olivedrab3"),linetype="dashed") +
  facet_grid(Species ~ .,
             scales = "free_x", # different scales per graph
             space = "free") + # graphs can occupy different spaces
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "none")
ggsave("DESeq/Output/Histograms_All_NLRs_logFC_per_Species_2.png", width=500, height=720, units='px', dpi = 125)
ggsave("DESeq/Output/Histograms_All_NLRs_logFC_per_Species_2.svg", width=500, height=720, units='px', dpi = 125)

table_NLRs <- table %>%
  filter(NLR == "NLR")
#table(table_NLRs$diffexpressed) # check how many NLRs were differentially expressed
#table(table_NLRs[table_NLRs$Species=="Pisum sativum",]$diffexpressed)
table_NLRs_Class <- full_join(table_NLRs, TotalClassNLRs, by=c("ID"="Gene_id"))
table_NLRs_Class <- table_NLRs_Class %>%
  filter(Class %in% c("CC-NLR", "TIR-NLR", "CCG10-NLR", "CCR-NLR")) %>%
  drop_na(Species.x)
ggplot(table_NLRs_Class, aes(x=log2FoldChange.x, fill=diffexpressed)) +
  geom_histogram(position="identity", #colour="grey40", 
                 #alpha=0.2, 
                 bins = 40) + xlim(c(-5,5)) + #ylim(c(0,1000)) +
  theme(strip.text.y = element_text(angle = 0)) + theme_bw() + coord_flip() +
  scale_fill_manual(values=c(NO="grey", Root="brown4", Leaf = "olivedrab3")) + 
  #geom_abline(slope = 1,intercept = c(-1,0,1),color=c("brown","black","olivedrab3"),linetype="dashed") +
  facet_grid(Species.x ~ Class,
             scales = "free_x", # different scales per graph
             space = "free") + # graphs can occupy different spaces
  theme(strip.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0),
        legend.position = "none")
ggsave("DESeq/Output/Histograms_All_NLRs_logFC_per_Species_class.png", width=1000, height=720, units='px', dpi = 125)
ggsave("DESeq/Output/Histograms_All_NLRs_logFC_per_Species_class.svg", width=1000, height=720, units='px', dpi = 125)



##### Venn Diagram Leaf vs Root - one species at a time ######################
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

##### Poisson Distance between samples and QC - one species at a time ############
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

##### Histogram based on Daniel Luedke's root-specificity paper - one per species #####
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

png(filename = paste("DESeq/Output/Histogram", species,"logFC.png", sep="_"),
    width = 1250,height = 720, units = "px"
    #, res = 600
)
hist(resNorm_Tissue_df$log2FoldChange, breaks=250, #ylim=c(0,1500), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border=F,#"darkgrey", 
     col = my_color#, show.value=TRUE
)
abline(v = c(-1,0,1), 
       col=c("brown","black","olivedrab3"), 
       lwd=1, lty=2)
dev.off()

##### Single gene plot counts based on Daniel Luedke's root-specificity paper - one per species ##########
#generates a plot for the representation normalized expression for root and shoot tissue
plotCounts(dds, gene="LATHSAT_LOCUS18167", intgroup = "Tissue", returnData = TRUE) %>% 
  ggplot() + aes(Tissue, count, color=Tissue) + 
  geom_jitter(width = 0.2, size = 4, shape = 4, stroke = 3) + ggtitle("Gene of Interest") +
  theme_bw()+scale_y_continuous(limits = c(0, 10), oob = scales::squish)

######### After doing all species: Total NLRs ########################
# Read file with all NLRs for each species. 
#The seqname does not have alternative splicing
TotalNLRs <- read.csv("Locus and Type Legume NLRs 28102024.csv", sep=",")

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
  filter(ID %in% TotalNLRs$Locus) %>%
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


##### NLRs logFC dotplot #####################

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
#Species == "Phaseolus vulgaris" ~ "#FF1100",
#Species == "Glycine max" ~ "#ff7400",
#Species == "Cicer arietinum" ~ "#FFBB00",
#Species == "Medicago truncatula" ~ "#488200",
#Species == "Lens culinaris" ~ "#00b67e",
#Species == "Pisum sativum" ~ "#003b84",
#Species == "Lathyrus sativus" ~ "#3f0051",
ggplot(NLRs, aes(Species, log2FoldChange)) +
  geom_boxplot(colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
                        "#FF1100", "#003b84"))+
  geom_jitter(colour=NLRs$cond, shape=1, size=1) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), 
             color = "black", linewidth=0.5)
ggsave("DESeq/Output/NLRs/Total_NLRs_logFC_per_Species_jitter3.png", width=1250, height=720, units='px', dpi = 125)
ggsave("DESeq/Output/NLRs/Total_NLRs_logFC_per_Species_jitter3.svg", width=1250, height=720, units='px', dpi = 125)

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

# Calculate histogram, but do not draw it
my_hist=hist(NLRs$log2FoldChange , breaks=250  , plot=F)
# Color vector according to log2FC
my_color= ifelse(my_hist$breaks < -1, "brown" , ifelse (my_hist$breaks >=1, "olivedrab3", rgb(0.2,0.2,0.2,0.2) ))

#Fold change distribution. Shows a histogram with the frequency of the fold change distribution for a subset of genes e.g. all NLRs as definde in the subsetting step
png(filename = "DESeq/Output/Histogram_All_NLRs_logFC_2.png",
    width = 1250,height = 720, units = "px"
    #, res = 600
    )
hist(NLRs$log2FoldChange, breaks=250, #ylim=c(0,1500), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border=F,#"darkgrey", 
     col = my_color#, show.value=TRUE
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


##### NLR Expression ##########################
deseq_All <- read.csv("DESeq/Raw/All_species_counts_without_MSTRG_not_filtered_or_normalized_112024_grouped_by_Gene_id.csv"#, 
                      #row.names =1
                      )
# sum the reads with the same gene ID 
deseq_All <- deseq_All %>%
  group_by(Species, Gene_id) %>% 
  summarise_all(.funs = sum,na.rm=T) #557146 to 557146
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
NLRs_expression <- read.csv("DESeq/Raw/All_species_NLR_expression.csv", row.names=1)

#NLRs <- NLRs %>%
#  mutate(cond = case_when(
#    log2FoldChange>1 ~ 'olivedrab3',
#    log2FoldChange<-1 ~ 'brown4',
#    TRUE ~ 'black'   #anything that does not meet the criteria above
#  ))

ggplot(NLRs, aes(Species, log2FoldChange)) +
  geom_jitter(colour=NLRs$cond) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), color = "black", linewidth=0.5)
#geom_point(data = ds, aes(y = mean), colour = 'red', size = 3)
ggsave("DESeq/Output/NLRs/Total_NLRs_logFC_jitter.png", width=1250, height=720, units='px', dpi = 125)
ggsave("DESeq/Output/NLRs/Total_NLRs_logFC_jitter.svg", width=1250, height=720, units='px', dpi = 125)

All_info_NLRs <- left_join(NLRs_expression, NLRs, by = c("Species", "Gene_id" = "ID"))
write.csv(All_info_NLRs, file = "DESeq/Raw/All_species_NLR_expression_logFC.csv")

##### NLR graphs per species #####
NLRs_expression <- read.csv("DESeq/Raw/All_species_NLR_expression.csv", row.names=1)

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
ggplot(data, aes(Sample, Gene_id, fill= log(Value+1))) + 
  geom_tile()  +
  scale_fill_gradient2(low="white", mid= "#4393c3", high="#67001f", 
                       midpoint = 
                         log((max(data$Value)+1))/2) + 
  # mid point is the logarithm of half of the max value
  theme(axis.text.y = element_text(size = 5)) #+
#theme_bw()
ggsave(paste("DESeq/Output/NLRs/All_NLRs", species_full, "CPM_Heatmap_log.png", sep="_"), width=1250, 
       height=720, units='px', dpi = 125)
ggsave(paste("DESeq/Output/NLRs/All_NLRs", species_full, "CPM_Heatmap_log.svg", sep="_"), width=1250, 
       height=720, units='px', dpi = 125)
ggplot(data, aes(Sample, Gene_id, fill= Value)) + 
  geom_tile()  +
  scale_fill_gradient2(low="white", mid= "#4393c3", high="#67001f", 
                       midpoint = max(data$Value)/2) + 
  # mid point is the logarithm of half of the max value
  theme(axis.text.y = element_text(size = 5)) #+
#theme_bw()
ggsave(paste("DESeq/Output/NLRs/All_NLRs", species_full, "CPM_Heatmap.png", sep="_"), width=1250, 
       height=720, units='px', dpi = 125)
ggsave(paste("DESeq/Output/NLRs/All_NLRs", species_full, "CPM_Heatmap.svg", sep="_"), width=1250, 
       height=720, units='px', dpi = 125)

# All_NLRs_Cicer_arietinum_CPM_DESeq_Heatmap_bw.png
ggplot(data, aes(Sample, Gene_id, fill= Value)) + 
  geom_tile()  +
  scale_fill_gradient2(low="white", mid= "darkblue", high="black", 
                       midpoint = max(data$Value)/2) + 
  theme(axis.text.y = element_text(size = 5))
ggsave(paste("DESeq/Output/NLRs/All_NLRs", species_full, "CPM_Heatmap_bw.png", sep="_"), width=1250, 
       height=720, units='px', dpi = 125)
ggsave(paste("DESeq/Output/NLRs/All_NLRs", species_full, "CPM_Heatmap_bw.svg", sep="_"), width=1250, 
       height=720, units='px', dpi = 125)
##### New counts ##############
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

##### specific NLR Venn ######
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
##### NLR Class Expression ###############
# Read file with all NLRs for each species. The seqname does not have alternative splicing
TotalNLRs <- read.csv("Locus and Type Legume NLRs 28102024.csv", sep=",", row.names=1)

TSSpecies <- c("Lathyrus sativus", "Pisum sativum", "Lens culinaris", 
               "Cicer arietinum", "Medicago truncatula", "Phaseolus vulgaris",
               "Glycine max")

# filter to have only the NLRs from the species we want
TotalNLRs <- TotalNLRs %>%
  filter(Species %in% TSSpecies) 

new_counts <- read.csv(file="DESeq/Raw/All_species_NLR_expression_Venn_table.csv", row.names=1) 
# Unfortunately, there are species whose annotation gene name is the Locus
# (in Lathyrus sativus and Pisum sativum). Therefore, we need to merge into two
# separate dataframes
TotalClassNLRs_1 <- inner_join(All_info_NLRs, TotalNLRs, by= c("Gene_id"="seqname", "Species"))
TotalClassNLRs_2 <- inner_join(All_info_NLRs, TotalNLRs, by = c("Gene_id"="Locus", "Species"))
TotalClassNLRs <- full_join(TotalClassNLRs_1,TotalClassNLRs_2)
# now we need to fetch the information we lost in the process and slice alternative splicing
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

#TotalClassNLRs <- full_join(TotalClassNLRs, new_counts)
TotalClassNLRs$averageLeaf <- rowMeans(TotalClassNLRs[, 1:3])
TotalClassNLRs$averageRoot <- rowMeans(TotalClassNLRs[, 4:6])
TotalClassNLRs$devLeaf <- apply(TotalClassNLRs[, 1:3], MARGIN =1, FUN = sd)
TotalClassNLRs$devRoot <- apply(TotalClassNLRs[, 4:6], MARGIN =1, FUN = sd)

TotalClassNLRs <- TotalClassNLRs %>%
  mutate(Leaf = case_when(
    averageLeaf>0.5 ~ "Yes", 
    averageLeaf<0.5 ~ "No",
    TRUE ~ "No"
  )) 

TotalClassNLRs <- TotalClassNLRs %>%
  mutate(Root = case_when(
    averageRoot>0.5 ~ "Yes", 
    averageRoot<0.5 ~ "No",
    TRUE ~ "No"
  )) 

write.csv(TotalClassNLRs, "DESeq/Output/NLRs/All_NLRs_Expression_logFC_Locus_Class.csv")

##### Venn Diagram all NLRs ##########################
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

TotalClassNLRs_filter <- TotalClassNLRs %>%
  filter(Class %in% c("CC-NLR", "TIR-NLR", "CCG10-NLR", "CCR-NLR"))

ggplot(TotalClassNLRs_filter, aes(x=log2FoldChange, fill=cond)) +
  geom_histogram(position="identity", #colour="grey40", 
                 #alpha=0.2, 
                 bins = 100) + 
  theme(strip.text.y = element_text(angle = 0)) + theme_bw() + coord_flip() +
  scale_fill_manual(values=c(black="black", brown4="brown4", olivedrab3 = "olivedrab3")) + 
  #geom_abline(slope = 1,intercept = c(-1,0,1),color=c("brown","black","olivedrab3"),linetype="dashed") +
  facet_grid(Species ~ Class,
             scales = "free", # different scales per graph
             space = "free") # graphs can occupy different spaces
  
ggsave("DESeq/Output/NLRs/Histogram_facet_grid_Species_Class_2.png", width=1250, height=720, units='px', dpi = 125)
ggsave("DESeq/Output/NLRs/Histogram_facet_grid_Species_Class_2.svg", width=1250, height=720, units='px', dpi = 125)


CC <- TotalClassNLRs %>%
  filter(Class=="CC-NLR")
# Calculate histogram, but do not draw it
my_hist=hist(CC$log2FoldChange , breaks=250  , plot=F)
# Color vector according to log2FC
my_color= ifelse(my_hist$breaks < -1, "brown" , ifelse (my_hist$breaks >=1, "olivedrab3", rgb(0.2,0.2,0.2,0.2) ))

#Fold change distribution. Shows a histogram with the frequency of the fold change distribution for a subset of genes e.g. all NLRs as definde in the subsetting step
png(filename = "DESeq/Output/Histogram_CC_NLRs_logFC_2.png",
    width = 1250,height = 720, units = "px"
    #, res = 600
)
hist(CC$log2FoldChange, breaks=250, #ylim=c(0,1500), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border=F,#"darkgrey", 
     col = my_color#, show.value=TRUE
)
abline(v = c(-1,0,1), 
       col=c("brown","black","olivedrab3"), 
       lwd=1, lty=2)
dev.off()


TIR <- TotalClassNLRs %>%
  filter(Class=="TIR-NLR")
# Calculate histogram, but do not draw it
my_hist=hist(TIR$log2FoldChange , breaks=250  , plot=F)
# Color vector according to log2FC
my_color= ifelse(my_hist$breaks < -1, "brown" , ifelse (my_hist$breaks >=1, "olivedrab3", rgb(0.2,0.2,0.2,0.2) ))

#Fold change distribution. Shows a histogram with the frequency of the fold change distribution for a subset of genes e.g. all NLRs as definde in the subsetting step
png(filename = "DESeq/Output/Histogram_TIR_NLRs_logFC_2.png",
    width = 1250,height = 720, units = "px"
    #, res = 600
)
hist(TIR$log2FoldChange, breaks=250, #ylim=c(0,1500), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border=F,#"darkgrey", 
     col = my_color#, show.value=TRUE
)
abline(v = c(-1,0,1), 
       col=c("brown","black","olivedrab3"), 
       lwd=1, lty=2)
dev.off()


CCR <- TotalClassNLRs %>%
  filter(Class=="CCR-NLR")
# Calculate histogram, but do not draw it
my_hist=hist(CCR$log2FoldChange , breaks=50  , plot=F)
# Color vector according to log2FC
my_color= ifelse(my_hist$breaks < -1, "brown" , ifelse (my_hist$breaks >=1, "olivedrab3", rgb(0.2,0.2,0.2,0.2) ))

#Fold change distribution. Shows a histogram with the frequency of the fold change distribution for a subset of genes e.g. all NLRs as definde in the subsetting step
png(filename = "DESeq/Output/Histogram_CCR_NLRs_logFC_2.png",
    width = 1250,height = 720, units = "px"
    #, res = 600
)
hist(CCR$log2FoldChange, breaks=50, #ylim=c(0,1500), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border=F,#"darkgrey", 
     col = my_color#, show.value=TRUE
)
abline(v = c(-1,0,1), 
       col=c("brown","black","olivedrab3"), 
       lwd=1, lty=2)
dev.off()


CCG <- TotalClassNLRs %>%
  filter(Class=="CCG10-NLR")
# Calculate histogram, but do not draw it
my_hist=hist(CCG$log2FoldChange , breaks=100  , plot=F)
# Color vector according to log2FC
my_color= ifelse(my_hist$breaks < -1, "brown" , ifelse (my_hist$breaks >=1, "olivedrab3", rgb(0.2,0.2,0.2,0.2) ))

#Fold change distribution. Shows a histogram with the frequency of the fold change distribution for a subset of genes e.g. all NLRs as definde in the subsetting step
png(filename = "DESeq/Output/Histogram_CCG10_NLRs_logFC_2.png",
    width = 1250,height = 720, units = "px"
    #, res = 600
)
hist(CCG$log2FoldChange, breaks=100, #ylim=c(0,1500), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border=F,#"darkgrey", 
     col = my_color#, show.value=TRUE
)
abline(v = c(-1,0,1), 
       col=c("brown","black","olivedrab3"), 
       lwd=1, lty=2)
dev.off()

TNP <- TotalClassNLRs %>%
  filter(Class=="TNP")
# Calculate histogram, but do not draw it
my_hist=hist(TNP$log2FoldChange , breaks=10  , plot=F)
# Color vector according to log2FC
my_color= ifelse(my_hist$breaks < -1, "brown" , ifelse (my_hist$breaks >=1, "olivedrab3", rgb(0.2,0.2,0.2,0.2) ))

#Fold change distribution. Shows a histogram with the frequency of the fold change distribution for a subset of genes e.g. all NLRs as definde in the subsetting step
png(filename = "DESeq/Output/Histogram_TNPs_logFC_2.png",
    width = 1250,height = 720, units = "px"
    #, res = 600
)
hist(TNP$log2FoldChange, breaks=10, #ylim=c(0,1500), 
     xlim=c(-15,15), lwd=2, 
     lend=100, border=F,#"darkgrey", 
     col = my_color#, show.value=TRUE
)
abline(v = c(-1,0,1), 
       col=c("brown","black","olivedrab3"), 
       lwd=1, lty=2)
dev.off()

##### Scatter plot with NLR classes ################
# Scatter plot with groups
class_color <- c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green")

TotalClassNLRs$Class = factor(TotalClassNLRs$Class, 
                                     levels =c("CC-NLR", "TIR-NLR", 
                                               "CCR-NLR", 
                                               "CCG10-NLR",
                                               "TNP", "Other"), ordered =TRUE) 
TotalClassNLRs$Species = factor(TotalClassNLRs$Species,
                                levels =c("Lathyrus sativus", "Pisum sativum",
                                          "Lens culinaris", "Medicago truncatula",
                                          "Cicer arietinum","Glycine max", 
                                          "Phaseolus vulgaris"), ordered =TRUE) 

ggplot(TotalClassNLRs, aes(x = averageLeaf, y = averageRoot, color = Class)) +
  geom_point() +  theme_bw() + ylim(0,60) + #xlim(0,179) + 
  # x=x line
  geom_abline(slope = 1,intercept = 0,color="blue",linetype="dashed") +
  scale_color_manual(values=class_color) + 
  facet_grid(Species~Class,
            scales = "free_x", # different scales per graph
            space = "free")+ # graphs can occupy different spaces
  theme(strip.text.y = element_text(angle = 0)) + stat_smooth()
  #method= lm)
ggsave("DESeq/Output/NLRs/Scatter_plot_facet_grid_Species_Class_2.png", width=1250, height=720, units='px', dpi = 125)


ggplot(TotalClassNLRs, aes(x = averageLeaf, y = averageRoot, color = Class)) +
  geom_point() +  theme_bw() + ylim(0,60) + #xlim(0,179) + 
  # x=x line
  geom_abline(slope = 1,intercept = 0,color="blue",linetype="dashed") +
  scale_color_manual(values=class_color) + 
  facet_grid(Species~Class,
             #scales = "free_x", # different scales per graph
             space = "free")+ # graphs can occupy different spaces
  theme(strip.text.y = element_text(angle = 0)) + stat_smooth()
#method= lm)
ggsave("DESeq/Output/NLRs/Scatter_plot_facet_grid_Species_Class_4.png", width=1250, height=720, units='px', dpi = 125)


NLRs_Lathyrus <- TotalClassNLRs %>%
  filter(Species=="Lathyrus sativus")
ggplot(NLRs_Lathyrus, aes(x = averageLeaf, y = averageRoot, color = Class)) +
  geom_point() + theme_minimal() +
  scale_color_manual(values=class_color)

##### Expression of Conserved NLRs ###############
BLASTp_NLRs <- read_xlsx("BLASTp_NLRs.xlsx")

# Separating the aminoacid sequence from the IDs
BLASTp_NLRs <- BLASTp_NLRs %>%
  separate(NLR_sequence, c("seqname","sequence"), sep=": ")

# Counting how many aminoacids the protein has into a new column called 'Length'
BLASTp_NLRs$sequence <- str_replace_all(BLASTp_NLRs$sequence," ","") # remove extra spaces
BLASTp_NLRs$sequence <- str_replace_all(BLASTp_NLRs$sequence,"-","") # remove gaps -
BLASTp_NLRs <- BLASTp_NLRs %>%
  mutate("Length" = nchar(BLASTp_NLRs$sequence))

# make seqname not have transcript info
BLASTp_NLRs <- BLASTp_NLRs %>%
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
    grepl("_(", seqname, fixed = TRUE) ~ str_sub(BLASTp_NLRs$seqname, end = -5),
    TRUE ~ seqname
  ))

TotalClassNLRs <- TotalClassNLRs %>%
  mutate(Gene_id = case_when(
    Species == "Pisum sativum" ~ seqname,
    TRUE ~ Gene_id
  ))
BLASTp_NLRs_expression <- full_join(BLASTp_NLRs, TotalClassNLRs, 
                                        by = c("seqname"="Gene_id"))
# Tidying the table
BLASTp_NLRs_expression <- BLASTp_NLRs_expression %>%
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
BLASTp_NLRs_expression <- subset(BLASTp_NLRs_expression, select=
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
colnames(BLASTp_NLRs_expression) <- columns

# Create a column with Tissue assignment (NE, L, R, Both)
BLASTp_NLRs_expression <- BLASTp_NLRs_expression %>%
  mutate(Tissue = case_when(
    Leaf_expressed == "Yes" & Root_expressed == "No" ~ "Leaf",
    Leaf_expressed == "No" & Root_expressed == "Yes" ~ "Root",
    Leaf_expressed == "Yes" & Root_expressed == "Yes" ~ "Both",
    Leaf_expressed == "No" & Root_expressed == "No" ~ "Not Expressed",
    TRUE ~ "Invalid"
  ))

write.csv(BLASTp_NLRs_expression, "BLASTp_NLRs_expression.csv")
##### Bar chart with presence/absence of NLRs with Class info ####
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
                                          "Lens culinaris", 
                                          "Cicer arietinum","Glycine max", 
                                          "Phaseolus vulgaris", "Medicago truncatula"), ordered =TRUE) 


ggplot(TotalClassNLRs, aes(Species, log2FoldChange)) +
  geom_boxplot(colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
                        "#FF1100", "#003b84"))+
  geom_jitter(colour=TotalClassNLRs$cond_logFC, shape=1, size = 1.5) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), 
             color = "black", linewidth=0.5)
ggsave("DESeq/Output/NLRs/Total_Class_NLRs_logFC_per_Species_jitter5.png", width=1250, height=720, units='px', dpi = 125)
ggsave("DESeq/Output/NLRs/Total_Class_NLRs_logFC_per_Species_jitter5.svg", width=1250, height=720, units='px', dpi = 125)

TotalClassNLRs$Class <- factor(TotalClassNLRs$Class,
                               levels=c("CC-NLR", "TIR-NLR",
                                        "CCR-NLR", "CCG10-NLR", 
                                        "TNP", "Other"))
ggplot(TotalClassNLRs, aes(Species, log2FoldChange, color=Class)) +
  geom_boxplot(#color=c("turquoise", "salmon", "thistle3", "lightgoldenrod2",
                        #"dodgerblue", "green", "black"
                        #"#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
                        #"#FF1100", "#003b84"
                        #)
)+ scale_color_manual(name='Species',
                     breaks=c('CC-NLR', 'TIR-NLR', 'CCR-NLR',
                              'CCG10-NLR',
                              'TNP', "Other"),
                     values=c('CC-NLR'= "turquoise3", 'TIR-NLR'="salmon3", 
                              'CCG10-NLR'="lightgoldenrod3",
                              'CCR-NLR'="thistle4",
                              'TNP'="dodgerblue",  
                              'Other'="green4")) +
  #geom_jitter(colour=TotalClassNLRs$cond_logFC, shape=1, size = 0.75) + 
  theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), 
             color = "black", linewidth=0.5) + 
  labs(x = "Species", y = "log(Fold Change) in leaves and roots")
ggsave("DESeq/Output/NLRs/Total_Class_NLRs_logFC_per_Species_jitter6.png", width=1250, height=720, units='px', dpi = 125)
ggsave("DESeq/Output/NLRs/Total_Class_NLRs_logFC_per_Species_jitter6.svg", width=1250, height=720, units='px', dpi = 125)

TotalClassNLRs$Tissue_logFC = factor(TotalClassNLRs$Tissue_logFC, 
                                            levels =c("Not Expressed", "Leaf-specific",
                                                      "Leaf tendency", "Same expression", 
                                                      "Root tendency", "Root-specific"), ordered =TRUE) 

ggplot(TotalClassNLRs, aes(x = averageLeaf, y = averageRoot, color = Tissue_logFC)) +
  geom_point() +  theme_bw() + ylim(0,60) + #xlim(0,179) + 
  # x=x line
  geom_abline(slope = 1,intercept = 0,color="blue",linetype="dashed") +
  scale_color_manual(values=c("black", "#488200", "#b6c630", "blue", "#FFAA00", "#834c3b")) + ####*ALTERAR AQUI*
  facet_grid(Species~Class,
             scales = "free_x", # different scales per graph
             space = "free")+ # graphs can occupy different spaces
  theme(strip.text.y = element_text(angle = 0)) + stat_smooth()
#method= lm)
ggsave("DESeq/Output/NLRs/Scatter_plot_facet_grid_Species_Class_3.png", width=1250, height=720, units='px', dpi = 125)


#TotalClassNLRs_Number <- TotalClassNLRs %>%
#  group_by(Species, Tissue) %>%
#  count(Class)
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
#TotalClassNLRs_Number$Tissue = factor(TotalClassNLRs_Number$Tissue, 
#                                            levels =c("Not Expressed", "Both", "Leaf", "Root"), ordered =TRUE) 
TotalClassNLRs_Number$Tissue_logFC = factor(TotalClassNLRs_Number$Tissue_logFC, 
                                      levels =c("Not Expressed", "Leaf-specific",
                                                "Leaf tendency", "Same expression", 
                                                 "Root tendency", "Root-specific"), ordered =TRUE) 

ggplot() +
  geom_bar(data=TotalClassNLRs_Number, aes(y = n, x = Tissue_logFC, fill = Class)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species) + theme_bw() + 
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_number_tendency.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_number_tendency.svg", width = 1250, height = 720, units="px", dpi = 120)


ggplot() +
  geom_bar(data=TotalClassNLRs_Number, aes(y = ClassPercentage, x = Tissue_logFC, fill = Class)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species) + theme_bw() + scale_y_continuous(minor_breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_percentage_tendency.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_percentage_tendency.svg", width = 1250, height = 720, units="px", dpi = 120)

TotalClassNLRs_Number <- TotalClassNLRs_Number %>%
  group_by(Species, Tissue_logFC) %>%
  mutate(Percentage_Tissue = sum(ClassPercentage)) %>%
  mutate(Number_Tissue_logFC = sum(n))

write.csv(TotalClassNLRs_Number, "Total_Class_NLRs_logFC_tendency_Summary.csv")
TotalClassNLRs_Number <- read.csv("Total_Class_NLRs_logFC_tendency_Summary.csv", row.names=1)

# make a table with the summarised NLR expression per species (no class info)
TotalClassNLRs_Number_by_species <- TotalClassNLRs_Number[, 1:5]
TotalClassNLRs_Number_by_species <- TotalClassNLRs_Number_by_species %>%
  group_by(Species, Tissue_logFC) %>%
  summarise(across(c(n, ClassPercentage), sum))
write.csv(TotalClassNLRs_Number_by_species, "TotalClassNLRs_Number_by_species.csv", row.names=FALSE)

TotalClassNLRs_Number_by_species_pivot <- TotalClassNLRs_Number_by_species[, 1:3] 
TotalClassNLRs_Number_by_species_pivot <- TotalClassNLRs_Number_by_species_pivot %>%
  pivot_wider(names_from=Tissue_logFC, values_from=n) %>%
  mutate(Both = `Leaf tendency` + `Root tendency` + `Same expression`)
write.csv(TotalClassNLRs_Number_by_species_pivot, "TotalClassNLRs_Number_by_species_pivot.csv", row.names=FALSE)

#my.lines = data.frame(x=1, y=seq(from=7, to=798, by = 7))
#Conserved_NLRs_pivot[,c(12:13)] <- NULL

ggplot(TotalClassNLRs_Number,   # Draw heatmap-like plot
       aes(Class, Tissue_logFC, fill = n)) +  geom_tile() + 
  #scale_y_discrete(breaks = TotalClassNLRs_Number$n 
                   # , labels = Conserved_NLRs$seqname) 
                   #+ #geom_hline(yintercept = 0.5 + 0:25908, colour = "black", linewidth = 0.5) +
  # scale_z_log10() +
  scale_fill_gradient2(low="white", mid = "#4393c3",high="#67001f", na.value="white"#, trans = "log" 
                       ) + 
  ggtitle("Expression of known NLR homologues", subtitle = "Leaves and Roots") +
  ylab("NLR homologue") + #geom_segment(data=my.lines, aes(x=1, y, xend=6.5, yend=y), linewidth=0.01, inherit.aes=F)
  theme_minimal() +
  # this line allows to divide the graph in NLR homologue groups
  geom_hline(yintercept = seq(from=0.5, to=798, by = 7), linewidth=0.00001)
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_heatmap_tendency.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_heatmap_tendency.svg", width = 1250, height = 720, units="px", dpi = 120)


ggplot() +
  geom_bar(data=TotalClassNLRs_Number, aes(y = n, x = Tissue_logFC, fill = Class)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species) + theme_bw() + 
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
# The NA is probably a truncated NLR, not considered for the NLR study
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_tendency_class.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/All_NLRs_expression_across_species_tendency_class.svg", width = 1250, height = 720, units="px", dpi = 120)

##### Table with NLR summary per species and bar charts - Fig 3 #####
Sum_All_NLRs <- read_xlsx("NLR expression patterns per species input.xlsx", sheet = "All_NLRs")

Sum_All_NLRs_2 <- Sum_All_NLRs %>%
  pivot_longer(cols=c("Expressed in both", "Leaf-specific", "Root-specific", "Not Expressed"),
               names_to = "Tissue expression")
Sum_All_NLRs_wo_NE <- Sum_All_NLRs %>%
  pivot_longer(cols=c("Expressed in both", "Leaf-specific", "Root-specific"),
               names_to = "Tissue expression")
Sum_All_NLRs_2$`Tissue expression` = factor(Sum_All_NLRs_2$`Tissue expression`, 
                      levels =c("Not Expressed", "Expressed in both",
                                "Leaf-specific", "Root-specific" 
                                ), ordered =TRUE)
Sum_All_NLRs_2$Species = factor(Sum_All_NLRs_2$Species, 
                                            levels =c("Medicago truncatula", "Lathyrus sativus", "Pisum sativum",
                                                      "Lens culinaris",
                                                      "Cicer arietinum",
                                                      "Glycine max",
                                                      "Phaseolus vulgaris"), ordered =TRUE)
ggplot(Sum_All_NLRs_2, aes(fill=`Tissue expression`, x=value, y=Species)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("lightgrey", "blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Number of NLRs")
ggsave("DESeq/Output/NLRs/NLR_expression_patterns_per_species_stacked_JK.png", width = 720, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/NLR_expression_patterns_per_species_stacked_JK.svg", width = 720, height = 720, units="px", dpi = 120)


ggplot(Sum_All_NLRs_2, aes(fill=`Tissue expression`, x=value, y=Species)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("lightgrey", "blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Number of NLRs")
ggsave("DESeq/Output/NLRs/NLR_expression_patterns_per_species_stacked_percentages_JK.png", width = 720, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/NLR_expression_patterns_per_species_stacked_percentages_JK.svg", width = 720, height = 720, units="px", dpi = 120)

Sum_All_NLRs_wo_NE$`Tissue expression` = factor(Sum_All_NLRs_wo_NE$`Tissue expression`, 
                                            levels =c("Expressed in both",
                                                      "Leaf-specific", "Root-specific" 
                                            ), ordered =TRUE)
Sum_All_NLRs_wo_NE$Species = factor(Sum_All_NLRs_wo_NE$Species, 
                                levels =c("Medicago truncatula", "Lathyrus sativus", "Pisum sativum",
                                          "Lens culinaris",
                                          "Cicer arietinum",
                                          "Glycine max",
                                          "Phaseolus vulgaris"), ordered =TRUE)
Sum_All_NLRs_wo_NE = ddply(Sum_All_NLRs_wo_NE, .(Species), transform, percent = value/sum(value) * 100)
Sum_All_NLRs_wo_NE = ddply(Sum_All_NLRs_wo_NE, .(Species), transform, pos = (cumsum(value) - 0.5 * value))
Sum_All_NLRs_wo_NE$label = paste0(sprintf("%.0f", Sum_All_NLRs_wo_NE$percent), "%")
ggplot(Sum_All_NLRs_wo_NE, aes(fill=`Tissue expression`, x=value, y=Species)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("blue", "#488200", "#834c3b")) +
  #geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Proportion of NLR expression classes")
ggsave("DESeq/Output/NLRs/NLR_expression_patterns_per_species_stacked_percentages_wo_NE_JK.png", width = 750, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/NLR_expression_patterns_per_species_stacked_percentages_wo_NE_JK.svg", width = 750, height = 720, units="px", dpi = 120)








Sum_Class_NLRs <- read_xlsx("USEFUL NLR expression per species and class.xlsx", sheet = 1)

Sum_Class_NLRs <- Sum_Class_NLRs %>%
  mutate(Tissue_logFC = case_when(
    Tissue_logFC=="Leaf tendency" ~ "Expressed in both",
    Tissue_logFC=="Same expression" ~ "Expressed in both",
    Tissue_logFC=="Root tendency" ~ "Expressed in both",
    TRUE ~ Tissue_logFC
  ))

Sum_Class_NLRs$Tissue_logFC = factor(Sum_Class_NLRs$Tissue_logFC, 
                                            levels =c("Not Expressed", "Expressed in both",
                                                      "Leaf-specific", "Root-specific" 
                                            ), ordered =TRUE)
Sum_Class_NLRs$Species = factor(Sum_Class_NLRs$Species, 
                                levels =c("Medicago truncatula", "Lathyrus sativus", "Pisum sativum",
                                          "Lens culinaris",
                                          "Cicer arietinum",
                                          "Glycine max",
                                          "Phaseolus vulgaris"), ordered =TRUE)
Sum_Class_NLRs_wo_NE <- Sum_Class_NLRs %>%
  filter(Tissue_logFC != "Not Expressed")

ggplot(Sum_Class_NLRs[Sum_Class_NLRs$Class=="CC-NLR",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("lightgrey", "blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Number of CC-NLRs")
ggsave("DESeq/Output/NLRs/CCNLR_expression_patterns_per_species_stacked.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/CCNLR_expression_patterns_per_species_stacked.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot(Sum_Class_NLRs_wo_NE[Sum_Class_NLRs_wo_NE$Class=="CC-NLR",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Percentages of CC-NLRs")
ggsave("DESeq/Output/NLRs/CCNLR_expression_patterns_per_species_stacked_percent.png", width = 550, height = 550, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/CCNLR_expression_patterns_per_species_stacked_percent.svg", width = 550, height = 550, units="px", dpi = 120)

ggplot(Sum_Class_NLRs[Sum_Class_NLRs$Class=="TIR-NLR",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("lightgrey", "blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Number of TIR-NLRs")
ggsave("DESeq/Output/NLRs/TIRNLR_expression_patterns_per_species_stacked.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/TIRNLR_expression_patterns_per_species_stacked.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot(Sum_Class_NLRs_wo_NE[Sum_Class_NLRs_wo_NE$Class=="TIR-NLR",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Percentages of TIR-NLRs")
ggsave("DESeq/Output/NLRs/TIRNLR_expression_patterns_per_species_stacked_percent.png", width = 550, height = 550, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/TIRNLR_expression_patterns_per_species_stacked_percent.svg", width = 550, height = 550, units="px", dpi = 120)

ggplot(Sum_Class_NLRs[Sum_Class_NLRs$Class=="CCR-NLR",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("lightgrey", "blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Number of CCR-NLRs")
ggsave("DESeq/Output/NLRs/CCRNLR_expression_patterns_per_species_stacked.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/CCRNLR_expression_patterns_per_species_stacked.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot(Sum_Class_NLRs_wo_NE[Sum_Class_NLRs_wo_NE$Class=="CCR-NLR",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Percentages of CCR-NLRs")
ggsave("DESeq/Output/NLRs/CCRNLR_expression_patterns_per_species_stacked_percent.png", width = 550, height = 550, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/CCRNLR_expression_patterns_per_species_stacked_percent.svg", width = 550, height = 550, units="px", dpi = 120)

ggplot(Sum_Class_NLRs[Sum_Class_NLRs$Class=="CCG10-NLR",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("lightgrey", "blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Number of CCG10-NLRs")
ggsave("DESeq/Output/NLRs/CCG10NLR_expression_patterns_per_species_stacked.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/CCG10NLR_expression_patterns_per_species_stacked.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot(Sum_Class_NLRs_wo_NE[Sum_Class_NLRs_wo_NE$Class=="CCG10-NLR",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Percentages of CCG10-NLRs")
ggsave("DESeq/Output/NLRs/CCG10NLR_expression_patterns_per_species_stacked_percent.png", width = 550, height = 550, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/CCG10NLR_expression_patterns_per_species_stacked_percent.svg", width = 550, height = 550, units="px", dpi = 120)

ggplot(Sum_Class_NLRs[Sum_Class_NLRs$Class=="TNP",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("lightgrey", "blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Number of TNPs")
ggsave("DESeq/Output/NLRs/TNP_expression_patterns_per_species_stacked.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/TNP_expression_patterns_per_species_stacked.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot(Sum_Class_NLRs_wo_NE[Sum_Class_NLRs_wo_NE$Class=="TNP",], aes(fill=Tissue_logFC, x=n, y=Species)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("blue", "#488200", "#834c3b")) +
  #scale_
  ggtitle("") +
  theme_bw() +
  ylab("") + xlab("Percentages of TNPs")
ggsave("DESeq/Output/NLRs/TNP_expression_patterns_per_species_stacked_percent.png", width = 550, height = 550, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/TNP_expression_patterns_per_species_stacked_percent.svg", width = 550, height = 550, units="px", dpi = 120)

CC_NLRs <- Sum_Class_NLRs %>%
  filter(Class=="CC-NLR")




#Sum_All_NLRs <- Sum_All_NLRs %>%
#  mutate(Expressed = `Expressed in both`+`Leaf-specific`+`Root-specific`)
#Sum_Tendency_NLRs <- read_xlsx("NLR expression patterns per species input.xlsx", sheet = "NLRs_Expressed_in_both")


#Sum_All_NLRs_pivot <- Sum_All_NLRs[,c(1,3,4,5,6)] %>%
#  pivot_longer(!c(Species, Expressed, `Leaf-specific`), 
#               names_to = "Tissue_Expression", values_to = "Number_negative") 
#Sum_All_NLRs_pivot <- Sum_All_NLRs_pivot %>%
#  pivot_longer(!Species, names_to = "Tissue_Expression", values_to = "Number_positive")

#ggplot(Sum_All_NLRs,aes(x=Species))+
  # Plotting the bar plot using geom_bar()
#  geom_bar(aes(y=`Not Expressed`*(-1), fill="Not Expressed"), stat='identity',width=0.8)+
#  geom_bar(aes(y=Expressed, fill="Expressed"), stat='identity',width=0.8)+coord_flip() +
#  scale_fill_manual(values = c("Not Expressed" = "black", "Expressed" = "blue")) +
#  # Specifying the x,y axis labels and tile of the plot
#  xlab("Letters")+ylab("Values")+labs(title="Customized Diverging Bar Plot") +
#  theme_bw() + facet_wrap(~ Species, ncol = 1)  # Separate by subcategory

# I had to manually input data from Sum_All_NLRs to fit the following graph
data <- read_xlsx("NLR expression patterns per species input.xlsx", sheet = "Histogram_ready")
data$Species = factor(data$Species, 
                      levels =c("Medicago truncatula","Lathyrus sativus", "Pisum sativum", 
                                "Lens culinaris", 
                                "Cicer arietinum",
                                "Glycine max", "Phaseolus vulgaris"), ordered =TRUE) 


ggplot(data, aes(x = Species)) +
  geom_bar(aes(y = Number_positive), fill= "blue", stat = "identity", width = 0.5) +
  geom_bar(aes(y = Number_negative), fill="black", stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("blue", "black")) + ylim(c(-350,450)) +
  geom_text(aes(y = Number_positive, 
                label = paste(Number_positive, "; ", round(Percentage_positive), " %", sep="")),  # Labels for Number_positive
            hjust = -0.2,  # Adjust horizontal position (right of the bar)
            size = 3,  # Adjust text size
            color = "black") +
  geom_text(aes(y = Number_negative, 
                label = paste(round(Percentage_negative), " %; ", Number_negative*(-1), sep="")),  # Labels for Number_negative
            hjust = 1.2,  # Adjust horizontal position (left of the bar)
            size = 3,  # Adjust text size
            color = "black") +
  coord_flip() +  # Flip to horizontal bars
  facet_wrap(~ Category, nrow = 1) +  # Separate by subcategory
  labs(y = "NLR Number", 
       fill = "Legend") + 
  theme_bw()
ggsave("DESeq/Output/NLRs/NLR_expression_patterns_per_species.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("DESeq/Output/NLRs/NLR_expression_patterns_per_species.svg", width = 1250, height = 720, units="px", dpi = 120)

##### Create line expression plot where which line is a species #####
# Overlaying two plots
# Review later
BLASTp_NLRs_expression$Clade <- paste(BLASTp_NLRs_expression$Class, 
                                         BLASTp_NLRs_expression$Node, sep="_")
BLASTp_NLRs_expression <- BLASTp_NLRs_expression[-which(BLASTp_NLRs_expression$Clade == "CCG10-NLR_12"), ]
write.csv(BLASTp_NLRs_expression, "BLASTp_NLRs_expression_no_problematic_CCG10.csv")
BLASTp_NLRs_expression <- read.csv("BLASTp_NLRs_expression_no_problematic_CCG10.csv", row.names=1)

BLASTp_NLRs_expression$Clade2 <- gsub("-", "_", BLASTp_NLRs_expression$Clade)
BLASTp_NLRs_expression_2 <- BLASTp_NLRs_expression %>%
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

anova_1 <- aov(log2FC ~ Clade2, data = BLASTp_NLRs_expression_2)
summary(anova_1)
tukey_1 <- TukeyHSD(anova_1)
print(tukey_1)
cld_1 <- multcompLetters4(anova_1, tukey_1) # compact letter display
print(cld_1)
# table with factors and 3rd quantile
Tk_1 <- group_by(BLASTp_NLRs_expression_2, Clade2) %>%
  summarise(mean=mean(log2FC), quant = quantile(log2FC, probs = 0.75)) %>%
  arrange(desc(mean))
# extracting the compact letter display and adding to the Tk table
cld_1 <- as.data.frame.list(cld_1$Clade2)
Tk_1$cld <- cld_1$Letters
print(Tk_1)
#BLASTp_NLRs_expression_2$Clade2 <- factor(BLASTp_NLRs_expression_2$Clade2,
#                                   levels=c("CC", "TIR",
#                                            "CCG10", "CCR", "TNP", "Other"))
library("scales")
ggplot(BLASTp_NLRs_expression_2, aes(x=Clade2, y=log2FC#, fill=Clade2
                                         )) + 
  geom_boxplot() +
  labs(x="NLR Clade", y="log2FC") +
  theme_bw() + 
  #scale_y_continuous(labels = scales::percent) +
  geom_jitter(color=BLASTp_NLRs_expression_2$Species_Color, size=1, alpha=0.9,
              width = 0.2, height = 0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  labs(title = "logFC of the expression of Conserved NLRs", subtitle = "In leaf and root tissue") +
  geom_text(data = Tk_1, aes(x = Clade2, y = quant, label = cld), size = 3, vjust=-1, hjust =-1) #+
  #scale_fill_brewer(c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
  #scale_fill_manual(values=c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
# saving the final figure
ggsave("DESeq/Output/NLRs/Conserved_NLR_nodes_logFC_jitterplot_2.png", width = 1250, height = 720, units="px", dpi = 120)


ggplot(subset(BLASTp_NLRs_expression, Species=="Phaseolus vulgaris"), 
       aes(Clade, log2FC)) + 
  geom_jitter(color = "#FF1100", size=2, width = 0.2, height = 0) + 
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Glycine max"), 
            color = "#ff7400", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Cicer arietinum"), 
            color = "#FFBB00", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Medicago truncatula"), 
            color = "#488200", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Lens culinaris"), 
            color = "#00b67e", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Pisum sativum"), 
            color = "#003b84", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Lathyrus sativus"), 
            color = "#3f0051", size=2, width = 0.2, height = 0) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  labs(title = "logFC of the expression of Conserved NLRs", subtitle = "In leaf and root tissue")
  
ggsave("DESeq/Output/NLRs/Conserved_NLR_nodes_logFC_jitterplot.png", width=1250, 
       height=720, units='px', dpi = 100)

ggplot(subset(BLASTp_NLRs_expression, Species=="Phaseolus vulgaris"), 
       aes(Clade, Leaf_mean_expression)) + 
  geom_jitter(color = "#FF1100", size=2, width = 0.2, height = 0) + 
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Glycine max"), 
            color = "#ff7400", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Cicer arietinum"), 
            color = "#FFBB00", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Medicago truncatula"), 
            color = "#488200", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Lens culinaris"), 
            color = "#00b67e", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Pisum sativum"), 
            color = "#003b84", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Lathyrus sativus"), 
            color = "#3f0051", size=2, width = 0.2, height = 0) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  labs(title = "Expression of Conserved NLRs", subtitle = "In leaf tissue") +
  ylim(c(0,55))


ggsave("DESeq/Output/NLRs/Conserved_NLR_nodes_TMM_Leaf_jitterplot.png", width=1250, 
       height=720, units='px', dpi = 100)

ggplot(subset(BLASTp_NLRs_expression, Species=="Phaseolus vulgaris"), 
       aes(Clade, Root_mean_expression)) + 
  geom_jitter(color = "#FF1100", size=2, width = 0.2, height = 0) + 
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Glycine max"), 
            color = "#ff7400", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Cicer arietinum"), 
            color = "#FFBB00", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Medicago truncatula"), 
            color = "#488200", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Lens culinaris"), 
            color = "#00b67e", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Pisum sativum"), 
            color = "#003b84", size=2, width = 0.2, height = 0) +
  geom_jitter(data = subset(BLASTp_NLRs_expression, Species=="Lathyrus sativus"), 
            color = "#3f0051", size=2, width = 0.2, height = 0) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  labs(title = "Expression of Conserved NLRs", subtitle = "In root tissue") +
  ylim(c(0,55))

ggsave("DESeq/Output/NLRs/Conserved_NLR_nodes_TMM_Root_jitterplot.png", width=1250, 
       height=720, units='px', dpi = 100)
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

##### Prep iTOL file ####
TotalClassNLRs <- read.csv("Total_Class_NLRs_logFC.csv", row.names = 1)

iTOL <- TotalClassNLRs[, c("Leaf.1","Leaf.2","Leaf.3",
                           "Root.1","Root.2","Root.3",
                           "Species","Class", "original","Tissue")]

iTOL <- iTOL %>%
  mutate(color_species = case_when(
    Species == "Phaseolus vulgaris" ~ "#FF1100",
    Species == "Glycine max" ~ "#ff7400",
    Species == "Cicer arietinum" ~ "#FFBB00",
    Species == "Medicago truncatula" ~ "#488200",
    Species == "Lens culinaris" ~ "#00b67e",
    Species == "Pisum sativum" ~ "#003b84",
    Species == "Lathyrus sativus" ~ "#3f0051",
    TRUE ~ "black"
  )) %>%
  mutate(color_class = case_when(
    Class == "CC-NLR" ~ "#40e0d0",
    Class == "TIR-NLR" ~ "#fa8072",
    Class == "CCG10-NLR" ~ "#FFBB00",
    Class == "CCR-NLR" ~ "#cdb5cd",
    Class == "TNP" ~ "#73a8d0",
    Class == "Other" ~ "#86b875",
    TRUE ~ "black"
  )) %>%
  mutate(color_tissue = case_when(
    Tissue == "Leaf" ~ "#488200",
    Tissue == "Root" ~ "#834c3b",
    Tissue == "Both" ~ "#0000ff",
    TRUE ~ "#000000"
  )) %>%
  arrange(original)



iTOL_Leaf <- iTOL %>%
  filter(Tissue == "Leaf")
iTOL_Leaf <- iTOL_Leaf[, c("Leaf.1","Leaf.2","Leaf.3",
                           "Root.1","Root.2","Root.3","Species","Class", 
                           "original","color_species","color_class")]

iTOL_Root <- iTOL %>%
  filter(Tissue == "Root")
iTOL_Root <- iTOL_Root[, c("Leaf.1","Leaf.2","Leaf.3",
                           "Root.1","Root.2","Root.3","Species","Class", 
                           "original","color_species","color_class")]

write.csv(iTOL, "DESeq/Output/iTOL.csv", row.names=FALSE)  
write.csv(iTOL_Leaf, "DESeq/Output/iTOL_Leaf.csv", row.names=FALSE)  
write.csv(iTOL_Root, "DESeq/Output/iTOL_Root.csv", row.names=FALSE)  
