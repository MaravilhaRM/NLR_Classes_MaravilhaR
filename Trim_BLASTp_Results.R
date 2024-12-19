library(tidyverse)
library(readxl)
library(tidyr)
library(ggplot2)

setwd("D:/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/NLR_seqs/BLASTp")
setwd("E:/KamounLab Dropbox/Kamounity folder/Rita/Tissue_Spec/NLR_seqs/BLASTp")

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

# insert headers in the tabular file:
# Headers: Query_sequence_ID	Target_sequence_ID	Percentage_of_identical_positions	Alignment_length	Number_of_mismatches	Number_of_gaps	Start_of_alignment_in_query	End_of_alignment_in_query	Start_of_alignment_in_target	End_of_alignment_in_target	E_value	Bit_score	All_target_sequence_Ids	Raw_score	Number_of_identical_matches	Number_of_positive_scoring_matches	Total_number_of_gaps	Percentage_of_positive_scoring_matches	Query_frame	Target_frame	Aligned_query_sequence	Align_target_sequence	Query_sequence_length	Target_sequence_length	Target_description
# read BLASTp output with headers
BLASTp <- read.csv(file=paste("BLASTp_", species, ".tabular", 
                              sep=""), 
                   header=TRUE, sep="\t")

BLASTp <- BLASTp %>%
  # group by NLR name
  group_by(Query_sequence_ID) %>%
  # remove all target hits with a description with an Arabidopsis gene or uncharacterized protein
  filter(!grepl('uncharacterized',
    #'At|uncharacterized', # when removing Arabidopsis genes, I lose NLRs and better hits 
    Target_description)) %>%
  # order by lowest e-value, then by largest alignment length and then by percent identity
  arrange(Query_sequence_ID, E_value, desc(Alignment_length), desc(Percentage_of_identical_positions))

# keep one hit per NLR
BLASTp <- BLASTp %>%
  group_by(Query_sequence_ID) %>%
  slice(1)

#test <- BLASTp$Query_sequence_ID %ni% BLASTp_Car$Query_sequence_ID
#test <- BLASTp %>%
#  filter(Query_sequence_ID %ni% BLASTp_Car$Query_sequence_ID)
# Car: 79913 -> 44135 (60735)-> 107 (110)
# Gma: 654392 -> (464022) -> (866)
# Lcu: 154892 -> 107644 (138337) -> 234 (238)
# Lsa: 126455 -> 92665 -> 210 (218)
# Mtr: 463870 -> 320238 (388584) -> 709 (724)
# Psa: 250323 -> 160618 (203710) -> 369 (380)
# Pvu: 312113 -> 197185 (250652) -> 442 (456)

write.csv(BLASTp, paste("BLASTp_", species, "_best_hits_w_At.csv", 
                        sep=""), row.names = FALSE)
# Look manually at the descriptions, insert a column with the NLR homologue;
# if there is no homologue, look at the BLASTp output table and take the next
# best output (provided e-value is very low and percent id is acceptable)
# if there still aren't no homologues, just put NA
All_BLASTp <- BLASTp
###### Join all BLASTp outputs ######
BLASTp_Car <- read_xlsx("BLASTp_Car_best_hits_w_At_edited.xlsx")
BLASTp_Gma <- read_xlsx("BLASTp_Gma_best_hits_w_At_edited.xlsx")
BLASTp_Lcu <- read_xlsx("BLASTp_Lcu_best_hits_w_At_edited.xlsx")
BLASTp_Lsa <- read_xlsx("BLASTp_Lsa_best_hits_w_At_edited.xlsx")
BLASTp_Mtr <- read_xlsx("BLASTp_Mtr_best_hits_w_At_edited.xlsx")
BLASTp_Psa <- read_xlsx("BLASTp_Psa_best_hits_w_At_edited.xlsx")
BLASTp_Pvu <- read_xlsx("BLASTp_Pvu_best_hits_w_At_edited.xlsx")

All_BLASTp <- rbind(BLASTp_Car, BLASTp_Gma, BLASTp_Lcu, BLASTp_Lsa, BLASTp_Mtr,
                    BLASTp_Psa, BLASTp_Pvu)

All_BLASTp <- All_BLASTp %>%
  # create a Like column
  mutate(Isoform = case_when(
    grepl("isoform", Target_description) ~ substr(sub(".*isoform ", "", Target_description), start = 1, 2),
    TRUE ~ Isoform
    )) %>%
  # create a Isoform column (extract what is after isoform)
  mutate(Like = case_when(
    grepl("like", Target_description) & !grepl("RPP13-like protein", Target_description) ~ "like",
    TRUE ~ Like
  )) 

# merge some redundant blastp outputs
All_BLASTp <- All_BLASTp %>%
  mutate(Best_hit = case_when(
    grepl("At1g61180", Best_hit) ~ "UNI",
    grepl("At1g12280", Best_hit) ~ "SUMM2",
    grepl("At1g15890", Best_hit) ~ "L3",
    grepl("At3g14460", Best_hit) ~ "LRRAC1",
    grepl("At5g66900", Best_hit) ~ "NRG1",
    grepl("At4g33300", Best_hit) ~ "ADR1-L1",
    grepl("At5g47280", Best_hit) ~ "ADR1-L3",
    grepl("At4g11170", Best_hit) ~ "RMG1",
    grepl("At5g43730", Best_hit) ~ "RSG2",
    grepl("At1g51480", Best_hit) ~ "RSG1",
    TRUE ~ Best_hit
  ))

write.csv(All_BLASTp, "All_BLASTp.csv")
All_BLASTp <- read.csv("All_BLASTp.csv", row.names=1)

# filter out homologues with few hits
All_BLASTp_count <- All_BLASTp %>%
  #group_by(Species) %>%
  count(Best_hit)  

# save the NLRs with few homologues (<7) in a variable
Few_homologues <- All_BLASTp_count[All_BLASTp_count$n < 7, ]$Best_hit

#All_BLASTp_count <- All_BLASTp_count %>%
#  filter()
###### Join BLASTp outputs with expression+loci data #####
Tissue_Class_NLRs_logFC <- read.csv("Total_Class_NLRs_logFC.csv", row.names=1)
Tissue_Class_NLRs_logFC_BLAST <- left_join(Tissue_Class_NLRs_logFC, All_BLASTp,
                                           by = c("original"="Query_sequence_ID"))
# Create a function that does the opposite of %in%
`%ni%` <- Negate(`%in%`)

Tissue_Class_NLRs_logFC_BLAST <- Tissue_Class_NLRs_logFC_BLAST %>%
  # remove NLRs without characterized BLAST outputs
  drop_na(Best_hit) %>%
  filter(Best_hit %ni% Few_homologues) %>%
  mutate(Species_colour = case_when(
    Species == "Phaseolus vulgaris" ~ "#FF1100",
    Species == "Glycine max" ~ "#ff7400",
    Species == "Cicer arietinum" ~ "#FFBB00",
    Species == "Medicago truncatula" ~ "#488200",
    Species == "Lens culinaris" ~ "#00b67e",
    Species == "Pisum sativum" ~ "#003b84",
    Species == "Lathyrus sativus" ~ "#3f0051",
    TRUE ~ "black"
  )) %>%
  mutate(Like_shape = case_when(
    Like == "like" ~ 1,
    TRUE ~ 19
  )) %>%
  mutate(dot_size=case_when(
    log2FoldChange > 1 ~ log(averageLeaf + 1),
    log2FoldChange < -1 ~ log(averageRoot + 1),
    TRUE ~ log(baseMean + 1)
    ))
  
Tissue_Class_NLRs_logFC_BLAST$Species <- factor(Tissue_Class_NLRs_logFC_BLAST$Species,
                                 levels=c("Phaseolus vulgaris", "Glycine max",
                                          "Cicer arietinum", "Medicago truncatula", 
                                          "Lens culinaris", "Pisum sativum",
                                          "Lathyrus sativus"))
write.csv(Tissue_Class_NLRs_logFC_BLAST, "Tissue_Class_NLRs_logFC_BLAST.csv")
# Idea: 
# make a combined jitter plot with NLR homologues by class (CC, TIR, CCG10, CCR), 
# coloring each dot according to species
# the value would be either logFC or expression in leaves/roots

# filter by class
Tissue_Class_NLRs_logFC_BLAST_CC <- Tissue_Class_NLRs_logFC_BLAST %>%
  filter(Class=="CC-NLR") %>%
  filter(Best_hit %ni% c("RUN1", "At5g47260")) %>%
  mutate(Like = case_when(
    grepl("RPP13-like", Best_hit) ~ "like",
    TRUE ~ Like
  )) %>%
  mutate(Like_shape = case_when(
    Like == "like" ~ 1,
    TRUE ~ 19
  )) %>%
  mutate(Best_hit = case_when(
    grepl("RPP13-like", Best_hit) ~ "RPP13",
    TRUE ~ Best_hit
  ))
Tissue_Class_NLRs_logFC_BLAST_TIR <- Tissue_Class_NLRs_logFC_BLAST %>%
  filter(Class=="TIR-NLR") %>%
  filter(Best_hit %ni% c("protein")) # these are not true homologues
Tissue_Class_NLRs_logFC_BLAST_CCR <- Tissue_Class_NLRs_logFC_BLAST %>%
  filter(Class=="CCR-NLR")
Tissue_Class_NLRs_logFC_BLAST_CCG10 <- Tissue_Class_NLRs_logFC_BLAST %>%
  filter(Class=="CCG10-NLR") %>%
  filter(Best_hit %ni% c("RGA1","RGA2", "RGA3", "TMV resistance protein N", 
                         "At5g05400","At5g47260", "RF9"))

##### CC-NLRs #####

ggplot(Tissue_Class_NLRs_logFC_BLAST_CC, aes(Best_hit, log2FoldChange#, color=Species
                                             )) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_CC$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_CC$Like_shape, 
              size=Tissue_Class_NLRs_logFC_BLAST_CC$dot_size) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), 
             color = "black", linewidth=0.5) + ylim(c(-10, 10)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  labs(x = "CC-NLR homologue", y = "log(Fold Change) in leaf and root tissue")
ggsave("Known_NLRs_logFC_per_Species_jitter_CC_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_logFC_per_Species_jitter_CC_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)

ggplot(Tissue_Class_NLRs_logFC_BLAST_CC, aes(Best_hit, averageLeaf#, color=Species
                                             )) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+ ylim(c(0,65)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_CC$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_CC$Like_shape, 
              size=log(Tissue_Class_NLRs_logFC_BLAST_CC$averageLeaf + 1)) + theme_bw() + 
  labs(x = "CC-NLR homologue", y = "Expression in leaf tissue (TPM)")
  #geom_hline(#yintercept=c(-1,0,1), 
             #linetype=c("dashed","solid","dashed"), 
             #color = "black", linewidth=0.5)
ggsave("Known_NLRs_averageLeaf_per_Species_jitter_CC_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_averageLeaf_per_Species_jitter_CC_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)

ggplot(Tissue_Class_NLRs_logFC_BLAST_CC, aes(Best_hit, averageRoot#, color=Species
                                             )) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+ ylim(c(0,65)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_CC$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_CC$Like_shape, 
              size=log(Tissue_Class_NLRs_logFC_BLAST_CC$averageRoot + 1)) + theme_bw() + 
  labs(x = "CC-NLR homologue", y = "Expression in root tissue (TPM)")
  #geom_hline(#yintercept=c(-1,0,1), 
             #linetype=c("dashed","solid","dashed"), 
             #color = "black", linewidth=0.5)
ggsave("Known_NLRs_averageRoot_per_Species_jitter_CC_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_averageRoot_per_Species_jitter_CC_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)

##### TIR-NLRs #####

ggplot(Tissue_Class_NLRs_logFC_BLAST_TIR, aes(Best_hit, log2FoldChange#, color=Species
)) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_TIR$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_TIR$Like_shape, 
              size=Tissue_Class_NLRs_logFC_BLAST_TIR$dot_size) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), 
             color = "black", linewidth=0.5) + ylim(c(-10, 10)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  labs(x = "TIR-NLR homologue", y = "log(Fold Change) in leaf and root tissue")
ggsave("Known_NLRs_logFC_per_Species_jitter_TIR_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_logFC_per_Species_jitter_TIR_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)

ggplot(Tissue_Class_NLRs_logFC_BLAST_TIR, aes(Best_hit, averageLeaf#, color=Species
)) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+ ylim(c(0,65)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_TIR$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_TIR$Like_shape, 
              size=log(Tissue_Class_NLRs_logFC_BLAST_TIR$averageLeaf + 1)) + theme_bw() + 
  labs(x = "TIR-NLR homologue", y = "Expression in leaf tissue (TPM)")
#geom_hline(#yintercept=c(-1,0,1), 
#linetype=c("dashed","solid","dashed"), 
#color = "black", linewidth=0.5)
ggsave("Known_NLRs_averageLeaf_per_Species_jitter_TIR_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_averageLeaf_per_Species_jitter_TIR_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)

ggplot(Tissue_Class_NLRs_logFC_BLAST_TIR, aes(Best_hit, averageRoot#, color=Species
)) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+ ylim(c(0,65)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_TIR$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_TIR$Like_shape, 
              size=log(Tissue_Class_NLRs_logFC_BLAST_TIR$averageRoot + 1)) + theme_bw() + 
  labs(x = "TIR-NLR homologue", y = "Expression in root tissue (TPM)")
#geom_hline(#yintercept=c(-1,0,1), 
#linetype=c("dashed","solid","dashed"), 
#color = "black", linewidth=0.5)
ggsave("Known_NLRs_averageRoot_per_Species_jitter_TIR_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_averageRoot_per_Species_jitter_TIR_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)
##### CCR-NLRs #####

ggplot(Tissue_Class_NLRs_logFC_BLAST_CCR, aes(Best_hit, log2FoldChange#, color=Species
)) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_CCR$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_CCR$Like_shape, 
              size=Tissue_Class_NLRs_logFC_BLAST_CCR$dot_size) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), 
             color = "black", linewidth=0.5) + ylim(c(-10, 10)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  labs(x = "CCR-NLR homologue", y = "log(Fold Change) in leaf and root tissue")
ggsave("Known_NLRs_logFC_per_Species_jitter_CCR_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_logFC_per_Species_jitter_CCR_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)

ggplot(Tissue_Class_NLRs_logFC_BLAST_CCR, aes(Best_hit, averageLeaf#, color=Species
)) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+ #ylim(c(0,25)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_CCR$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_CCR$Like_shape, 
              size=log(Tissue_Class_NLRs_logFC_BLAST_CCR$averageLeaf + 1)) + theme_bw() + 
  labs(x = "CCR-NLR homologue", y = "Expression in leaf tissue (TPM)")
#geom_hline(#yintercept=c(-1,0,1), 
#linetype=c("dashed","solid","dashed"), 
#color = "black", linewidth=0.5)
ggsave("Known_NLRs_averageLeaf_per_Species_jitter_CCR_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_averageLeaf_per_Species_jitter_CCR_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)

ggplot(Tissue_Class_NLRs_logFC_BLAST_CCR, aes(Best_hit, averageRoot#, color=Species
)) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+ ylim(c(0,65)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_CCR$Species_colour, 
              shape=Tissue_Class_NLRs_logFC_BLAST_CCR$Like_shape, 
              size=log(Tissue_Class_NLRs_logFC_BLAST_CCR$averageRoot + 1)) + theme_bw() + 
  labs(x = "CCR-NLR homologue", y = "Expression in root tissue (TPM)")
#geom_hline(#yintercept=c(-1,0,1), 
#linetype=c("dashed","solid","dashed"), 
#color = "black", linewidth=0.5)
ggsave("Known_NLRs_averageRoot_per_Species_jitter_CCR_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_averageRoot_per_Species_jitter_CCR_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)
##### CCG10-NLRs #####

ggplot(Tissue_Class_NLRs_logFC_BLAST_CCG10, aes(Best_hit, log2FoldChange#, color=Species
)) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_CCG10$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_CCG10$Like_shape, 
              size=Tissue_Class_NLRs_logFC_BLAST_CCG10$dot_size) + theme_bw() + 
  geom_hline(yintercept=c(-1,0,1), linetype=c("dashed","solid","dashed"), 
             color = "black", linewidth=0.5) + ylim(c(-10, 10)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  labs(x = "CCG10-NLR homologue", y = "log(Fold Change) in leaf and root tissue")
ggsave("Known_NLRs_logFC_per_Species_jitter_CCG10_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_logFC_per_Species_jitter_CCG10_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)

ggplot(Tissue_Class_NLRs_logFC_BLAST_CCG10, aes(Best_hit, averageLeaf#, color=Species
)) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+ ylim(c(0,65)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_CCG10$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_CCG10$Like_shape, 
              size=log(Tissue_Class_NLRs_logFC_BLAST_CCG10$averageLeaf + 1)) + theme_bw() + 
  labs(x = "CCG10-NLR homologue", y = "Expression in leaf tissue (TPM)")
#geom_hline(#yintercept=c(-1,0,1), 
#linetype=c("dashed","solid","dashed"), 
#color = "black", linewidth=0.5)
ggsave("Known_NLRs_averageLeaf_per_Species_jitter_CCG10_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_averageLeaf_per_Species_jitter_CCG10_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)

ggplot(Tissue_Class_NLRs_logFC_BLAST_CCG10, aes(Best_hit, averageRoot#, color=Species
)) +
  geom_boxplot(#colour=c("#FFBB00", "#ff7400", "#3f0051", "#00b67e", "#488200",
    #        "#FF1100", "#003b84")
  )+ ylim(c(0,65)) +
  scale_color_manual(name='Species',
                     breaks=c('Phaseolus vulgaris', 'Glycine max', 'Cicer arietinum',
                              'Medicago truncatula',
                              'Lens culinaris',  
                              'Pisum sativum', 'Lathyrus sativus'),
                     values=c('Phaseolus vulgaris'= "#FF1100", 'Glycine max'="#ff7400", 
                              'Cicer arietinum'="#FFBB00",
                              'Medicago truncatula'="#488200",
                              'Lens culinaris'="#00b67e",  
                              'Pisum sativum'="#003b84", 'Lathyrus sativus'="#3f0051")) +
  # colour by species, shape by like, size by basemean
  geom_jitter(colour=Tissue_Class_NLRs_logFC_BLAST_CCG10$Species_colour, 
              #shape=Tissue_Class_NLRs_logFC_BLAST_CCG10$Like_shape, 
              size=log(Tissue_Class_NLRs_logFC_BLAST_CCG10$averageRoot + 1)) + theme_bw() + 
  labs(x = "CCG10-NLR homologue", y = "Expression in root tissue (TPM)")
#geom_hline(#yintercept=c(-1,0,1), 
#linetype=c("dashed","solid","dashed"), 
#color = "black", linewidth=0.5)
ggsave("Known_NLRs_averageRoot_per_Species_jitter_CCG10_3_no_like.png", width=1250, height=720, units='px', dpi = 125)
ggsave("Known_NLRs_averageRoot_per_Species_jitter_CCG10_3_no_like.svg", width=1250, height=720, units='px', dpi = 125)
# Idea
##### make stacked bar plot with frequencies per NLR homologue and with colors per species ####

Tissue_Class_NLRs_logFC_BLAST<- read.csv("Tissue_Class_NLRs_logFC_BLAST.csv", row.names=1)

Tissue_Class_NLRs_logFC_BLAST_Number_1 <- Tissue_Class_NLRs_logFC_BLAST %>%
  group_by(Species, Class) %>%
  count(Best_hit)

Tissue_Class_NLRs_logFC_BLAST_Number_1 <- Tissue_Class_NLRs_logFC_BLAST_Number %>%
  filter(Class %ni% c("TNP", "Other")) %>%
  group_by(Species) %>%
  mutate(ClassPercentage = n/(sum(n))*100) %>%#Percentage of NLRs with that class
  mutate(Total = sum(n)) %>%
  filter(Class %ni% c("Other"))

write.csv(Tissue_Class_NLRs_logFC_BLAST_Number_1, "Tissue_Class_NLRs_logFC_BLAST_Number_1.csv")

ggplot() +
  geom_bar(data=Tissue_Class_NLRs_logFC_BLAST_Number_1, aes(y = n, x = Class, fill = Species)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  #facet_grid( ~ Species) + 
  theme_bw() + 
  #scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  scale_fill_manual(values=c("#FF1100","#ff7400","#FFBB00", "#488200", "#00b67e", 
                              "#003b84","#3f0051"))
ggsave("Known_NLRs_logFC_per_Species_bar_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_logFC_per_Species_bar_2.svg", width = 1250, height = 720, units="px", dpi = 120)

##### make stacked bar plot with frequencies per NLR homologue and with colors per species ####

Tissue_NLRs_logFC_BLAST_Number <- Tissue_Class_NLRs_logFC_BLAST %>%
  group_by(Species, Tissue_logFC) %>%
  count(Best_hit)

Tissue_NLRs_logFC_BLAST_Number <- Tissue_NLRs_logFC_BLAST_Number %>%
  #filter(Class %ni% c("TNP", "Other")) %>%
  group_by(Species) %>%
  mutate(Tissue_logFCPercentage = n/(sum(n))*100) %>%#Percentage of NLRs with that Tissue_logFC
  mutate(Total = sum(n)) #%>%
  #filter(Tissue_logFC %ni% c("Other"))
  
Tissue_NLRs_logFC_BLAST_Number <- Tissue_NLRs_logFC_BLAST_Number %>%
  group_by(Tissue_logFC) %>%
  mutate(Total_Tissue = sum(n))

Tissue_NLRs_logFC_BLAST_Number <- Tissue_NLRs_logFC_BLAST_Number %>%
  group_by(Species, Tissue_logFC) %>%
  mutate(Total_Tissue_Species=sum(n))

Tissue_NLRs_logFC_BLAST_Number$Tissue_logFC <- factor(Tissue_NLRs_logFC_BLAST_Number$Tissue_logFC,
                                                levels=c("Not Expressed", "Leaf-specific",
                                                         "Leaf tendency", "Same expression", 
                                                         "Root tendency", "Root-specific"))
write.csv(Tissue_NLRs_logFC_BLAST_Number, "Tissue_NLRs_logFC_BLAST_Number.csv")

ggplot() +
  geom_bar(data=Tissue_NLRs_logFC_BLAST_Number, aes(y = n, x = Tissue_logFC, fill = Species)
           , stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  #facet_grid( ~ Species) + 
  theme_bw() + 
  #scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  scale_fill_manual(values=c("#FF1100","#ff7400","#FFBB00", "#488200", "#00b67e", 
                             "#003b84","#3f0051"))
ggsave("Known_NLRs_Tissue_logFC_per_Species_bar_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Tissue_logFC_per_Species_bar_2.svg", width = 1250, height = 720, units="px", dpi = 120)

##### make stacked bar plot with frequencies per NLR homologue and with colors per species ####

Tissue_NLRs_logFC_Class_BLAST_Number <- Tissue_Class_NLRs_logFC_BLAST %>%
  group_by(Species, Class, Tissue_logFC) %>%
  count(Best_hit)

Tissue_NLRs_logFC_Class_BLAST_Number <- Tissue_NLRs_logFC_Class_BLAST_Number %>%
  group_by(Species, Class) %>%
  mutate(Tissue_logFCPercentage = n/(sum(n))*100) %>%#Percentage of NLRs with that Tissue_logFC
  mutate(Total = sum(n)) #%>%
#filter(Tissue_logFC %ni% c("Other"))

Tissue_NLRs_logFC_Class_BLAST_Number$Tissue_logFC <- factor(Tissue_NLRs_logFC_Class_BLAST_Number$Tissue_logFC,
                                                      levels=c("Not Expressed", "Leaf-specific",
                                                               "Leaf tendency", "Same expression", 
                                                               "Root tendency", "Root-specific"))
Tissue_NLRs_logFC_Class_BLAST_Number$Class <- factor(Tissue_NLRs_logFC_Class_BLAST_Number$Class,
                                                            levels=c("CC-NLR", "TIR-NLR",
                                                                     "CCR-NLR", "CCG10-NLR", 
                                                                     "TNP", "Other"))
write.csv(Tissue_NLRs_logFC_Class_BLAST_Number, "Tissue_NLRs_logFC_Class_BLAST_Number.csv")

ggplot() +
  geom_bar(data=Tissue_NLRs_logFC_Class_BLAST_Number, aes(y = n, x = Tissue_logFC, fill = Class),
           stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Species) + 
  theme_bw() + 
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  #scale_fill_manual(values=c("#FF1100","#ff7400","#FFBB00", "#488200", "#00b67e", 
  #                           "#003b84","#3f0051"))
ggsave("Known_NLRs_Tissue_logFC_per_Species_Class_bar_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Tissue_logFC_per_Species_Class_bar_2.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot() +
  geom_bar(data=Tissue_NLRs_logFC_Class_BLAST_Number, aes(y = n, x = Tissue_logFC, fill = Species),
           stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Class) + 
  theme_bw() + 
  #scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  scale_fill_manual(values=c("#FF1100","#ff7400","#FFBB00", "#488200", "#00b67e", 
                             "#003b84","#3f0051"))
ggsave("Known_NLRs_Tissue_logFC_per_Class_Species_bar_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Tissue_logFC_per_Class_Species_bar_2.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot() +
  geom_bar(data=Tissue_NLRs_logFC_Class_BLAST_Number, aes(y = n, x = Species, fill = Class),
           stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Tissue_logFC) + 
  theme_bw() + 
  scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  #scale_fill_manual(values=c("#FF1100","#ff7400","#FFBB00", "#488200", "#00b67e", 
  #                           "#003b84","#3f0051"))
ggsave("Known_NLRs_Species_per_Tissue_logFC_Class_bar_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Species_per_Tissue_logFC_Class_bar_2.svg", width = 1250, height = 720, units="px", dpi = 120)

#### Bar graph with NLR homologues x-axis and frequency of Tissue_logFC ####

# remove stuff that is not from that class
Tissue_Class_NLRs_logFC_BLAST_Number_CC <- Tissue_NLRs_logFC_Class_BLAST_Number %>%
  filter(Class=="CC-NLR") %>%
  filter(Best_hit %ni% c("RUN1", "At5g47260")) %>%
  mutate(Best_hit = case_when(
    grepl("RPP13-like", Best_hit) ~ "RPP13",
    TRUE ~ Best_hit
  ))
Tissue_Class_NLRs_logFC_BLAST_Number_TIR <- Tissue_NLRs_logFC_Class_BLAST_Number %>%
  filter(Class=="TIR-NLR") %>%
  filter(Best_hit %ni% c("protein"))
Tissue_Class_NLRs_logFC_BLAST_Number_CCR <- Tissue_NLRs_logFC_Class_BLAST_Number %>%
  filter(Class=="CCR-NLR")
Tissue_Class_NLRs_logFC_BLAST_Number_CCG10 <- Tissue_NLRs_logFC_Class_BLAST_Number %>%
  filter(Class=="CCG10-NLR") %>%
  filter(Best_hit %ni% c("RGA1","RGA2", "RGA3", "TMV resistance protein N", 
                         "At5g05400","At5g47260", "RF9"))

Tissue_Class_NLRs_logFC_BLAST_Number_all <- rbind(Tissue_Class_NLRs_logFC_BLAST_Number_CC,Tissue_Class_NLRs_logFC_BLAST_Number_TIR,
      Tissue_Class_NLRs_logFC_BLAST_Number_CCR,Tissue_Class_NLRs_logFC_BLAST_Number_CCG10)

ggplot() +
  geom_bar(data=Tissue_Class_NLRs_logFC_BLAST_Number_all, aes(y = n, x = Best_hit, fill = Tissue_logFC),
           stat="identity",
           position='stack') +  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  facet_grid( ~ Class) + 
  theme_bw() + coord_flip() +
  #scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  scale_fill_manual(values=c("#000000","#488200","#b6c630", "#003b84", "#FFAA00", 
                             "#834c3b"))
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_all_NLR_class_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_all_NLR_class_2.svg", width = 1250, height = 720, units="px", dpi = 120)

# reorder all hits by class and expression
Tissue_Class_NLRs_logFC_BLAST_Number_all$Best_hit <- factor(Tissue_Class_NLRs_logFC_BLAST_Number_all$Best_hit,
                                                     levels=c("RFL1", "RPS2", "At1g61310", "At1g61300", "UNI", "At1g52660", "SUMM2", 
                                                              "At5g63020", "At4g19050", "At4g27190", "At4g27220",
                                                              "ADR1-L1", "ADR1-L3", "NRG1",
                                                              "RMG1", "SNC1", "RML1A", "Roq1", "DSC1", "RUN1", "RPV1", "TMV resistance protein N",
                                                              "RPP8","RF9", "R1B-8", "At1g50180", "RGA1", "RGA2", "RGA3", "RGA4", "LRRAC1", "RPM1", "RPP13"))
# order
#c("RFL1", "RPS2", "At1g61310", "At1g61300", "At1g61180", "At1g52660", "SUMM2", 
#"At5g63020", "At4g19050", "At4g27190", "At4g27220",
#"At4g33300", "At5g47280", "At5g66900",
#"At4g11170", "SNC1", "RML1A", "Roq1", "DSC1", "RUN1", "RPV1", "TMV resistance protein N",
#"RPP8","RF9", "R1B-8", "At1g50180", "RGA1", "RGA2", "RGA3", "RGA4", "At3g14460", "RPM1", "RPP13")
write.csv(Tissue_Class_NLRs_logFC_BLAST_Number_all, "Tissue_Class_NLRs_logFC_BLAST_Number_all.csv")
Tissue_Class_NLRs_logFC_BLAST_Number_all <- read.csv("Tissue_Class_NLRs_logFC_BLAST_Number_all.csv", row.names=1)

# testing which NLR homologues did not pass filtering
#test <- Tissue_NLRs_logFC_BLAST_Number %>%
#  filter(Best_hit %ni% unique(Tissue_Class_NLRs_logFC_BLAST_Number_all$Best_hit))

ggplot() +
  geom_bar(data=Tissue_Class_NLRs_logFC_BLAST_Number_all, aes(y = n, x = Best_hit, fill = Tissue_logFC),
           stat="identity",
           position='stack') +  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  #facet_grid( ~ Class) + 
  theme_bw() + coord_flip() +
  #scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  scale_fill_manual(values=c("#000000","#488200","#b6c630", "#003b84", "#FFAA00", 
                             "#834c3b"))
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_all_NLR_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_all_NLR_2.svg", width = 1250, height = 720, units="px", dpi = 120)


# x=Species must change to Best_hit
ggplot() +
  geom_bar(data=Tissue_Class_NLRs_logFC_BLAST_Number_CC, aes(y = n, x = Best_hit, fill = Tissue_logFC),
           stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  #facet_grid( ~ Class) + 
  theme_bw() + 
  #scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  scale_fill_manual(values=c("#000000","#488200","#b6c630", "#003b84", "#FFAA00", 
                             "#834c3b"))
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_CC_NLR_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_CC_NLR_2.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot() +
  geom_bar(data=Tissue_Class_NLRs_logFC_BLAST_Number_TIR, aes(y = n, x = Best_hit, fill = Tissue_logFC),
           stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  #facet_grid( ~ Class) + 
  theme_bw() + 
  #scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  scale_fill_manual(values=c("#000000","#488200","#b6c630", "#003b84", "#FFAA00", 
                             "#834c3b"))
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_TIR_NLR_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_TIR_NLR_2.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot() +
  geom_bar(data=Tissue_Class_NLRs_logFC_BLAST_Number_CCR, aes(y = n, x = Best_hit, fill = Tissue_logFC),
           stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  #facet_grid( ~ Class) + 
  theme_bw() + 
  #scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  scale_fill_manual(values=c("#000000","#488200","#b6c630", "#003b84", "#FFAA00", 
                             "#834c3b"))
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_CCR_NLR_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_CCR_NLR_2.svg", width = 1250, height = 720, units="px", dpi = 120)

ggplot() +
  geom_bar(data=Tissue_Class_NLRs_logFC_BLAST_Number_CCG10, aes(y = n, x = Best_hit, fill = Tissue_logFC),
           stat="identity",
           position='stack') +  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +# ylim(0,460) +
  #facet_grid( ~ Class) + 
  theme_bw() + 
  #scale_fill_manual(values=c("turquoise","salmon","thistle3","lightgoldenrod2","dodgerblue","green"))
  scale_fill_manual(values=c("#000000","#488200","#b6c630", "#003b84", "#FFAA00", 
                             "#834c3b"))
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_CCG10_NLR_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Known_NLRs_Tissue_logFC_per_NLR_Homologue_CCG10_NLR_2.svg", width = 1250, height = 720, units="px", dpi = 120)

#Tissue_logFC == "Not Expressed" ~ "#000000",
#Tissue_logFC == "Same expression" ~ "blue",
#Tissue_logFC == "Root-specific" ~ "#834c3b",
#Tissue_logFC == "Leaf-specific" ~ "#488200",
#Tissue_logFC == "Leaf tendency" ~ "#b6c630",
#Tissue_logFC == "Root tendency" ~ "#FFAA00",
Tissue_Class_NLRs_logFC_BLAST_Number_all$Tissue_logFC <- factor(Tissue_Class_NLRs_logFC_BLAST_Number_all$Tissue_logFC,
                                                            levels=c("Not Expressed", "Leaf-specific",
                                                                     "Leaf tendency", "Same expression", 
                                                                     "Root tendency", "Root-specific"))
Tissue_Class_NLRs_logFC_BLAST_Number_all$Best_hit <- factor(Tissue_Class_NLRs_logFC_BLAST_Number_all$Best_hit,
                                                            levels=c("RFL1", "RPS2", "At1g61310", "At1g61300", "UNI", "At1g52660", "SUMM2", 
                                                                     "At5g63020", "At4g19050", "At4g27190", "At4g27220",
                                                                     "ADR1-L1", "ADR1-L3", "NRG1",
                                                                     "RMG1", "SNC1", "RML1A", "Roq1", "DSC1", "RUN1", "RPV1", "TMV resistance protein N",
                                                                     "RPP8","RF9", "R1B-8", "At1g50180", "RGA1", "RGA2", "RGA3", "RGA4", "LRRAC1", "RPM1", "RPP13"))

ggplot(Tissue_Class_NLRs_logFC_BLAST_Number_all,   # Draw heatmap-like plot
       aes(Tissue_logFC, Best_hit, fill = log(n))) +  geom_tile() + 
  #scale_y_discrete(breaks = TotalClassNLRs_Number$n 
  # , labels = Conserved_NLRs$seqname) 
  #+ #geom_hline(yintercept = 0.5 + 0:25908, colour = "black", linewidth = 0.5) +
  # scale_z_log10() +
  scale_fill_gradient2(low="white", mid = "#4393c3",high="#67001f", na.value="white"#, trans = "log" 
  ) + 
  ggtitle("Expression of known NLR homologues", subtitle = "Leaves and Roots") +
  ylab("NLR homologue") + #geom_segment(data=my.lines, aes(x=1, y, xend=6.5, yend=y), linewidth=0.01, inherit.aes=F)
  theme_minimal() #+
  # this line allows to divide the graph in NLR homologue groups
  #geom_hline(yintercept = c((11),(8),(11),(3)), linewidth=0.00001)
ggsave("Heatmap_NLR_homologues_expression_pattern_tendency_2.png", width = 1250, height = 720, units="px", dpi = 120)
ggsave("Heatmap_NLR_homologues_expression_pattern_tendency_2.svg", width = 1250, height = 720, units="px", dpi = 120)


# Idea
# make a table to aid in the discussion stating the name of the NLR and where
# it has been studied