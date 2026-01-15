# Author: Rita Maravilha Marques
# Script to process NLR class variants predicted by NLRtracker per species

# Select working directory
setwd("C://Users/Rita/Dropbox (KamounLab)/Kamounity folder/Rita/NLRtracker/AnalysisNLRtracker131/")
setwd("D:/Dropbox (KamounLab)/Kamounity folder/Rita/NLRtracker/AnalysisNLRtracker131")
setwd("D:/KamounLab Dropbox/Kamounity folder/Rita/NLRtracker/AnalysisNLRtracker131")
setwd("C:/Users/RitaML/KamounLab Dropbox/Kamounity folder/Rita/NLRtracker/AnalysisNLRtracker131")

#Load in the required packages
library(readxl)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(stringr)

# 1. Align NB ARC files from the NLR tracker output with the RefPlantNLR
# datasets to correctly assign each NLR to its respective class.
# 2. Create the NLR Classes w RefPlantNLR file with the following columns:
# Species seqname Class
#generated an excel file labelled “NLR Classes with RefPlantNLR.xlsx” (Table S3) 
#by pasting the sequence names of the NLRs from each clade and creating a 
#column “Class” with the respective assigned class
df1 <- read_xlsx("NLR Classes Legumes w RefPlantNLR 032024.xlsx")
df1 <- df1 %>% 
  mutate(Species = case_when(
    Species == "Lathyrus sativus" ~ "Lathyrus sativus Rbp",
    Species == "Lathyrus sativus 2024" ~ "Lathyrus sativus",
    Species == "Pisum sativum" ~ "Pisum sativum Caméor",
    Species == "Pisum sativum ZW6" ~ "Pisum sativum",
    Species == "Cajanus cajan" ~ "Cajanus cajan contaminated",
    Species == "Cajanus cajan v1.1" ~ "Cajanus cajan",
    Species == "Lotus japonicus" ~ "Lotus japonicus v1.0",
    Species == "Lotus japonicus Gifu" ~ "Lotus japonicus",
    TRUE ~ Species
  ) ) 

RefPlantNLRdelete <- read_xlsx("RefPlantNLR.xlsx") #to remove RefPlantNLR NLR references
LOCUS <- read_tsv("Fabales_LOCUS_06032024.tsv") #for deduplicating

#listing the species
LegumeSpecies <- c("Abrus precatorius","Aeschynomene evenia",
                   "Arabidopsis thaliana",
                   "Arachis duranensis",
                   "Arachis hypogaea","Arachis ipaensis",
                   "Cajanus cajan",
                   "Cercis canadensis", 
                   "Cicer arietinum",
                   "Glycine max","Glycine soja",
                   "Lens culinaris","Lupinus albus",
                   "Lupinus angustifolius","Lotus japonicus","Lathyrus sativus",
                   "Medicago truncatula", 
                   "Prosopis alba","Prosopis cineraria",  
                   "Pisum sativum",
                   "Phaseolus acutifolius", "Phaseolus lunatus", "Phaseolus vulgaris",
                   "Oryza sativa","Solanum lycopersicum",
                   "Trifolium pratense",
                   "Vicia faba", "Vicia sativa", 
                   "Vigna angularis",
                   "Vigna radiata","Vigna unguiculata", "Zea mays")

#ordering the legumes according to phylogeny
#Source https://link.springer.com/chapter/10.1007/13836_2021_89
LegumePhylogeny <- c(
  "Oryza sativa","Zea mays",
  "Arabidopsis thaliana","Solanum lycopersicum", 
  "Cercis canadensis",
  "Prosopis alba", 
  "Prosopis cineraria", 
  "Lupinus angustifolius","Lupinus albus",
  "Aeschynomene evenia","Arachis duranensis","Arachis ipaensis","Arachis hypogaea",
  "Cajanus cajan","Abrus precatorius","Vigna radiata",
  "Vigna angularis","Vigna unguiculata", "Phaseolus acutifolius",
  "Phaseolus vulgaris",
  "Phaseolus lunatus",
  "Glycine soja",
  "Glycine max", "Lotus japonicus","Cicer arietinum", "Medicago truncatula",
  "Trifolium pratense","Lens culinaris", 
  "Vicia sativa", "Vicia faba", 
  "Pisum sativum", 
  "Lathyrus sativus")

# filter out every row whose
# Seqname matches with RefPlantNLR's seqname
df2 <- df1 %>%
  filter(!seqname %in% RefPlantNLRdelete$seqname) 
# https://www.datasciencemadesimple.com/delete-or-drop-rows-in-r-with-conditions-2/
# Dropping rows in R
df2<-df2[!(df2$seqname=="AfsR_STRGR 2"),] #For some reason this ID is not in
df2<-df2[!(df2$seqname=="PH0952 2"),] #For some reason this ID is not in 
#RefPlantNLR
write_csv(df2,"NLR Classes Legumes Filtered 06032024.csv")
df2 <- read_csv("NLR Classes Legumes Filtered 06032024.csv")

# Count how many NLRs are in each class, for each Species
NLRClasses <- df2 %>%
  group_by(Species) %>%
  count(Class) 
write_csv(NLRClasses,"NLR Number per Class Not deduplicated 06032024.csv")

##### Generating graphs from non-deduplicated data ##################
NLRClasses <- NLRClasses %>%
  group_by(Species) %>%
  mutate(ClassPercentage = n/(sum(n))*100) #Percentage of NLRs with that class
View(NLRClasses)
write_csv(NLRClasses, "NLR Classes Not Deduplicated 06032024.csv")

NLRClasses %>% #Total number of NLRs
  ggplot(aes(x = Species, y = n, fill = Class)) +
  geom_bar(position="stack",stat="identity") + coord_flip() + theme_classic(base_size = 8) +
  scale_fill_manual(values=c("turquoise","lightgoldenrod2","thistle3","palegreen3","salmon","dodgerblue")) +
  labs(title="Classes of Legume NLRs (not Deduplicated)", y="NLR Number", x="Species", fill="Classes") +
  scale_x_discrete(limits = LegumePhylogeny)
#colors for ggplot http://sape.inf.usi.ch/quick-reference/ggplot2/colour

NLRClasses %>% #Percentages
  ggplot(aes(x = Species, y = ClassPercentage, fill = Class)) +
  geom_bar(position="stack",stat="identity") + coord_flip() +
  theme_classic() +
  scale_fill_manual(values=c("turquoise","lightgoldenrod2","thistle3","palegreen3","salmon","dodgerblue")) +
  labs(title="Percentages of NLR Classes in Legumes (not Deduplicated)", y="NLR Class Percentage", x="Species", fill="Classes") +
  scale_x_discrete(limits = LegumePhylogeny)

#arrange the table according to phylogeny and then NLR class
NLRClasses <- NLRClasses %>%
  arrange(factor(Species,levels=LegumePhylogeny),Class) 

# Filtering the main table according to NLR class
NLRCCNumber <- NLRClasses %>%
  filter(Class=="CC-NLR") %>%
  filter(Species!="Zea mays" & Species!="Oryza sativa")
NLRTIRNumber <- NLRClasses %>%
  filter(Class=="TIR-NLR")
NLRCCG10Number <- NLRClasses %>%
  filter(Class=="CCG10-NLR") %>%
  filter(Species!="Zea mays" & Species!="Oryza sativa")
NLRCCRNumber <- NLRClasses %>%
  filter(Class=="CCR-NLR") %>%
  filter(Species!="Zea mays" & Species!="Oryza sativa")
NLRTNPNumber <- NLRClasses %>%
  filter(Class=="TNP") %>%
  filter(Species!="Zea mays" & Species!="Oryza sativa")

# Fetching each number of NLRs in each class
TIR <- NLRTIRNumber$n
CC <- NLRCCNumber$n
CCG10 <- NLRCCG10Number$n
CCR <- NLRCCRNumber$n
TNP <- NLRTNPNumber$n

# Make a data frame with legume classes as columns
NLRs <- data.frame(Species = LegumePhylogeny[3:length(LegumePhylogeny)],TIR,CC,CCG10,CCR,TNP)

# Making ratios between NLR class numbers
NLRs$TIRCC <- NLRTIRNumber$n/NLRCCNumber$n
NLRs$CCG10CCR <- NLRCCG10Number$n/NLRCCRNumber$n
NLRs$TIRCCR <- NLRTIRNumber$n/NLRCCRNumber$n

NLRTIRCC <- data.frame(Species=c(
  #"Oryza sativa","Zea mays",
  "Arabidopsis thaliana","Solanum lycopersicum", 
  "Cercis canadensis",
  "Prosopis alba", 
  "Prosopis cineraria", 
  "Lupinus angustifolius","Lupinus albus",
  "Aeschynomene evenia","Arachis duranensis","Arachis ipaensis","Arachis hypogaea",
  "Cajanus cajan","Abrus precatorius","Vigna radiata",
  "Vigna angularis","Vigna unguiculata", "Phaseolus acutifolius",
  "Phaseolus vulgaris",
  "Phaseolus lunatus",
  "Glycine soja",
  "Glycine max", "Lotus japonicus","Cicer arietinum", "Medicago truncatula",
  "Trifolium pratense","Lens culinaris", 
  "Vicia sativa", "Vicia faba", 
  "Pisum sativum",
  "Pisum sativum ZW6", 
  "Lathyrus sativus"), TIR, CC, NLRs$TIRCC)
NLRTIRCC <- NLRTIRCC %>% 
  rename("TIRCC" = "NLRs.TIRCC")

NLRTIRCC$Ratio <- 1- NLRTIRCC$TIRCC
NLRTIRCC <- NLRTIRCC %>%
  mutate(cond = case_when(
    Ratio<0 ~ 'salmon',
    Ratio>0 ~ 'turquoise',
    TRUE ~ 'yellow'   #anything that does not meet the criteria above
  ))
write.csv(NLRTIRCC, "NLRTIRCC not Deduplicated.csv")
NLRTIRCC <- read.csv("NLRTIRCC not Deduplicated.csv")

annotation <- data.frame(
  x = c(1:31),
  y = c(rep(-3.5,31),rep(2,31)),
  label = c(NLRTIRCC$TIR, NLRTIRCC$CC),
  color = c(rep("salmon",31),rep("turquoise",31))
)

#Ratio between TIR-NLRs and CC-NLRs
NLRTIRCC %>% #Total number of NLRs
  ggplot(aes(x = Species, y = Ratio)) +
  geom_bar(position="stack",stat="identity", fill=NLRTIRCC$cond) + coord_flip() + geom_hline(yintercept=c(-1,0,1), linetype="dashed",
                                                                                             color = "black", size=0.5) +
  scale_fill_identity() + theme_classic() + geom_text(data=annotation, aes(x=x, y=y, label=label),
                                                      color=annotation$color,
                                                      size=3,
                                                      angle=0, fontface="bold") +
  labs(title="Ratio between TIR-NLRs and CC-NLRs in Legumes", y="1 - TIR-NLR/CC-NLR Ratio", x="Species") +
  scale_y_continuous(limits=c(-3.5,2.5)) +
  scale_x_discrete(limits = LegumePhylogeny[3:33])

##### Deduplicating alternate splicing variants #####################

# Removing the Splicing alternatives to get the Loci 
# (i.e. remove every character after the last dot)
#df2 <- df2 %>%
#  select(df2$Species=="Glycine max"|df2$Species=="Lens culinaris") %>%
#  mutate(df2$seqname, gsub("[.][^.]+$", "", df2$seqname))

# Join the NLR data with the locus data
#df2 <- df2 %>%
#  filter(Species!="Pisum sativum ZW6" & Species!="Cercis canadensis" & Species!="Vicia faba"
#         & Species!="Vicia sativa" & Species!="Prosopis cineraria")
df2 <- read_csv("NLR Classes Legumes Filtered 06032024.csv")
# Group species with the same transcript behavior
# Trifolium and Cercis do not have .1 but they have random '.' in the sequence
Species_list_1 <- LegumeSpecies[! LegumeSpecies %in% c('Trifolium pratense', 'Cercis canadensis')]
# add the previous lathyrus assembly because it behaves similarly
Species_list_1 <- c(Species_list_1, "Lathyrus sativus Rbp", "Lotus japonicus v1.0")
# all these species have _1 and _(1) transcripts
Species_list_2 <- LegumeSpecies[LegumeSpecies %in% c("Lotus japonicus","Prosopis cineraria", "Zea mays", "Solanum lycopersicum",
                                                     "Oryza sativa", "Arabidopsis thaliana", "Vigna unguiculata",
                                                    "Prosopis alba", "Medicago truncatula", "Glycine soja", 
                                                    "Cajanus cajan", "Arachis ipaensis", "Arachis hypogaea",
                                                    "Arachis duranensis", "Abrus precatorius", "Pisum sativum")]
Species_list_2 <- c(Species_list_2, "Cajanus cajan contaminated")
# has both . and _ in the sequence
Species_list_3 <- c("Lupinus albus", "Lathyrus sativus Rbp", "Cicer arietinum")

# Make the seqname ready to be a single locus. Keep the original transcript info
df2 <- df2 %>%
  mutate(original = seqname) %>%
  mutate(seqname = case_when(
    # remove everything after the last .
    Species %in% Species_list_1 ~ gsub("[.][^.]+$", "", seqname),
    TRUE ~ seqname
   )) %>%
  mutate(seqname = case_when(
    # prosopis cineraria & Co. have an extra _
    Species %in% Species_list_2 ~ sub('_', '.', seqname),
    Species %in% Species_list_3 ~ sub('_', '-', seqname),
    TRUE ~ seqname # sub only replaces the first occurrence
   )) %>%
  mutate(seqname = case_when(
    Species == "Lathyrus sativus Rbp" ~ sub('_', '-', seqname),
    TRUE ~ seqname # Rbp has 2 flipping _
  )) %>%
  mutate(seqname = gsub("[_][^_]+$", "", seqname)) %>% 
  mutate(seqname = case_when( 
    # Phaseolus lunatus has an extra .v1 after the transcript
    Species == "Phaseolus lunatus" ~ gsub("[.][^.]+$", "", seqname),
    TRUE ~ seqname
  )) %>% 
  # remove everything after the last _
  mutate(seqname = case_when(# put back the _ in prosopis cineraria & Co.
    Species %in% Species_list_2 ~ sub('\\.', '_', seqname), # sub only replaces the first occurrence   
    # the . must be escaped, or it just means "any character"
    Species %in% Species_list_3 ~ sub("-", "_", seqname),
    Species == "Lathyrus sativus Rbp" ~ sub("-", "_", seqname),
    TRUE ~ seqname)) 

#LOCUS <- LOCUS %>%
#  mutate(original = seqname) %>%
#  mutate(seqname = gsub("[.][^.]+$", "", seqname))

LOCUS <- read_tsv("Fabales_LOCUS_06032024.tsv") #for deduplicating

# ensure LOCUS always has a Locus
LOCUS$Locus <- ifelse(is.na(LOCUS$Locus), LOCUS$seqname, LOCUS$Locus)

# LOCUS has 2x faba and phaseolus vulgaris entries?
# remove phaseolus vulgaris locus = 5 characters
# remove faba locus = 
# problem: faba, Vicia sativa, phaseolus, glycine max
LOCUS <- LOCUS %>%
  # copy the info in coded_by when the phaseolus vulgaris or cercis canadensis
  #locus is incomplete
  mutate(Locus = case_when(
    Locus == "Phvul" ~ str_sub(LOCUS$coded_by, end = -3),
    Locus == "cerca.ISC453364.gnm1.ann1" ~ seqname,
    # in faba, the Locus must be standardized
    Organism == "Vicia faba" & nchar(Locus) %in% c(nchar("Vfaba.Hedin2.R1.5g105680.1"), nchar("Vfaba.Hedin2.R1.Ung118840.1")) ~ str_sub(LOCUS$Locus, end = -3), 
    Organism == "Vicia sativa" & grepl(".t", Locus, fixed = TRUE) #find .t in locus 
    ~ str_sub(LOCUS$Locus, end = -4), # remove the transcript part
    Organism == "Glycine max" & grepl("Glyma.", Locus, fixed = TRUE)
    ~ str_sub(LOCUS$Locus, end = -3), # remove the transcript part
    TRUE ~ Locus
  )) %>%
  mutate(coded_by = case_when(
    coded_by == "cerca.ISC453364.gnm1.ann1" ~ seqname,
    TRUE ~ coded_by
  )  )
LOCUS <- LOCUS %>%
  group_by(Genome, seqname) %>% # we don't group by locus because we want to keep them
  arrange(desc(seqname)) %>%
  arrange(desc(Length)) %>% # keep the larger
  slice(1)
#2363945 to 2204250 lines -> maybe it got rid of some weird cases like cercis

write_tsv(LOCUS, "Fabales_LOCUS_06032024_edited.tsv")
LOCUS <- read_tsv("Fabales_LOCUS_06032024_edited.tsv")
  
df2 <- df2 %>%
  mutate(original=case_when(
    # find _() in the original seqname. These will not have a match in LOCUS,
    # we need to correct for that before joining the NLRs with LOCUS information
    grepl("_(", original, fixed = TRUE) ~ str_sub(df2$original, end = -5),
    TRUE ~ original
  ))
write_csv(df2, "NLR Classes Legumes Filtered 06032024 edited.csv")
df2 <- read_csv("NLR Classes Legumes Filtered 06032024 edited.csv")
df2_1 <- left_join(df2, LOCUS, by=c("original" = "seqname"))
df2_1 <- df2_1 %>% # 15075
  group_by(Genome, seqname) %>%
  arrange(desc(seqname)) %>%
  arrange(desc(Length)) %>%
  arrange(desc(Locus)) %>%
  slice(1)
df2_2 <- left_join(df2, LOCUS, by=c("seqname" = "Locus" )) # 19418
# got rid of some repetitive matches
df3 <- full_join(df2_1, df2_2)

df3 <- df3 %>% 
  # if seqname matches, left join with seqname, if it does not match, left join 
  #with locus
  # this way, all rows must have a Locus given by LOCUS. 
  # Remove the rows that don't because that means they were not corresponded on both dataframes.
  drop_na(Locus) %>%
  mutate(coded_by = case_when(
    is.na(coded_by) ~ Locus,
    TRUE ~ coded_by # fill in the coded_by column
  )) %>%
  mutate(Length = case_when(
    is.na(Length) ~ nchar(sequence),
    TRUE ~ Length # fill in the length column
    )) # 29877 to 15329

write_csv(df3,"df3 06032024.csv")
df3 <- read_csv("df3 06032024.csv")

# Order by seqname and Length, keep the most complete transcript (maximum length)

NLR_LOCUS <- df3 %>%
  group_by(Genome, Locus) %>%
  arrange(desc(seqname)) %>%
  arrange(desc(Length)) %>%
  slice(1)
#11806
View(NLR_LOCUS)
write_csv(NLR_LOCUS,"Locus and Type Legume NLRs 06032024.csv")
NLR_LOCUS <- read_csv("Locus and Type Legume NLRs 06032024.csv")


##### Get how many deduplicated NLRs are in each type ######################################
NLRClassesDeduplicatedNumber <- NLR_LOCUS %>%
  group_by(Species) %>%
  count(Class)
# Get the percentages of NLRs in each class
NLRClassesDeduplicatedNumber <- NLRClassesDeduplicatedNumber %>%
  group_by(Species) %>%
  mutate(ClassPercentage = n/(sum(n))*100) #Percentage of NLRs with that class
View(NLRClassesDeduplicatedNumber)
write_csv(NLRClassesDeduplicatedNumber, "NLRClassesDeduplicatedNumber_06032024.csv")

NLRClassesDeduplicatedNumber %>% #Total number of NLRs
  ggplot(aes(x = Species, y = n, fill = Class)) +
  geom_bar(position="stack",stat="identity") + coord_flip() +
  theme_classic(base_size = 7) +
  scale_fill_manual(values=c("turquoise","lightgoldenrod2","thistle3","palegreen3","salmon","dodgerblue")) +
  labs(title="Classes of Legume NLRs", y="NLR Number", x="Species", fill="Classes") +
  scale_x_discrete(limits = LegumePhylogeny) #+ 
  #theme(axis.text=element_text(size=7))
#colors for ggplot http://sape.inf.usi.ch/quick-reference/ggplot2/colour
ggsave("Output/Classes of Legume NLRs deduplicated 06032024.png", width = 1920, 
       height = 1080, units = "px")


NLRClassesDeduplicatedNumber %>% #Percentages
  ggplot(aes(x = Species, y = ClassPercentage, fill = Class)) +
  geom_bar(position="stack",stat="identity") + coord_flip() +
  theme_classic(base_size = 7) +
  scale_fill_manual(values=c("turquoise","lightgoldenrod2","thistle3","palegreen3","salmon","dodgerblue")) +
  labs(title="Percentages of NLR Classes in Legumes", y="NLR Class Percentage", x="Species", fill="Classes") +
  scale_x_discrete(limits = LegumePhylogeny)
ggsave("Output/Class Percentages of Legume NLRs deduplicated 06032024.png", width = 1920, 
       height = 1080, units = "px")

##Normalize the NLR Number with the number of different loci in Fabales_LOCUS ##
NLRClassesDeduplicatedNumber <- read_csv("NLRClassesDeduplicatedNumber_06032024.csv")
LOCUS <- read_tsv("Fabales_LOCUS_06032024_edited.tsv")

LOCUS <- LOCUS %>%
  #add CAAS_Psat_ZW6_1.0 to the Genome for Pisum sativum ZW6 
  mutate(Genome = case_when(
    Organism == "Pisum sativum ZW6" ~ "CAAS_Psat_ZW6_1.0",
    TRUE ~ Genome) )

#count the amount of loci each genome has
LOCI_Number <- LOCUS %>%
  group_by(Genome) %>%
  summarize(count = n_distinct(Locus))

LOCI_Genome <- LOCI_Number$Genome #get the genome IDs
View(LOCI_Genome)
# Genome IDs
#"Abrus_2018", "Aeschynomene_evenia_v1.0.protein", "Arabidopsis thaliana TAIR 10.1", "Aradu1.1",
#"arahy.Tifrunner.gnm1.KYV3", "Araip1.1", "ASM33114v1", "ASM411807v1", "ASM419377v2", "ASM479914v1", "C.cajan_V1.0",
#"C.cajan_V1.1", "CAAS_Psat_ZW6_1.0", "cerca.ISC453364.gnm1.ann1.B05Z",
#"cicar.CDCFrontier.gnm1.ann1.nRhs.protein", "Glycine max Williams 82", "Glycine_max_v4.0",
#"LATSA3860_EIv1.0.annotation.gff3.pep", "Lens culinaris Redberry v2.0", "LjGifu_v1.2", "lotja.MG20.gnm3.ann1.WF9B.protein",
#"lupal.Amiga.gnm1.ann1.3GKS.protein", "lupan.Tanjil.gnm1.ann1.nnV9.protein", "LupAngTanjil_v1.0", "MtrunA17r5.0-ANR",
#"Oryza sativa Nipponbare IRGSP 1.0", "Pacutifolius_580_v1.0.protein", "PacutifoliusWLD_581_v2.0.protein", "PC_final"
#"phavu.G19833.gnm2.ann1.PB8d.protein", "PhaVulg1_0", "pissa.Cameor.gnm1.ann1.7SZR.protein", "Pisum_sativum_v1a_prot", 
#"Plunatus_563_V1.protein", "Solanum lycopersicum SL3.0", "tripr.MilvusB.gnm2.ann1.DFgp.protein",
#"Vicia faba HEDIN_TMP5", "Vicia sativa",
#"vigan.Gyeongwon.gnm3.ann1.3Nz5.protein", "Vigan1.1", "vigra.VC1973A.gnm6.ann1.M1Qs.protein", "Vradiata_ver6",
#"Zea mays B73 v5"

# NEW Order - 032023
#"ASM33114v1", "ASM411807v1", "ASM419377v2", "ASM479914v1", 
#"Abrus_2018", "Aeschynomene_evenia_v1.0.protein", "Arabidopsis thaliana TAIR 10.1", "Aradu1.1",
#"Araip1.1", "C.cajan_V1.0", "CAAS_Psat_ZW6_1.0", "Glycine max Williams 82",
#"Glycine_max_v4.0", "LATSA3860_EIv1.0.annotation.gff3.pep", "Lens culinaris Redberry v2.0",
#"LjGifu_v1.2"
#"LupAngTanjil_v1.0", "MtrunA17r5.0-ANR", "Oryza sativa Nipponbare IRGSP 1.0", "PC_final",
#"PacutifoliusWLD_581_v2.0.protein", "Pacutifolius_580_v1.0.protein", "PhaVulg1_0",
#"Pisum_sativum_v1a_prot", "Plunatus_563_V1.protein", "Solanum lycopersicum SL3.0"
#"Vicia faba HEDIN_TMP5", "Vicia sativa", "Vigan1.1", "Vradiata_ver6", "Zea mays B73 v5",
#"arahy.Tifrunner.gnm1.KYV3", "cerca.ISC453364.gnm1.ann1.B05Z", "cicar.CDCFrontier.gnm1.ann1.nRhs.protein",
#"lotja.MG20.gnm3.ann1.WF9B.protein", "lupal.Amiga.gnm1.ann1.3GKS.protein", "lupan.Tanjil.gnm1.ann1.nnV9.protein",
#"phavu.G19833.gnm2.ann1.PB8d.protein", "pissa.Cameor.gnm1.ann1.7SZR.protein", "tripr.MilvusB.gnm2.ann1.DFgp.protein"
#"vigan.Gyeongwon.gnm3.ann1.3Nz5.protein", "vigra.VC1973A.gnm6.ann1.M1Qs.protein"

#LOCI_Organism <- c("Abrus precatorius", "Aeschynomene evenia", "Arabidopsis thaliana", "Arachis duranensis",
#"Arachis hypogaea", "Arachis ipaensis", "", "Vigna unguiculata", "Glycine soja", "Prosopis alba", "Cajanus cajan",
#"Pisum sativum ZW6", "Cercis canadensis", "Cicer arietinum", "Glycine max", "",
#"Lathyrus sativus", "Lens culinaris", "Lotus japonicus",
#"Lupinus albus", "Lupinus angustifolius", "", "Medicago truncatula",
#"Oryza sativa", "", "Phaseolus acutifolius", "Prosopis cineraria",
#"Phaseolus vulgaris", "", "", "Pisum sativum", "Phaseolus lunatus", "Solanum lycopersicum", "Trifolium pratense",
#                   "Vicia faba", "Vicia sativa","Vigna angularis", "", "Vigna radiata", "", "Zea mays")

LOCI_Organism <- c("", "Vigna unguiculata", "Glycine soja", "Prosopis alba", 
                   "Abrus precatorius", "Aeschynomene evenia", "Arabidopsis thaliana", "Arachis duranensis",
                   "Arachis ipaensis", "", "Cajanus cajan", "Pisum sativum", "Glycine max",
                   "", "", "Lens culinaris", "Lotus japonicus",  
                   "", "", "Medicago truncatula", "Oryza sativa", "Prosopis cineraria",
                   "Phaseolus acutifolius", "", "",
                   "", "Phaseolus lunatus", "Solanum lycopersicum",
                   "Vicia faba", "Vicia sativa", "", "", "Zea mays",
                   "Arachis hypogaea", "Cercis canadensis", "Cicer arietinum",
                   "Lathyrus sativus",
                   "", "Lupinus albus", "Lupinus angustifolius",
                   "Phaseolus vulgaris", "", "Trifolium pratense", 
                   "Vigna angularis", "Vigna radiata")

LOCI_Number <- LOCI_Number %>%
  mutate(Organism = LOCI_Organism) %>% #add a new column with the Organism name
  filter(Organism != "") %>% #get rid of the genomes which weren't used
  arrange(Organism) #arrange the table by Organism name

# Remove the other Pisum sativum genome
#LOCI_Number <- LOCI_Number[-c(23), ]

# the protein-coding loci number does not correspond to the one reported for

write.csv(LOCI_Number, "Loci Number for Each Species_06032024.csv")
LOCI_Number <- read.csv("Loci Number for Each Species_06032024.csv")
colnames(NLRClassesDeduplicatedNumber) <- c("Species","Class","NLRNumber","ClassPercentage") #rename the columns
#join the deduplicated numbers with the LOCI_Number
NLRClassesDeduplicatedNumber <- left_join(NLRClassesDeduplicatedNumber,LOCI_Number, by=c("Species" = "Organism"))

NLRClassesDeduplicatedNumber <- NLRClassesDeduplicatedNumber %>%
  group_by(Species) %>%
  mutate(NormalizedNumber = NLRNumber/count) #divide the NLR Numbers for the total protein coding loci of each genome
View(NLRClassesDeduplicatedNumber)

#arrange the table according to phylogeny and then NLR class
NLRClassesDeduplicatedNumber <- NLRClassesDeduplicatedNumber %>%
  arrange(factor(Species,levels=LegumePhylogeny),Class) 

#get the total amount of NLRs in each species
NLRClassesDeduplicatedNumber <- NLRClassesDeduplicatedNumber %>%
  group_by(Species) %>%
  mutate(TotalNLRNumber = sum(NLRNumber)) %>%
  drop_na(Genome) # remove old genomes

#get a list of the total number of NLRs
NLR_number <- unique(NLRClassesDeduplicatedNumber$TotalNLRNumber) 
#Phaseolus vulgaris + Glycine soja and Cercis + Araip 
#had the same total number of NLRs, so I had to manually repeat (rep) 328 twice 
#in the list
NLR_number <- append(NLR_number, 328, after = 20)
NLR_number <- append(NLR_number, 435, after = 11)

#NLR_number <- c(NLR_number[1:4], rep(NLR_number[5], 2), NLR_number[6:18], 
#                rep(NLR_number[19], 2), NLR_number[20:30]) 
write.csv(NLRClassesDeduplicatedNumber, "NLR Classes Normalized 06032024.csv") #save the table
NLRClassesDeduplicatedNumber <- read.csv("NLR Classes Normalized 06032024.csv")

NLRClassesDeduplicatedNumber %>% #Normalized NLRs
  ggplot(aes(x = Species, y = NormalizedNumber*100, fill = Class)) +
  #by multiplying NLR Number/Number of Loci by 100, we get the percentage of NLR genes!
  theme_classic(base_size = 7) +
  geom_bar(position="stack",stat="identity") + coord_flip() +
  scale_fill_manual(values=c("turquoise","lightgoldenrod2","thistle3","palegreen3","salmon","dodgerblue")) +
  labs(title="Normalized Legume NLRs", y="Percentage of Protein coding loci encoding for NLRs", x="Species", fill="Classes") +
  scale_x_discrete(limits = LegumePhylogeny) +
  annotate("text", label=paste(#"TN: ", NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="TN Other"],
    #", TIR: ", NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="TIR-NLR"],
    #", Other: ", NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="Other"],
    #", CCR: ", NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="CCR-NLR"],
    #", CCG10: ", NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="CCG10-NLR"],
    #", CC: ", NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="CC-NLR"],
    #",
    " ", NLR_number,
    #NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="TN Other"] +
    # NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="TIR-NLR"] +
    #NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="Other"] +
    #NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="CCR-NLR"] +
    #NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="CCG10-NLR"] +
    #NLRClassesDeduplicatedNumber$NLRNumber[NLRClassesDeduplicatedNumber$Class=="CC-NLR"],
    sep=""), fontface = "bold", #might be kind of spaghetti code but it works!
           #annotate the BUSCO summary on the graph: fetch the number of BUSCOs in each category for each species
           #https://ggplot2-book.org/annotations.html annotations in ggplot2
           y = 2.3, x = seq_along(LegumePhylogeny),
           size = 0.75*3*1,
           colour = "black",
           hjust=0,
           family="sans"
  )
ggsave("Output/Normalized Legume NLRs 04122024.png", width = 1920, 
       height = 1080, units = "px")
ggsave("Output/Normalized Legume NLRs 04122024.svg", width = 1920, 
       height = 1080, units = "px")
# colors for ggplot http://sape.inf.usi.ch/quick-reference/ggplot2/colour
# try to annotate the NLR Numbers in the graph

##### Boxplots for Normalized Number of NLRs #######################
#https://r-charts.com/distribution/add-points-boxplot/
NLRClassesDeduplicatedNumber <- read.csv("NLR Classes Normalized 06032024.csv")

# Create a file with info exclusively for legume species (exclude the outgroups)
LegumeNLRClasses <- subset(NLRClassesDeduplicatedNumber, Species!="Oryza sativa"
  & Species!="Arabidopsis thaliana" & Species!="Solanum lycopersicum" & Species!="Zea mays")

# Vertical box plot by group
LegumeNLRClasses$Class <- factor(LegumeNLRClasses$Class,
                                 levels=c("CC-NLR", "TIR-NLR",
                                          "CCG10-NLR", "CCR-NLR", "TNP", "Other"))
View(LegumeNLRClasses)
write.csv(LegumeNLRClasses, "Legume NLR Classes 06032024.csv")
LegumeNLRClasses <- read.csv("Legume NLR Classes 06032024.csv", row.names=1)

# Normalized Number boxplots
png(filename="Output/Percentage of NLR-coding loci in Legumes 04122024.png",
    width = 854, height = 480)
boxplot(LegumeNLRClasses$NormalizedNumber*100 ~ LegumeNLRClasses$Class,
        data = LegumeNLRClasses, col = "white", theme(size=11),
        ylab="Percentage of Protein-coding loci encoding for NLRs",
        xlab="NLR Class", main="Percentage of NLR-coding loci in Legumes", notch = FALSE, yaxt="n")
#You can represent the 95% confidence intervals for the median in a R boxplot, setting the notch argument to TRUE
axis(2, at=pretty(LegumeNLRClasses$NormalizedNumber*100),
     lab=paste0(pretty(LegumeNLRClasses$NormalizedNumber*100), "%"), las=TRUE)
# Add grey grid
grid(nx = NA, ny = NULL, col = "grey59", lty = 'dotted',
     lwd = par("lwd"), equilogs = TRUE)
# Points
stripchart(LegumeNLRClasses$NormalizedNumber*100 ~ LegumeNLRClasses$Class,
           data = LegumeNLRClasses,
           method = "jitter",
           pch = 19,
           col = c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"),
           vertical = TRUE,
           add = TRUE)

dev.off()

##### Boxplots with significant letters #####
library("multcompView")
LegumeNLRClasses_1 <- LegumeNLRClasses %>%
  mutate(Class = case_when(
    Class == "CC-NLR" ~ "CC",
    Class == "CCG10-NLR" ~ "CCG10",
    Class == "CCR-NLR" ~ "CCR",
    Class == "TIR-NLR" ~ "TIR",
    TRUE ~ Class
  ))
anova_1 <- aov(NormalizedNumber ~ Class, data = LegumeNLRClasses_1)
summary(anova_1)
tukey_1 <- TukeyHSD(anova_1)
print(tukey_1)
cld_1 <- multcompLetters4(anova_1, tukey_1) # compact letter display
print(cld_1)
# table with factors and 3rd quantile
Tk_1 <- group_by(LegumeNLRClasses_1, Class) %>%
  summarise(mean=mean(NormalizedNumber), quant = quantile(NormalizedNumber, probs = 0.75)) %>%
  arrange(desc(mean))
# extracting the compact letter display and adding to the Tk table
cld_1 <- as.data.frame.list(cld_1$Class)
Tk_1$cld <- cld_1$Letters
print(Tk_1)
LegumeNLRClasses_1$Class <- factor(LegumeNLRClasses_1$Class,
                                 levels=c("CC", "TIR",
                                          "CCG10", "CCR", "TNP", "Other"))
library("scales")
ggplot(LegumeNLRClasses_1, aes(x=Class, y=NormalizedNumber, fill=Class)) + 
  geom_boxplot() +
  labs(x="Class", y="Percentage of Protein-coding loci encoding for NLRs") +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Percentage of NLR-coding loci in Legumes") +
  geom_text(data = Tk_1, aes(x = Class, y = quant, label = cld), size = 3, vjust=-1, hjust =-1) +
  #scale_fill_brewer(c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
  scale_fill_manual(values=c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
# saving the final figure
ggsave("Output/Percentage of NLR-coding loci in Legumes 04122024_1.png", width = 1250, height = 720, units="px", dpi = 200)
ggsave("Output/Percentage of NLR-coding loci in Legumes 04122024_1.svg", width = 1250, height = 720, units="px", dpi = 200)

# NLR Number
png(filename="Output/Legume NLR Number across Classes 04122024.png",
    width = 854, height = 480)

boxplot(LegumeNLRClasses$NLRNumber ~ LegumeNLRClasses$Class,
        data = LegumeNLRClasses, col = "white",
        ylab="Number of NLRs",
        xlab="NLR Class", main="Legume NLR Number across Classes", notch = FALSE)
# Add white grid
grid(nx = NA, ny = NULL, col = "grey59", lty = 'dotted',
     lwd = par("lwd"), equilogs = TRUE)
# Points
stripchart(LegumeNLRClasses$NLRNumber ~ LegumeNLRClasses$Class,
           data = LegumeNLRClasses,
           method = "jitter",
           pch = 19,
           col = c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"),
           vertical = TRUE,
           add = TRUE)
dev.off()

anova_2 <- aov(NLRNumber ~ Class, data = LegumeNLRClasses_1)
summary(anova_2)
tukey_2 <- TukeyHSD(anova_2)
print(tukey_2)
cld_2 <- multcompLetters4(anova_2, tukey_2) # compact letter display
print(cld_2)
# table with factors and 3rd quantile
Tk_2 <- group_by(LegumeNLRClasses_1, Class) %>%
  summarise(mean=mean(NLRNumber), quant = quantile(NLRNumber, probs = 0.75)) %>%
  arrange(desc(mean))
# extracting the compact letter display and adding to the Tk table
cld_2 <- as.data.frame.list(cld_2$Class)
Tk_2$cld <- cld_2$Letters
print(Tk_2)
LegumeNLRClasses_1$Class <- factor(LegumeNLRClasses_1$Class,
                                   levels=c("CC", "TIR",
                                            "CCG10", "CCR", "TNP", "Other"))
library("scales")
ggplot(LegumeNLRClasses_1, aes(x=Class, y=NLRNumber, fill=Class)) + 
  geom_boxplot() +
  labs(x="NLR Class", y="Number of NLRs") +
  theme_bw() + 
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Legume NLR Number across Classes") +
  geom_text(data = Tk_2, aes(x = Class, y = quant, label = cld), size = 3, vjust=-1, hjust =-1) +
  #scale_fill_brewer(c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
  scale_fill_manual(values=c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
# saving the final figure
ggsave("Output/Legume NLR Number across Classes 04122024_2.png", width = 1250, height = 720, units="px", dpi = 200)
ggsave("Output/Legume NLR Number across Classes 04122024_2.svg", width = 1250, height = 720, units="px", dpi = 200)


# NLR Percentage
png(filename="Output/Distribution of Legume NLRs across Classes 04122024.png",
    width = 854, height = 480)
boxplot(LegumeNLRClasses$ClassPercentage ~ LegumeNLRClasses$Class,
        data = LegumeNLRClasses, col = "white",
        ylab="NLR Class Percentage",
        xlab="NLR Class", main="Distribution of Legume NLRs across Classes", notch = FALSE, yaxt="n")
# You can represent the 95% confidence intervals for the median in a R boxplot, setting the notch argument to TRUE
axis(2, at=pretty(LegumeNLRClasses$ClassPercentage), lab=paste0(pretty(LegumeNLRClasses$ClassPercentage), "%"), las=TRUE)
# Add white grid
grid(nx = NA, ny = NULL, col = "grey59", lty = 'dotted',
     lwd = par("lwd"), equilogs = TRUE)
# Points
stripchart(LegumeNLRClasses$ClassPercentage ~ LegumeNLRClasses$Class,
           data = LegumeNLRClasses,
           method = "jitter",
           pch = 19,
           col = c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"),
           vertical = TRUE,
           add = TRUE)
dev.off()

anova_3 <- aov(ClassPercentage ~ Class, data = LegumeNLRClasses_1)
summary(anova_3)
tukey_3 <- TukeyHSD(anova_3)
print(tukey_3)
cld_3 <- multcompLetters4(anova_3, tukey_3) # compact letter display
print(cld_3)
# table with factors and 3rd quantile
Tk_3 <- group_by(LegumeNLRClasses_1, Class) %>%
  summarise(mean=mean(ClassPercentage), quant = quantile(ClassPercentage, probs = 0.75)) %>%
  arrange(desc(mean))
# extracting the compact letter display and adding to the Tk table
cld_3 <- as.data.frame.list(cld_3$Class)
Tk_3$cld <- cld_3$Letters
print(Tk_3)
LegumeNLRClasses_1$Class <- factor(LegumeNLRClasses_1$Class,
                                   levels=c("CC", "TIR",
                                            "CCG10", "CCR", "TNP", "Other"))
library("scales")
ggplot(LegumeNLRClasses_1, aes(x=Class, y=ClassPercentage, fill=Class)) + 
  geom_boxplot() +
  labs(x="NLR Class", y="NLR Class Percentage") +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent_format(scale=1)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Distribution of Legume NLRs across Classes") +
  geom_text(data = Tk_3, aes(x = Class, y = quant, label = cld), size = 3, vjust=-1, hjust =-1) +
  #scale_fill_brewer(c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
  scale_fill_manual(values=c("turquoise","salmon","lightgoldenrod2","thistle3","dodgerblue","palegreen3"))
# saving the final figure
ggsave("Output/Distribution of Legume NLRs across Classes 04122024_3.png", width = 1250, height = 720, units="px", dpi = 200)
ggsave("Output/Distribution of Legume NLRs across Classes 04122024_3.svg", width = 1250, height = 720, units="px", dpi = 200)

##### Ratio between TIR and CC NLRs ############################################
NLRClassesDeduplicatedNumber <- read.csv("NLR Classes Normalized 06032024.csv",row.names = 1)

# Filtering the main table according to NLR class
NLRCCNumber <- NLRClassesDeduplicatedNumber %>%
  filter(Class=="CC-NLR") %>%
  filter(Species!="Zea mays" & Species!="Oryza sativa")
NLRTIRNumber <- NLRClassesDeduplicatedNumber %>%
  filter(Class=="TIR-NLR")
NLRCCG10Number <- NLRClassesDeduplicatedNumber %>%
  filter(Class=="CCG10-NLR") %>%
  filter(Species!="Zea mays" & Species!="Oryza sativa")
NLRCCRNumber <- NLRClassesDeduplicatedNumber %>%
  filter(Class=="CCR-NLR") %>%
  filter(Species!="Zea mays" & Species!="Oryza sativa")
NLRTNPNumber <- NLRClassesDeduplicatedNumber %>%
  filter(Class=="TNP") %>%
  filter(Species!="Zea mays" & Species!="Oryza sativa")

# Fetching each number of NLRs in each class
TIR <- NLRTIRNumber$NLRNumber
CC <- NLRCCNumber$NLRNumber
CCG10 <- NLRCCG10Number$NLRNumber
CCR <- NLRCCRNumber$NLRNumber
TNP <- NLRTNPNumber$NLRNumber

# Make a data frame with legume classes as columns
NLRs <- data.frame(Species = LegumePhylogeny[3:length(LegumePhylogeny)],TIR,CC,CCG10,CCR,TNP)

# Making ratios between NLR class numbers
NLRs$TIRCC <- NLRCCNumber$NLRNumber/NLRTIRNumber$NLRNumber
NLRs$CCG10CCR <- NLRCCG10Number$NLRNumber/NLRCCRNumber$NLRNumber
NLRs$TIRCCR <- NLRTIRNumber$NLRNumber/NLRCCRNumber$NLRNumber

NLRTIRCC <- data.frame(Species=c("Arabidopsis thaliana","Solanum lycopersicum",
  "Cercis canadensis", "Prosopis alba", "Prosopis cineraria", "Lupinus angustifolius","Lupinus albus",
  "Aeschynomene evenia","Arachis duranensis","Arachis ipaensis","Arachis hypogaea",
  "Cajanus cajan","Abrus precatorius","Vigna radiata",
  "Vigna angularis","Vigna unguiculata", "Phaseolus acutifolius",
  "Phaseolus vulgaris",
  "Phaseolus lunatus",
  "Glycine soja",
  "Glycine max", "Lotus japonicus","Cicer arietinum", "Medicago truncatula",
  "Trifolium pratense","Lens culinaris","Vicia sativa" ,"Vicia faba", 
  "Pisum sativum", "Lathyrus sativus"), TIR, CC, NLRs$TIRCC)
NLRTIRCC <- NLRTIRCC %>% 
  rename("TIRCC" = "NLRs.TIRCC")
NLRTIRCC$Ratio <- 1- NLRTIRCC$TIRCC
NLRTIRCC <- NLRTIRCC %>%
  mutate(cond = case_when(
  Ratio<0 ~ 'turquoise',
  Ratio>0 ~ 'salmon',
  TRUE ~ 'yellow'   #anything that does not meet the criteria above
  ))
write.csv(NLRTIRCC, "NLRTIRCC 06032024.csv")
NLRTIRCC <- read.csv("NLRTIRCC 06032024.csv")

#FIX this
#annotation <- data.frame(
#  x = c(1:30),
#  y = c(rep(-3.5,30),rep(2,30)),
#  label = c(NLRTIRCC$TIR, NLRTIRCC$CC),
#  color = c(rep("salmon",30),rep("turquoise",30))
#)

#Ratio between TIR-NLRs and CC-NLRs
NLRTIRCC %>% #Total number of NLRs
  ggplot(aes(x = Species, y = Ratio, fill=cond)) +
  geom_bar(position="stack", stat="identity") + 
  coord_flip() + geom_hline(yintercept=c(-1,0,1), linetype="dashed",
                            color = "black", linewidth=0.5) +
  scale_fill_identity() + theme_classic() + #geom_text(data=annotation, aes(x=x, y=y, label=label),
                                            #                               color=annotation$color,
                                            #          size=3,
                                            #          angle=0, fontface="bold") +
  labs(title="Ratio between CC-NLRs and TIR-NLRs in Legumes", y="1 - CC-NLR/TIR-NLR Ratio", x="Species") +
  #scale_y_continuous(limits=c(-3.5,2.5)) +
  scale_x_discrete(limits = LegumePhylogeny[3:32])
ggsave("Output/CC vs TIR NLRs 06032024.png", width = 1920, 
       height = 1080, units = "px")


#Ratio between CCG10-NLRs and CCR-NLRs
NLRs %>% #Total number of NLRs
  ggplot(aes(x = Species, y = CCG10CCR, fill="lightgoldenrod2")) +
  geom_bar(position="stack",stat="identity"#, fill=NLRTIRCC$cond
           ) + coord_flip() + geom_hline(yintercept=c(-1,0,1), linetype="dashed",
                                                                                             color = "black", size=0.5) +
  scale_fill_identity() + theme_classic() + #geom_text(data=annotation, aes(x=x, y=y, label=label),
                                                      #color=annotation$color,
                                                      #size=3,
                                                      #angle=0, fontface="bold") +
  labs(title="Ratio between CCG10-NLRs and CCR-NLRs in Legumes", y="CCG10-NLR/CCR-NLR Ratio", x="Species") +
  #scale_y_continuous(limits=c(-3.5,2.5)) +
  scale_x_discrete(limits = LegumePhylogeny[3:length(LegumePhylogeny)])
ggsave("Output/CCG10 vs CCR NLRs 06032024.png", width = 1920, 
       height = 1080, units = "px")

#Ratio between TIR and CCR-NLRs
NLRs %>% #Total number of NLRs
  ggplot(aes(x = Species, y = TIRCCR, fill="salmon")) +
  geom_bar(position="stack",stat="identity"#, fill=NLRTIRCC$cond
  ) + coord_flip() + geom_hline(yintercept=c(-1,0,1), linetype="dashed",
                                color = "black", size=0.5) +
  scale_fill_identity() + theme_classic() + #geom_text(data=annotation, aes(x=x, y=y, label=label),
  #color=annotation$color,
  #size=3,
  #angle=0, fontface="bold") +
  labs(title="Ratio between TIR-NLRs and CCR-NLRs in Legumes", y="TIR-NLR/CCR-NLR Ratio", x="Species") +
  #scale_y_continuous(limits=c(-3.5,2.5)) +
  scale_x_discrete(limits = LegumePhylogeny[3:length(LegumePhylogeny)])
ggsave("Output/TIR vs CCR NLRs 06032024.png", width = 1920, 
       height = 1080, units = "px")

library("ggExtra")
library("ggpubr")

ggscatter(NLRs, x = "TIR", y = "CC", color = "darkgreen", shape = 21, size = 3,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = NLRs$Species,
          #repel=TRUE,
          #font.label(c(5,"plain")),
          xlab = "Number of TIR-NLRs", ylab = "Number of CC-NLRs",
          title = "Correlation between the number of TIR-NLRs and CC-NLRs") #+
  #xlim(c(15000, 105000)) +
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)
ggsave("Output/Correlation between number of TIR and CC-NLRs 06032024.png", width = 1920, 
       height = 1080, units = "px")

ggscatter(NLRs, x = "CCR", y = "CCG10", color = "darkgreen", shape = 21, size = 3,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = NLRs$Species,
          #repel=TRUE,
          #font.label(c(5,"plain")),
          xlab = "Number of CCR-NLRs", ylab = "Number of CCG10-NLRs",
          title = "Correlation between the number of CCR-NLRs and CCG10-NLRs")
ggsave("Output/Correlation between number of CCG10 and CCR 06032024.png", width = 1920, 
       height = 1080, units = "px")

ggscatter(NLRs, x = "CCR", y = "TIR", color = "darkgreen", shape = 21, size = 3,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = NLRs$Species,
          #repel=TRUE,
          #font.label(c(5,"plain")),
          xlab = "Number of CCR-NLRs", ylab = "Number of TIR-NLRs",
          title = "Correlation between the number of CCR-NLRs and TIR-NLRs")
ggsave("Output/Correlation between number of CCR and TIR 06032024.png", width = 1920, 
       height = 1080, units = "px")

ggscatter(NLRs, x = "CCG10", y = "CC", color = "darkgreen", shape = 21, size = 3,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = NLRs$Species,
          #repel=TRUE,
          #font.label(c(5,"plain")),
          xlab = "Number of CCG10-NLRs", ylab = "Number of CC-NLRs",
          title = "Correlation between the number of CCG10-NLRs and CC-NLRs")
ggsave("Output/Correlation between number of CCG10 and CC 06032024.png", width = 1920, 
       height = 1080, units = "px")
##### Number of NLRs vs Number of Protein-coding loci ##############
LegumeNLRClasses <- read.csv("Legume NLR Classes 06032024.csv",row.names = 1)
#library(ggpmisc)
library("ggExtra")
library("ggpubr")


png(filename = "Output/Correlation between NLR Number and Number of Protein-coding Loci 04122024 2.png", 
    width = 1280, height = 720, units = "px", res =120)
#scatter_TotalNumberNLR_vs_ProteinCodingLoci <- 
  LegumeNLRClasses %>% #Normalized NLRs
  ggplot(aes(x = count , y = TotalNLRNumber,
             #color="darkgreen"
  )) +
  #by multiplying NLR Number/Number of Loci by 100, we get the percentage of NLR genes!
  theme_classic(base_size = 8)  +
  geom_point() + geom_smooth(method='lm',
                             se=TRUE, 
                             #formula = my.formula, 
                             fullrange=TRUE) + #scale_y_continuous(trans = 'log10') +
  #scale_x_continuous(trans = 'log2') + 
  #stat_poly_eq(formula = my.formula, 
  #             aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #             parse = TRUE) +
  xlim(c(15000,80000)) + ylim(c(0,1100)) +geom_text(aes(
    #colour = factor(wild)
    ), label=LegumeNLRClasses$Species, check_overlap = FALSE, 
    position = position_dodge(width=0.9),
    vjust = -0.75,
    hjust = 1, #"inward", 
    angle = -30, #size = 10,
    #nudge_y= 15
    ) +
  labs(title="Correlation between NLR Number and Number of Protein-coding Loci", y="NLR Number", x="Number of Protein-coding Loci") +
  #scale_color_brewer(palette="Dark2")
  annotate(geom="text",x = 25000, y = c(900,800), label = c("R^2 == 0.35", "p == 3.3*e^{-15}"), fontface = 'italic', parse = TRUE)
scatter_TotalNumberNLR_vs_ProteinCodingLoci 
                                                       # )
#p1 <- 
  ggMarginal(scatter_TotalNumberNLR_vs_ProteinCodingLoci, type="boxplot", size=10)
ggsave("Output/Correlation between NLR Number and Number of Protein-coding Loci 04122024.png", width = 1920, 
       height = 1080, units = "px", dpi = 125)
ggsave("Output/Correlation between NLR Number and Number of Protein-coding Loci 04122024.svg", width = 1920, 
       height = 1080, units = "px", dpi = 125)
dev.off()

p <- ggscatter(LegumeNLRClasses, x = "count", y = "TotalNLRNumber", color = "darkgreen", shape = 21, size = 3,
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson", label = LegumeNLRClasses$Species,
               #repel=TRUE,
               #font.label(c(5,"plain")),
               xlab = "Number of Protein-coding Loci", ylab = "Total NLR Number",
               title = "Total NLR Number according to Protein-coding Loci") +
               xlim(c(15000, 80000)) #+
               #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)
p
ggsave("Output/Total NLR Number according to Protein-coding Loci 04122024.png", width = 1920, 
       height = 1080, units = "px")
ggsave("Output/Total NLR Number according to Protein-coding Loci 04122024.svg", width = 1920, 
       height = 1080, units = "px")
#"R^2 = 0.33, p = 2,8e-12"

