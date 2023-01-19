# Select working directory
setwd("C://Users/Rita/Dropbox (KamounLab)/Kamounity folder/Rita/NLRtracker/AnalysisNLRtracker131/")
setwd("D:/Dropbox (KamounLab)/Kamounity folder/Rita/NLRtracker/AnalysisNLRtracker131")
#Load in the required packages
library(readxl)
library(tidyverse)
library(RColorBrewer)

# 1. Align NB ARC files from the NLR tracker output with the RefPlantNLR
# datasets to correctly assign each NLR to its respective class.
# 2. Create the NLR Classes w RefPlantNLR file with the following columns:
# Species seqname Class
df1 <- read_xlsx("NLR Classes Legumes w RefPlantNLR.xlsx")
RefPlantNLRdelete <- read_xlsx("RefPlantNLR.xlsx") #to remove RefPlantNLR NLR references
LOCUS <- read_tsv("Fabales_LOCUS_09012023.tsv") #for deduplicating

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
                   "Prosopis alba","Prosopis cineraria", "Pisum sativum", 
                   "Pisum sativum ZW6",
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
  "Pisum sativum ZW6", 
  "Lathyrus sativus")

# filter out every row whose
# Seqname matches with RefPlantNLR's seqname
df2 <- df1 %>%
  filter(!seqname %in% RefPlantNLRdelete$seqname) 
# https://www.datasciencemadesimple.com/delete-or-drop-rows-in-r-with-conditions-2/
# Dropping rows in R
write_csv(df2,"NLR Classes Legumes Filtered.csv")
df2 <- read_csv("NLR Classes Legumes Filtered.csv")

# Count how many NLRs are in each class, for each Species
NLRClasses <- df2 %>%
  group_by(Species) %>%
  count(Class) 
write_csv(NLRClasses,"NLR Number per Class.csv")

############### Generating graphs from non-deduplicated data ##################
NLRClasses <- NLRClasses %>%
  group_by(Species) %>%
  mutate(ClassPercentage = n/(sum(n))*100) #Percentage of NLRs with that class
View(NLRClasses)
write_csv(NLRClasses, "NLRClassesNotDeduplicated.csv")

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

############### Deduplicating alternate splicing variants #####################

# Removing the Splicing alternatives to get the Loci 
# (i.e. remove every character after the last dot)
#df2 <- df2 %>%
#  select(df2$Species=="Glycine max"|df2$Species=="Lens culinaris") %>%
#  mutate(df2$seqname, gsub("[.][^.]+$", "", df2$seqname))

# Join the NLR data with the locus data
#df2 <- df2 %>%
#  filter(Species!="Pisum sativum ZW6" & Species!="Cercis canadensis" & Species!="Vicia faba"
#         & Species!="Vicia sativa" & Species!="Prosopis cineraria")
df3 <- left_join(df2,LOCUS, by="seqname")
write_csv(df3,"df3.csv")
df3 <- read_csv("df3.csv")
NLR_LOCUS <- df3 %>%
  group_by(Genome,Locus) %>%
  arrange(desc(seqname)) %>%
  slice(1)
  #case_when(Species != "Pisum sativum ZW6" ~ slice(1))

View(NLR_LOCUS)
write_csv(NLR_LOCUS,"Locus and Type Legume NLRs.csv")
NLR_LOCUS <- read_csv("Locus and Type Legume NLRs.csv")


######################## Get how many deduplicated NLRs are in each type ######################################
NLRClassesDeduplicatedNumber <- NLR_LOCUS %>%
  group_by(Species) %>%
  count(Class)
# Get the percentages of NLRs in each class
NLRClassesDeduplicatedNumber <- NLRClassesDeduplicatedNumber %>%
  group_by(Species) %>%
  mutate(ClassPercentage = n/(sum(n))*100) #Percentage of NLRs with that class
View(NLRClassesDeduplicatedNumber)
write_csv(NLRClassesDeduplicatedNumber, "NLRClassesDeduplicatedNumber.csv")

NLRClassesDeduplicatedNumber %>% #Total number of NLRs
  ggplot(aes(x = Species, y = n, fill = Class)) +
  geom_bar(position="stack",stat="identity") + coord_flip() + theme_classic(base_size = 8) +
  scale_fill_manual(values=c("turquoise","lightgoldenrod2","thistle3","palegreen3","salmon","dodgerblue")) +
  labs(title="Classes of Legume NLRs", y="NLR Number", x="Species", fill="Classes") +
  scale_x_discrete(limits = LegumePhylogeny)
#colors for ggplot http://sape.inf.usi.ch/quick-reference/ggplot2/colour

NLRClassesDeduplicatedNumber %>% #Percentages
  ggplot(aes(x = Species, y = ClassPercentage, fill = Class)) +
  geom_bar(position="stack",stat="identity") + coord_flip() +
  theme_classic() +
  scale_fill_manual(values=c("turquoise","lightgoldenrod2","thistle3","palegreen3","salmon","dodgerblue")) +
  labs(title="Percentages of NLR Classes in Legumes", y="NLR Class Percentage", x="Species", fill="Classes") +
  scale_x_discrete(limits = LegumePhylogeny)

##Normalize the NLR Number with the number of different loci in Fabales_LOCUS ##
NLRClassesDeduplicatedNumber <- read_csv("NLRClassesDeduplicatedNumber.csv")
LOCUS <- read_tsv("Fabales_LOCUS_09012023.tsv")

LOCUS <- LOCUS %>%
  #add CAAS_Psat_ZW6_1.0 to the Genome for Pisum sativum ZW6 
  mutate(Genome = case_when(Organism == "Pisum sativum ZW6" ~ "CAAS_Psat_ZW6_1.0",
                            TRUE ~ Genome) )

LOCI_Number <- count(LOCUS,Genome) #count the amount of loci each genome has

LOCI_Genome <- LOCI_Number$Genome #get the genome IDs
View(LOCI_Genome)
# Genome IDs
#"Abrus_2018", "Aeschynomene_evenia_v1.0.protein", "Arabidopsis thaliana TAIR 10.1", "Aradu1.1",
#"arahy.Tifrunner.gnm1.KYV3", "Araip1.1", "ASM33114v1", "ASM411807v1", "ASM419377v2", "ASM479914v1", "C.cajan_V1.0",
#"CAAS_Psat_ZW6_1.0", "cerca.ISC453364.gnm1.ann1.B05Z",
#"cicar.CDCFrontier.gnm1.ann1.nRhs.protein", "Glycine max Williams 82", "Glycine_max_v4.0",
#"LATSA3860_EIv1.0.annotation.gff3.pep", "Lens culinaris Redberry v2.0", "lotja.MG20.gnm3.ann1.WF9B.protein",
#"lupal.Amiga.gnm1.ann1.3GKS.protein", "lupan.Tanjil.gnm1.ann1.nnV9.protein", "LupAngTanjil_v1.0", "MtrunA17r5.0-ANR",
#"Oryza sativa Nipponbare IRGSP 1.0", "Pacutifolius_580_v1.0.protein", "PacutifoliusWLD_581_v2.0.protein", "PC_final"
#"phavu.G19833.gnm2.ann1.PB8d.protein", "PhaVulg1_0", "pissa.Cameor.gnm1.ann1.7SZR.protein", "Pisum_sativum_v1a_prot", 
#"Plunatus_563_V1.protein", "Solanum lycopersicum SL3.0", "tripr.MilvusB.gnm2.ann1.DFgp.protein",
#"Vicia faba HEDIN_TMP5", "Vicia sativa",
#"vigan.Gyeongwon.gnm3.ann1.3Nz5.protein", "Vigan1.1", "vigra.VC1973A.gnm6.ann1.M1Qs.protein", "Vradiata_ver6",
#"Zea mays B73 v5"

LOCI_Organism <- c("Abrus precatorius", "Aeschynomene evenia", "Arabidopsis thaliana", "Arachis duranensis",
"Arachis hypogaea", "Arachis ipaensis", "", "Vigna unguiculata", "Glycine soja", "Prosopis alba", "Cajanus cajan",
"Pisum sativum ZW6", "Cercis canadensis", "Cicer arietinum", "Glycine max", "",
"Lathyrus sativus", "Lens culinaris", "Lotus japonicus",
"Lupinus albus", "Lupinus angustifolius", "", "Medicago truncatula",
"Oryza sativa", "", "Phaseolus acutifolius", "Prosopis cineraria",
"Phaseolus vulgaris", "", "", "Pisum sativum", "Phaseolus lunatus", "Solanum lycopersicum", "Trifolium pratense",
                   "Vicia faba", "Vicia sativa","Vigna angularis", "", "Vigna radiata", "", "Zea mays")

LOCI_Number <- LOCI_Number %>%
  mutate(Organism = LOCI_Organism) %>% #add a new column with the Organism name
  filter(Organism != "") %>% #get rid of the genomes which weren't used
  arrange(Organism) #arrange the table by Organism name

write.csv(LOCI_Number, "Loci Number for Each Species.csv")
colnames(NLRClassesDeduplicatedNumber) <- c("Species","Class","NLRNumber","ClassPercentage") #rename the columns
NLRClassesDeduplicatedNumber <- left_join(NLRClassesDeduplicatedNumber,LOCI_Number, by=c("Species" = "Organism"))
#join the deduplicated numbers with the LOCI_Number

NLRClassesDeduplicatedNumber <- NLRClassesDeduplicatedNumber %>%
  group_by(Species) %>%
  mutate(NormalizedNumber = NLRNumber/n) #divide the NLR Numbers for the total protein coding loci of each genome
View(NLRClassesDeduplicatedNumber)

NLRClassesDeduplicatedNumber <- NLRClassesDeduplicatedNumber %>%
  arrange(factor(Species,levels=LegumePhylogeny),Class) #arrange the table according to phylogeny and then NLR class

NLRClassesDeduplicatedNumber <- NLRClassesDeduplicatedNumber %>%
  group_by(Species) %>%
  mutate(TotalNLRNumber = sum(NLRNumber)) #get the total amount of NLRs in each species

write.csv(NLRClassesDeduplicatedNumber, "NLR Classes Normalized.csv") #save the table
NLRClassesDeduplicatedNumber <- read.csv("NLR Classes Normalized.csv")

NLR_number <- unique(NLRClassesDeduplicatedNumber$TotalNLRNumber) #get a list of the total number of NLRs
NLR_number <- c(NLR_number[1:2],rep(NLR_number[3], 2),NLR_number[4:32]) #Arabidopsis thaliana and Solanum lycopersicum
#had the same total number of NLRs, so I had to manually repeat (rep) 170 twice in the list
write.csv(NLRClassesDeduplicatedNumber, "NLR Classes Normalized.csv") #save the table


NLRClassesDeduplicatedNumber %>% #Normalized NLRs
  ggplot(aes(x = Species, y = NormalizedNumber*100, fill = Class)) +
  #by multiplying NLR Number/Number of Loci by 100, we get the percentage of NLR genes!
  theme_classic(base_size = 8) +
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
           y = 1.65, x = seq_along(LegumePhylogeny),
           size = 0.75*4*1,
           colour = "black",
           hjust=0,
           family="sans"
  )
# colors for ggplot http://sape.inf.usi.ch/quick-reference/ggplot2/colour
# try to annotate the NLR Numbers in the graph

############# Boxplots for Normalized Number of NLRs #######################
#https://r-charts.com/distribution/add-points-boxplot/
NLRClassesDeduplicatedNumber <- read.csv("NLR Classes Normalized.csv")

# Create a file with info exclusively for legume species (exclude the outgroups)
LegumeNLRClasses <- subset(NLRClassesDeduplicatedNumber, Species!="Oryza sativa"
  & Species!="Arabidopsis thaliana" & Species!="Solanum lycopersicum" & Species!="Zea mays")

# Vertical box plot by group
LegumeNLRClasses$Class <- factor(LegumeNLRClasses$Class,
                                 levels=c("CC-NLR", "TIR-NLR",
                                          "CCG10-NLR", "CCR-NLR", "TNP", "Other"))
View(LegumeNLRClasses)
write.csv(LegumeNLRClasses, "Legume NLR Classes.csv")

# Normalized Number boxplots
boxplot(LegumeNLRClasses$NormalizedNumber*100 ~ LegumeNLRClasses$Class,
        data = LegumeNLRClasses, col = "white",
        ylab="Percentage of Protein-coding loci encoding for NLRs",
        xlab="NLR Class", main="Percentage of NLR-coding loci in Legumes", notch = FALSE, yaxt="n")
#You can represent the 95% confidence intervals for the median in a R boxplot, setting the notch argument to TRUE
axis(2, at=pretty(LegumeNLRClasses$NormalizedNumber*100),
     lab=paste0(pretty(LegumeNLRClasses$NormalizedNumber*100), "%"), las=TRUE)
# Add white grid
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

# NLR Number
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

# NLR Percentage
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

################################ Ratio between TIR and CC NLRs ############################################
NLRClassesDeduplicatedNumber <- read.csv("NLR Classes Normalized.csv")

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
NLRs$TIRCC <- NLRTIRNumber$NLRNumber/NLRCCNumber$NLRNumber
NLRs$CCG10CCR <- NLRCCG10Number$NLRNumber/NLRCCRNumber$NLRNumber
NLRs$TIRCCR <- NLRTIRNumber$NLRNumber/NLRCCRNumber$NLRNumber

NLRTIRCC <- data.frame(Species=c("Arabidopsis thaliana","Solanum lycopersicum",
  "Cercis canadensis",
                                 "Prosopis alba", "Prosopis cineraria", "Lupinus angustifolius","Lupinus albus",
  "Aeschynomene evenia","Arachis duranensis","Arachis ipaensis","Arachis hypogaea",
  "Cajanus cajan","Abrus precatorius","Vigna radiata",
  "Vigna angularis","Vigna unguiculata", "Phaseolus acutifolius",
  "Phaseolus vulgaris",
  "Phaseolus lunatus",
  "Glycine soja",
  "Glycine max", "Lotus japonicus","Cicer arietinum", "Medicago truncatula",
  "Trifolium pratense","Lens culinaris","Vicia sativa" ,"Vicia faba", 
  "Pisum sativum","Pisum sativum ZW6", "Lathyrus sativus"), TIR, CC, NLRs$TIRCC)
NLRTIRCC <- NLRTIRCC %>% 
  rename("TIRCC" = "NLRs.TIRCC")
NLRTIRCC$Ratio <- 1- NLRTIRCC$TIRCC
NLRTIRCC <- NLRTIRCC %>%
  mutate(cond = case_when(
  Ratio<0 ~ 'salmon',
  Ratio>0 ~ 'turquoise',
  TRUE ~ 'yellow'   #anything that does not meet the criteria above
  ))
write.csv(NLRTIRCC, "NLRTIRCC.csv")
NLRTIRCC <- read.csv("NLRTIRCC.csv")

#FIX this
annotation <- data.frame(
  x = c(1:30),
  y = c(rep(-3.5,30),rep(2,30)),
  label = c(NLRTIRCC$TIR, NLRTIRCC$CC),
  color = c(rep("salmon",30),rep("turquoise",30))
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
  scale_x_discrete(limits = LegumePhylogeny[3:28])

#Ratio between CCG10-NLRs and CCR-NLRs
NLRs %>% #Total number of NLRs
  ggplot(aes(x = Species, y = CCG10CCR)) +
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

#Ratio between TIR and CCR-NLRs
NLRs %>% #Total number of NLRs
  ggplot(aes(x = Species, y = TIRCCR)) +
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

ggscatter(NLRs, x = "CCR", y = "CCG10", color = "darkgreen", shape = 21, size = 3,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = NLRs$Species,
          #repel=TRUE,
          #font.label(c(5,"plain")),
          xlab = "Number of CCR-NLRs", ylab = "Number of CCG10-NLRs",
          title = "Correlation between the number of CCR-NLRs and CCG10-NLRs")

ggscatter(NLRs, x = "CCR", y = "TIR", color = "darkgreen", shape = 21, size = 3,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", label = NLRs$Species,
          #repel=TRUE,
          #font.label(c(5,"plain")),
          xlab = "Number of CCR-NLRs", ylab = "Number of TIR-NLRs",
          title = "Correlation between the number of CCR-NLRs and TIR-NLRs")
####################### Wild vs Cultivated ##############################
LegumeNLRClasses <- read.csv("Legume NLR Classes.csv")
wild <- c("Prosopis alba", "Lupinus angustifolius",
          "Aeschynomene evenia","Arachis duranensis","Arachis ipaensis",
          "Abrus precatorius",
          "Glycine soja","Lotus japonicus","Medicago truncatula",
          "Trifolium pratense")
cultivated <- c("Lupinus albus",
                "Arachis hypogaea",
                "Cajanus cajan","Vigna radiata",
                "Vigna angularis","Vigna unguiculata", "Phaseolus acutifolius",
                "Phaseolus vulgaris", "Phaseolus lunatus",
                "Glycine max", "Cicer arietinum",
                "Lens culinaris","Pisum sativum","Lathyrus sativus")
LegumeNLRClasses <- LegumeNLRClasses %>%
  mutate(wild=case_when(
    Species %in% wild ~ "Wild",
    Species %in% cultivated ~ "Cultivated"
  ))

boxplot(LegumeNLRClasses$TotalNLRNumber ~ LegumeNLRClasses$wild,
        data = LegumeNLRClasses, col = "white",
        ylab="Total NLR Number",
        xlab="Species Wildness", main="Number of NLRs according to Species Wildness", notch = FALSE, yaxt="n")
#You can represent the 95% confidence intervals for the median in a R boxplot, setting the notch argument to TRUE
axis(2, at=pretty(LegumeNLRClasses$TotalNLRNumber), lab=paste0(pretty(LegumeNLRClasses$TotalNLRNumber), ""), las=TRUE)
# Add white grid
grid(nx = NA, ny = NULL, col = "grey59", lty = 'dotted',
     lwd = par("lwd"), equilogs = TRUE)
# Points
stripchart(LegumeNLRClasses$TotalNLRNumber ~ LegumeNLRClasses$wild,
           data = LegumeNLRClasses,
           method = "jitter",
           pch = 19,
           col = c("dodgerblue","palegreen3"),
           vertical = TRUE,
           add = TRUE)

############## Number of NLRs vs Number of Protein-coding loci ##############
LegumeNLRClasses <- read.csv("Legume NLR Classes.csv")
#library(ggpmisc)
library("ggExtra")
library("ggpubr")

#lm_eqn <- function(df){
#  m <- lm(y ~ x, df);
#  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                   list(a = format(unname(coef(m)[1]), digits = 2),
#                        b = format(unname(coef(m)[2]), digits = 2),
#                        r2 = format(summary(m)$r.squared, digits = 3)))
#  as.character(as.expression(eq));
#}

#my.formula <- LegumeNLRClasses$TotalNLRNumber ~ LegumeNLRClasses$n

scatter_TotalNumberNLR_vs_ProteinCodingLoci <- LegumeNLRClasses %>% #Normalized NLRs
  ggplot(aes(x = n , y = TotalNLRNumber,
             #color="darkgreen"
  )) +
  #by multiplying NLR Number/Number of Loci by 100, we get the percentage of NLR genes!
  theme_classic(base_size = 8)  +
  geom_point() + geom_smooth(method='lm',
                             #aes(fill=wild),
                             se=TRUE, 
                             #formula = my.formula, 
                             fullrange=TRUE) + #scale_y_continuous(trans = 'log10') +
  #scale_x_continuous(trans = 'log2') + 
  #stat_poly_eq(formula = my.formula, 
  #             aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #             parse = TRUE) +
  xlim(c(15000,105000)) + ylim(c(0,1100)) +geom_text(aes(
    #colour = factor(wild)
    ), label=LegumeNLRClasses$Species, check_overlap = FALSE, position = position_dodge(width=0.9),
    vjust = -0.75,
    hjust = 1, #"inward", 
    angle = -30, #size = 10,
    #nudge_y= 15
    ) +
  labs(title="Correlation between NLR Number and Number of Protein-coding Loci", y="NLR Number", x="Number of Protein-coding Loci") #+
  #scale_color_brewer(palette="Dark2")
scatter_TotalNumberNLR_vs_ProteinCodingLoci + annotate(geom="text",x = 25000, y = c(900,800), label = c("R^2 == 0.33", "p == 2.8*e^{-12}"), fontface = 'italic', parse = TRUE)
                                                       # )
p1 <- ggMarginal(scatter_TotalNumberNLR_vs_ProteinCodingLoci, type="boxplot", size=10)
p1

p <- ggscatter(LegumeNLRClasses, x = "n", y = "TotalNLRNumber", color = "darkgreen", shape = 21, size = 3,
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson", label = LegumeNLRClasses$Species,
               #repel=TRUE,
               #font.label(c(5,"plain")),
               xlab = "Number of Protein-coding Loci", ylab = "Total NLR Number",
               title = "Total NLR Number according to Protein-coding Loci") +
               xlim(c(15000, 105000)) +
               stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)
p
#"R^2 = 0.33, p = 2,8e-12"
#Scatter plot of nr NLRs vs Nr of Protein-coding Loci
#Highlight species with fewer NLR coding loci
#Plot data with BUSCO scores
#Analyse genomes with other tools besides BUSCO
#Check how many NLRs would be predicted using NLR Annotator vs the predicted Number of NLRs with NLRtracker
#Proportion of truncated NLRs (identified by NLRtracker but not considered) vs whole NLRs
#Bootstrap wild vs cultivated to assess significance of NLR number vs wildness
#Check if the conserved TIR-NLR is present in the legume genomes
#Think about other interesting questions I could be asking the data besides wildness/artificial selection
#QC tools for genomes
#Relate duplicated BUSCOs with ploidy
#Add two new columns to the dataset: ploidy (relate to BUSCO scores) and rhizobia type (relate to NLR number)
#Forget about PRRtracker for now, start thinking of ideas for my data. Colombus didn't think of discovering a new continent after discovering America

############################################# Patristic Distances NLRs ###########################################
# This script is to filter out NLRs from specific species to pinpoint the most conserved NLRs
NLR_LOCUS <- read_csv("Locus and Type Legume NLRs.csv")
Patristic_Distances <- NLR_LOCUS %>%
  filter(Species=="Glycine max" | Species=="Medicago truncatula" | Species=="Arachis hypogaea" |
           Species=="Vigna unguiculata" |
           Species=="Lathyrus sativus" | Species=="Lupinus albus" | Species=="Phaseolus vulgaris") %>%
  arrange(Species, Class)
Glycine_max_NLRs <- NLR_LOCUS %>%
  filter(Species=="Glycine max") %>%
  arrange(Class)
write_csv(Patristic_Distances,"Patristic Distances/Patristic_Distances_NLRs.csv")
write_csv(Glycine_max_NLRs,"Patristic Distances/Glycine_max_NLRs.csv")