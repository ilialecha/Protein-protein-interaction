library(ggplot2)

setwd("~/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files")

data <- read.delim("~/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files/ala_interaction_energies.tsv",
                                       header=TRUE, stringsAsFactors=TRUE)

data2 <- read.delim("~/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files/interaction_energies.tsv",
                   header=TRUE, stringsAsFactors=TRUE)

print(c(data2$RES_ID[data2$CHAIN_ID=="A"]),colapse=" or ")
data2$RES_ID[data2$CHAIN_ID=="E"]

m0j <- read.delim("~/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files/6m0j_biopython.tsv", stringsAsFactors=TRUE)

m0j_A <- m0j[m0j$CHAIN_ID == "A", ]
m0j_E <- m0j[m0j$CHAIN_ID == "E", ]


ggplot(m0j, aes(m0j$RES, fill = m0j$CHAIN_ID)) + geom_bar(position = 'identity', alpha = 0.3) + labs(fill = "Chain") + labs(x="Residue num (3 Code nomenclature)") + labs(y="Count")


data_A <- data[data$CHAIN.ID == "A", ]
data_E <- data[data$CHAIN.ID == "E", ]


library(gghighlight)

ggplot(data_A,aes(paste(data_A$RES.ID,data_A$RES.Num),data_A$Total.AAG.Change, fill=data_A$Total.AAG.Change))+gghighlight(max(data_A$Total.AAG.Change) > 0)+labs(fill = "Difference")+ geom_hline(yintercept=mean(data_A$Total.AAG.Change), linetype="dashed", color = "red") +geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "Residue treated as Alanine.") + labs(y = "Relative change on ΔΔG") + labs(title = "Alanine-Scanning on interaction surface resdiues: Chain A")


ggplot(data_E,aes(paste(data_E$RES.ID,data_E$RES.Num),data_E$Total.AAG.Change,fill=data_E$Total.AAG.Change))+labs(fill = "Difference")+ geom_hline(yintercept=mean(data_E$Total.AAG.Change), linetype="dashed", color = "red") +geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "Residue treated as Alanine.") + labs(y = "Relative change on ΔΔG") + labs(title = "Alanine-Scanning on interaction surface resdiues: Chain E")

