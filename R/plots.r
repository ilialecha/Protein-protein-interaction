library(ggplot2)

setwd("~/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files")

ala_interaction_energies <- read.delim("~/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files/ala_interaction_energies.tsv",
                                       header=FALSE, stringsAsFactors=TRUE)

data_GLY <- ala_interaction_energies[ala_interaction_energies$V1 in "ALA"]
attach(ala_interaction_energies)
ggplot(ala_interaction_energies,aes(ala_interaction_energies$V1,ala_interaction_energies$V6))+geom_bar(stat="identity")
