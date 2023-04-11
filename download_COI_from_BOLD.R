# install packages (if necessary)
# install.packages("bold")
# install.packages("taxize")
# install.packages("seqinr")
# load packages
library(tidyverse)
library(bold)    # API interface to BOLD
library(taxize)  # for NCBI taxonomy lookup
library(seqinr)  # for FASTA output

# get class-level taxa within "Animalia" from NCBI taxonomy
# taxa <- downstream("Animalia", db = "ncbi", downto = "class")
taxa <- downstream("Anisakidae", db="ncbi", downto="species")
taxa <- downstream("Anisakidae", db="ncbi", downto="genus")
# check if taxa present in BOLD
checks <- bind_rows(bold_tax_name(taxa$Anisakidae$childtaxa_name[1:71]),
                    bold_tax_name(taxa$Anisakidae$childtaxa_name[74:78]))
checks <- bold_tax_name(taxa$Anisakidae$childtaxa_name)

taxa_bold <- checks[!is.na(checks$taxon),]$taxon

# download sequences from BOLD for each class-level taxon
sequences_species <- map(taxa_bold, bold_seq, marker='COI-5P') %>%
        flatten() %>%
        bind_rows() %>% 
        mutate(sequence=gsub("-", "", sequence))

# write sequences to file
write.fasta(
        sequences = as.list(sequences_species$sequence), 
        names = paste(as.list(sequences_species$id), as.list(sequences_species$name)), 
        nbchar = 100, 
        file.out = 'anisakidae_species_coi5p.fasta')
