library(msa)
library(ape)
library(seqinr)
library(tidyverse)

# read in fasta file
fcontent <- readLines("fasta/subseq_anisakidae_coi5p.fasta")

seqStart <- grep('>', fcontent)
seqEnd <- seqStart - 1
seqEnd <- seqEnd[-1]
seqEnd <- c(seqEnd, length(fcontent))

mSeqPos <- cbind(seqStart, seqEnd)
rownames(mSeqPos) <- fcontent[seqStart]
mSeqPos[, 1] <- mSeqPos[, 1] + 1
mSeqPos <- cbind(mSeqPos, sequence=NA)

for (x in 1:nrow(mSeqPos)){
        sequ <- ''
        for (l in mSeqPos[x, 'seqStart']:mSeqPos[x, 'seqEnd']){
                sequ <- paste0(sequ, fcontent[l])
        }
        mSeqPos[x, 'sequence'] <- sequ
}

seqPos <- mSeqPos %>% 
        as.data.frame() %>% 
        rownames_to_column('id') %>% 
        filter(grepl('Anisakis|Pseudoterranova|Phocascaris|Contracaecum', id)) %>%
        mutate(id=gsub('cf. ', '', id)) %>% 
        mutate(id=paste(gsub('>', '', word(id, 1)), word(id, 2), word(id, 3))) %>% 
        mutate(genus=word(id, 2)) %>% 
        mutate(species=paste(word(id, 2), word(id,3))) %>% 
        select(id, genus, species, sequence) %>% 
        column_to_rownames('id')

# Step 1: take the DNA sequence
DNASEQ <- seqPos[, 'sequence']
names(DNASEQ) <- rownames(seqPos)
# Step 2: convert to DNAStringSet
DNASS <- DNAStringSet(DNASEQ)
# or directly read in fasta using readDNAStringSet
# DNASS <- readDNAStringSet("fasta/subseq_anisakidae_18s.fasta")
# DNASS <- readDNAStringSet("fasta/subseq_anisakidae_coi5p.fasta")

# Step 3: multiple sequence alignment
alignment <- msa(DNASS, 'ClustalOmega')
# Step 4: clustering and distances
alignment <- msaConvert(alignment, "seqinr::alignment")
distMatrix <- dist.alignment(alignment, 'similarity')
clustering <- hclust(distMatrix)
# Step 5: APE to make cladogram
phylotree <- as.phylo(clustering)
# plot(phylotree, type="cladogram")

ggtree(phylotree) +
        geom_tiplab(size=2.8) +
        ggexpand(ratio=1.2, side='h')


