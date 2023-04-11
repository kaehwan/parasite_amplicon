library(msa)
library(ape)
library(seqinr)
library(tidyverse)


# PSA: Targeted Anisakidae 18S sequence -----------------------------------

# 18S sequences were downloaded from NCBI
# Pre-filter to exclude entries without species info
# Filter to include 'Anisakis|Pseudoterranova|Phocascaris|Contracaecum'


# read in fasta file
fcontent <- readLines("fasta/subseq_anisakidae_18s.fasta")

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
        mutate(id=paste(gsub('>', '', word(id, 1)), word(id, 2), word(id, 3))) %>% 
        mutate(genus=word(id, 2)) %>% 
        mutate(species=paste(word(id, 2), word(id,3))) %>% 
        select(id, genus, species, sequence) %>% 
        column_to_rownames('id')

seqPos[1, 'species'] <- 'Anisakis pegreffii'
seqPos %>% 
        group_by(genus, species) %>% 
        count() %>% 
        ungroup() %>% 
        ggplot(aes(x=species, y=n)) + 
        geom_bar(stat='identity') + geom_text(aes(label=n), vjust=-0.25) +
        facet_grid(. ~ genus, scales='free') + 
        theme_bw() + 
        theme(axis.text.x=element_text(angle=45, hjust=1)) + 
        labs(title='18S')

mPairAlign <- matrix(data=NA, nrow=8, ncol=8)
rownames(mPairAlign) <- rownames(seqPos)
colnames(mPairAlign) <- rownames(seqPos)

for (i in 1:dim(seqPos)[1]){
        for (j in 1:dim(seqPos)[1]){
                align <- Biostrings::pairwiseAlignment(
                        seqPos[i, 'sequence'],
                        seqPos[j, 'sequence'],
                        type='local'
                )
                identity <- Biostrings::pid(align)
                mPairAlign[i, j] <- identity
        }
}

PairAlign <- mPairAlign %>% 
        as.matrix() %>% 
        as_tibble(rownames="A") %>% 
        pivot_longer(-A, names_to="B", values_to="identity")

PairAlign %>% 
        ggplot(aes(x=A, y=B, fill=identity)) +
        geom_tile() +
        geom_text(aes(label=format(round(identity, 2), nsmall=2)), 
                color=case_when(
                        PairAlign$identity == max(PairAlign$identity) ~ "red",
                        PairAlign$identity < max(PairAlign$identity) ~ "white"
                        ), 
                size=3) +
        coord_fixed() +
        theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('subseq_anisakidae_18s_pairwise.svg', width=10, height=10)



# PSA: Targeted Anisakidae COI5P sequence ---------------------------------


# COI5P sequences were downloaded from BOLD systems using "download_COI_from_BOLD.R"
# Pre-filter to exclude entries without species info
# Filter to include 'Anisakis|Pseudoterranova|Phocascaris|Contracaecum'


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
        mutate(id=paste(gsub('>', '', word(id, 1)), word(id, 2), word(gsub('cf. ', '', id), 3))) %>% 
        mutate(genus=word(id, 2)) %>% 
        mutate(species=paste(word(id, 2), word(id,3))) %>% 
        select(id, genus, species, sequence) %>% 
        column_to_rownames('id')

seqPos %>% 
        group_by(genus, species) %>% 
        count() %>% 
        ungroup() %>% 
        ggplot(aes(x=species, y=n)) + 
        geom_bar(stat='identity') + geom_text(aes(label=n), vjust=-0.25) +
        facet_grid(. ~ genus, scales='free') + 
        theme_bw() + 
        theme(axis.text.x=element_text(angle=45, hjust=1)) + 
        labs(title='COI')

n <- dim(seqPos)[1]
mPairAlign <- matrix(data=NA, nrow=n, ncol=n)
rownames(mPairAlign) <- rownames(seqPos)
colnames(mPairAlign) <- rownames(seqPos)

for (i in 1:dim(seqPos)[1]){
        for (j in 1:dim(seqPos)[1]){
                align <- Biostrings::pairwiseAlignment(
                        seqPos[i, 'sequence'],
                        seqPos[j, 'sequence'],
                        type='local'
                )
                identity <- Biostrings::pid(align)
                mPairAlign[i, j] <- identity
        }
}

PairAlign <- mPairAlign %>% 
        as.matrix() %>% 
        as_tibble(rownames="A") %>% 
        pivot_longer(-A, names_to="B", values_to="identity")

PairAlign %>% 
        ggplot(aes(x=A, y=B, fill=identity)) +
        geom_tile() +
        geom_text(aes(label=round(identity, 2)), 
                  color=case_when(
                          PairAlign$identity == max(PairAlign$identity) ~ "red",
                          PairAlign$identity < max(PairAlign$identity) ~ "white"
                  ), 
                  size=1) +
        coord_fixed() +
        theme(axis.text=element_text(size=5),
              axis.text.x=element_text(angle=90, hjust=1))
ggsave('subseq_anisakidae_coi5p_pairwise.svg', width=10, height=10)




