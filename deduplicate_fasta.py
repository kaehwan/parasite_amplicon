from __future__ import print_function
import sys
import argparse
import time
import re
from Bio import SeqIO
from Bio.Seq import reverse_complement

# define functions
# to parse fasta file


def parsing_fasta(input_file):
    # read fasta entries and sequences
    entries, sequences = [], []
    records = SeqIO.parse(input_file, 'fasta')
    for record in records:
        entries += ['>%s' % str(record.description)]
        sequences += ['%s' % str(record.seq).replace('-', '')]
    return entries, sequences

# to write output fasta file


def write_out(output_file, starting_entries, starting_sequences, curated_sequences):
    print('\n-------\nfiltered sequences = %d from %d starting sequences\nresulting sequence = %d sequences\n-------\n' % (
        len(starting_entries) - len(curated_sequences), len(starting_entries), len(curated_sequences)))
    with open(output_file, 'w') as out:
        entry = starting_entries[starting_sequences.index(
            curated_sequences[0])]
        curated_entries = [entry]
        out.write('%s\n%s' % (entry, curated_sequences[0]))
        for i in range(1, len(curated_sequences)):
            entry = starting_entries[starting_sequences.index(
                curated_sequences[i])]
            curated_entries.append(entry)
            out.write('\n%s\n%s' % (entry, curated_sequences[i]))
    deleted = list(
        set([entry for entry in starting_entries if entry not in curated_entries]))
    if len(deleted) > 0:
        with open('removed_fasta_entries.txt', 'w') as out:
            out.write('%s\n' % deleted[0])
            for i in range(1, len(deleted)):
                out.write('\n%s\n' % deleted[i])

# to generate one k-mer


def kmer_gen(sequence, k, start=0):
    # sequence is a fasta sequence
    # k is the length of kmer
    # start is the position to start generating, which is 0
    kmer = sequence[start:(start + k + 1)]
    return kmer

# to cleanup empty sequence(s) in list after curation


def cleanup(_list):
    # _list is the list of sequences generated during curating
    while '' in _list:
        _list.remove('')
    return _list


def derep_longest(input_file, sequence_type='n'):
    starting_entries, starting_sequences = parsing_fasta(input_file)
    if len(starting_entries) < 1:
        sys.exit('\n\nInvalid input file\n\n')

    editing = list(set(starting_sequences))
    editing.sort(key=len)

    for i in range((len(editing) - 1)):
        if len(editing[i]) >= minimum:
            comparing = editing[i]
            kmer = kmer_gen(comparing, k=length, start=0)
            for seq in editing[(i + 1):]:
                if kmer in seq:
                    if comparing in seq:
                        editing[i] = ''
                        break
                elif sequence_type == 'n' and reverse_complement(kmer) in seq:
                    if reverse_complement(comparing) in seq:
                        editing = ''
                        break
        else:
            editing = ''
    editing = cleanup(editing)
    write_out(output_file, starting_entries, starting_sequences, editing)


parser = argparse.ArgumentParser(prog='\n\nSequence Database Dereplicator program',
                                 usage='\n%(prog)s : dereplicates nucleotide database'
                                 'from a list of sequences (by exact match).\n\n')
parser.add_argument('-mode', dest='mode', required=True, choices=['derep'],
                    help='dereplicate a fasta file with multiple entries')
parser.add_argument('-in', dest='input_file', type=argparse.FileType('r'), required=True,
                    help='Input file containing fasta entries')
parser.add_argument('-out', dest='output_file', required=True,
                    help='Your output file')
parser.add_argument('-n', dest='database', action='store_const', const='n',
                    help='nucleotide sequences')
parser.add_argument('-len', dest='length', default=30, type=int,
                    help='length of k-mer')
parser.add_argument('-min_length', dest='minimum', default=1, type=int,
                    help='minimum sequence length in your data (for dereplication)')
args = parser.parse_args()

# Parsing the argument !!!

mode = args.mode
input_file = args.input_file
output_file = args.output_file
database = args.database
length = args.length
minimum = args.minimum

if mode == 'derep':
    derep_longest(input_file, sequence_type=database)
    print('Author\t: %s\nDate\t: %s' % (
        'Sim Kae Hwan', 'April 11, 2023'
    ))
