#! /usr/bin/env python3

from collections import defaultdict

d1 = {
'Aenictus_a' : 'ATTATGCGCACGTTCG',
'Aenictus_b' : 'ATTAAGCGCTCGTTCG',
'Aenictus_c' : 'ATTTTGCGCTCGTTCG',
'Dorylus_z'  : 'TATATGCGGTATCG--'
}

d2 = {
'Aenictus_z' : '?????GCGCGCGGGCG',
'Aenictus_b' : 'GCGCGGCGCGCGGGCG',
'Aenictus_c' : 'GCGCGGCGCGCGGGCG',
'Dorylus_z'  : '--GCGCGCCGCCGGGC'
}

d3 ={
'Aenictus_a' : 'ATATATATATAT',
'Aenictus_b' : 'ATTAAATATTAT',
'Aenictus_x' : 'ATTTTAATTATT'
}

in_files = [d1,d2,d3] 

# create empty dictionary of lists
concatenated = defaultdict(list)

# first create list of taxa in all alignments
# you need this to insert empty seqs in
# the concatenated alignment
all_taxa = []

for alignment in in_files:
    for taxon in alignment.keys():
        if taxon not in all_taxa:
            all_taxa.append(taxon) 
print(all_taxa)

for alignment in in_files:
    # get empty sequence if there is missing taxon
    # getting length from first element of list of keys
    # created from the original dict for this alignment
    empty_seq = '?' * len(alignment[list(alignment.keys())[0]])

    for taxon in all_taxa:
        if taxon not in alignment.keys():
            concatenated[taxon].append(empty_seq)
        else:
            concatenated[taxon].append(alignment[taxon])

for taxon, seqs in concatenated.items():

    seqs = ''.join(seqs)
    concatenated[taxon] = seqs    

print(concatenated)
