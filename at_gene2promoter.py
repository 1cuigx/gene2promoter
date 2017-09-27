#!/usr/bin/env python3

'''
this script is to extract the location of a gene from gff3 file
then extract its promoter sequence from the genome fasta file
'''

import re


gene_loc = {}
with open('Arabidopsis_thaliana.TAIR10.37.gff3', 'r') as gff:
    for line in gff:
        if line[0] != '#':
            if line.split("\t")[2] == "gene":
                gene_id = re.search("AT[1-5MC]G\d+", line).group(0)
                gene_loc[gene_id] = [(int(line.split("\t")[3]), int(line.split("\t")[4])), line.split("\t")[6]]

promoter_loc = {}
for gene in gene_loc:
    if gene_loc[gene][1] == "+":
        promoter_loc[gene] = [min(gene_loc[gene][0]) - 2000, min(gene_loc[gene][0])]
    elif gene_loc[gene][1] == "-":
        promoter_loc[gene] = [max(gene_loc[gene][0]), max(gene_loc[gene][0]) + 2000]

chr_seq = {}
with open('Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', 'r') as fa:
    for line in fa:
        line = line.strip()
        if line:
            if line[0] == ">":
                chr_id = line[1:].split(" ")[0]
                chr_seq[chr_id] = ''
            else:
                chr_seq[chr_id] = line


def rev_com(dna):
    dna = dna[::-1]
    dna_rc = ''
    for a in dna:
        if a == 'A':
            dna_rc += 'T'
        elif a == 'T':
            dna_rc += 'A'
        elif a == 'G':
            dna_rc += 'C'
        else:
            dna_rc += 'G'
    return dna_rc

promoter_seq = {}
for id in promoter_loc:
    start = promoter_loc[id][0]
    end = promoter_loc[id][1]
    if gene_loc[id][1] == '+':
        if id[2] in ['1', '2', '3', '4', '5']:
            promoter_seq[id] = chr_seq[id[2]][start:end]
        elif id[2] == 'M':
            promoter_seq[id] = chr_seq['Mt'][start:end]
        elif id[2] == 'C':
            promoter_seq[id] = chr_seq['Pt'][start:end]
    else:
        if id[2] in ['1', '2', '3', '4', '5']:
            promoter_seq[id] = rev_com(chr_seq[id[2]][start:end])
        elif id[2] == 'M':
            promoter_seq[id] = rev_com(chr_seq['Mt'][start:end])
        elif id[2] == 'C':
            promoter_seq[id] = rev_com(chr_seq['Pt'][start:end])

with open('Arabidopsis_thaliana.TAIR10.promoter.fa', 'w') as prom:
    for t in promoter_seq:
        print('>' + t, promoter_seq[t], sep='\n', file=prom)
