from primer3 import bindings
from pprint import pprint
from Bio import SeqIO
from Bio.Seq import reverse_complement
import argparse, operator, re, os.path, collections, traceback, pickle
parser = argparse.ArgumentParser(description="Make fusion RT-qPCR primers for all genes with one of the primers \
overlapping an exon-exon junction if the gene has multiple exons. Requires primer3_core, primer3-py, and biopython")
parser.add_argument("gff_file")
parser.add_argument("fasta_file")
args = parser.parse_args()
gff = [line.strip().split('\t') for line in open(args.gff_file,'r') if len(line.strip().split('\t')) > 1]
gff_raw = open(args.gff_file,'r').read()
gene_count = gff_raw.count('gene\t')
merged_cdna = ''
class gene_loc():
    """Class that records location information and sequence for gene"""
    def __init__(self, name, chrom, min_coord, max_coord, strand, introns, exons):
        self.name = name
        self.chrom = chrom
        self.min_coord = min_coord
        self.max_coord = max_coord
        self.strand = strand
        self.introns = introns
        self.exons = exons
        self.seq = str(fasta_dict[self.chrom].seq[self.min_coord - 1: self.max_coord])
        def get_cdna_and_exon_locs(dna, exons):
            #Takes dna sequence and list of exon coordinates and returns tuple of cdna and list of exon locations within the cdna
            global merged_cdna
            normalizer = exons[0] - 1
            exons = [exon - normalizer for exon in exons] #shifts exon coordinates so they are relative to gene instead of contig
            cdna = ''
            exon_juncs = []
            intron_sum = 0
            for i in range(0,len(exons),2):
                cdna += dna[exons[i]-1:exons[i+1]]
                if i > 0:
                    intron_sum += exons[i] - exons[i-1]
                exon_juncs.append(exons[i+1] - intron_sum)
            exon_juncs.pop()
            merged_cdna += 'n' + cdna + 'n' + reverse_complement(cdna)
            return (cdna, exon_juncs)
        if len(self.exons) > 2:
            self.cdna, self.exon_juncs = get_cdna_and_exon_locs(self.seq, self.exons)
        else:
            self.cdna = self.seq
            self.exon_juncs = []
def get_cds_coord(gff):
    #Get coordinates for gff
    coord = []
    temp = []
    for idx, line in enumerate(gff):
        if 'gene' in line[2]:
            if temp:
                coord.append([gff[gene_idx][8][3:gff[gene_idx][8].index(';')], gff[gene_idx][0], temp, gff[gene_idx][6]])
            gene_idx = idx
            temp = []
        if 'CDS' in line[2]:
            if gene_count -1 == len(coord):
                coord.append([gff[gene_idx][8][3:gff[gene_idx][8].index(';')], gff[gene_idx][0], temp, gff[gene_idx][6]])
            temp.append(line[3])
            temp.append(line[4])
    coord_sorted = []
    for values in coord:
        exons = sorted(map(int, values[2]))
        introns = list(zip(exons[1::2], exons[2::2]))
        coord_sorted.append(gene_loc(name=values[0], chrom=values[1], min_coord=min(exons), max_coord=max(exons), strand=values[3], introns=introns, exons=exons))
    with open('coord','w') as myfile:
        pprint(coord, myfile)
    return coord_sorted
fasta_dict = SeqIO.to_dict(SeqIO.parse(open(args.fasta_file,'r'), 'fasta'))
merged_contigs = ''
for fasta in fasta_dict:
    merged_contigs += 'n' + str(fasta_dict[fasta].seq) + 'n' + str(fasta_dict[fasta].reverse_complement())
coord = get_cds_coord(gff)
genes_with_errors = []
num_fail = 0
overlap_map = {'+': 'PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION', '-': 'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION'}
f_primer_map = {'+': 'F', '-': 'R'}
r_primer_map = {'+': 'R', '-': 'F'}
with open('Af3357_qPCR_primers.csv','w') as f, \
open('genes_with_errors','w') as genes_with_errors_f, open('primers.gff','w') as gff_f:
    temp = []
    for idx, gene in enumerate(coord,1):
        if len(gene.seq) < 100:
            continue
        if len(gene.exons) > 2:
            print(gene.name)
            goi = bindings.designPrimers({
                'SEQUENCE_ID': gene.name,
                'SEQUENCE_TEMPLATE': gene.cdna,
                'SEQUENCE_OVERLAP_JUNCTION_LIST': gene.exon_juncs,
                },
                {'PRIMER_PRODUCT_SIZE_RANGE': [[75, 150]],
                'PRIMER_EXPLAIN_FLAG':1,
                'PRIMER_MAX_TM':68,
                'PRIMER_MIN_TM':52,
                'PRIMER_PICK_INTERNAL_OLIGO':0,
                overlap_map[gene.strand]: 9})
        else:
            print(gene.name)
            goi = bindings.designPrimers({
                'SEQUENCE_ID': gene.name,
                'SEQUENCE_TEMPLATE': gene.cdna
                },
                {'PRIMER_PRODUCT_SIZE_RANGE': [[75, 150]],
                'PRIMER_EXPLAIN_FLAG':1,
                'PRIMER_MAX_TM':68,
                'PRIMER_MIN_TM':52,
                'PRIMER_PICK_INTERNAL_OLIGO':0})
        if 'PRIMER_LEFT_0' not in goi:
            print('LEFT:' + goi['PRIMER_LEFT_EXPLAIN'])
            genes_with_errors_f.write(gene.name + ' ' +  goi['PRIMER_LEFT_EXPLAIN'])
            continue
        if 'PRIMER_RIGHT_0' not in goi:
            print('RIGHT:' + goi['PRIMER_RIGHT_EXPLAIN'])
            genes_with_errors_f.write(gene.name + ' ' +  goi['PRIMER_RIGHT_EXPLAIN'])
            continue
        temp.append(','.join([gene.name+'-q' + f_primer_map[gene.strand], goi['PRIMER_LEFT_0_SEQUENCE'], 
            str(goi['PRIMER_LEFT_0_PENALTY']), str(goi['PRIMER_PAIR_0_PRODUCT_SIZE']), str(goi['PRIMER_LEFT_0_SEQUENCE'] not in gene.seq),
            str(merged_cdna.count(goi['PRIMER_LEFT_0_SEQUENCE'][-12:])), str(merged_contigs.count(goi['PRIMER_LEFT_0_SEQUENCE'][-12:]))])+'\n')
        temp.append(','.join([gene.name+'-q'+ r_primer_map[gene.strand], goi['PRIMER_RIGHT_0_SEQUENCE'], 
            str(goi['PRIMER_RIGHT_0_PENALTY']), str(goi['PRIMER_PAIR_0_PRODUCT_SIZE']), str(goi['PRIMER_RIGHT_0_SEQUENCE'] not in gene.seq),
            str(merged_cdna.count(goi['PRIMER_RIGHT_0_SEQUENCE'][-12:])), str(merged_contigs.count(goi['PRIMER_RIGHT_0_SEQUENCE'][-12:]))])+'\n')
    temp.sort()
    f.write(','.join(['primer name', 'primer sequence', 'penalty', 'product size (bp)', 'split_primer', '# of times 3\' 12bp of primer is present in cDNA', '# of times 3\' 12bp of primer is present in gDNA'])+'\n')
    for line in temp:
        f.write(line)
