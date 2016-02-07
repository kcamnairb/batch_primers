from primer3 import bindings
from pprint import pprint
from Bio import SeqIO
import argparse, operator, re, os.path, collections, traceback
parser = argparse.ArgumentParser(description="Make fusion PCR primers for making KOs of all genes.")
parser.add_argument("gff_file")
parser.add_argument("fasta_file")
args = parser.parse_args()
gff = [line.strip().split('\t') for line in open(args.gff_file,'r') if len(line.strip().split('\t')) > 1]
gff_raw = open(args.gff_file,'r').read()
gene_count = gff_raw.count('gene\t')
genome_name = os.path.splitext(os.path.basename(args.fasta_file))[0]
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
    #Keeps only lowest and highest CDS coordinates
    coord_sorted = []
    gene_loc = collections.namedtuple('gene_loc', ['name', 'chrom', 'min_coord', 'max_coord', 'strand', 'introns', 'exons'])
    for values in coord:
        exons = sorted(map(int, values[2]))
        introns = list(zip(exons[1::2], exons[2::2]))
        coord_sorted.append(gene_loc(values[0], values[1], min(exons), max(exons), values[3], introns, exons))
    
    
    with open('coord','w') as myfile:
        pprint(coord, myfile)
    #with open('gff','w') as myfile:
    #    for line in gff:
    #        myfile.write(str(line) + '\n')
    return coord_sorted
fasta_dict = SeqIO.to_dict(SeqIO.parse(open(args.fasta_file,'r'), 'fasta'))
coord = get_cds_coord(gff)
pyrg5R = 'gcaataagcccaaccctatcggc'
pyrg3F = 'tttgtaccggagtgtctgaagg'
ptra5R = 'gggatcccgtaatcaattgccc'
ptra3F = 'caagagcggctcatcgtcaccc'
forward_primer, reverse_primer = 'SO:0000121', 'SO:0000132'
genes_with_errors = []

with open(genome_name+'_KO_primers_pyrG_selection.csv','w') as pyrg_f, open(genome_name+'_KO_primers_ptrA_selection.csv','w') as ptra_f, \
open('genes_with_errors','w') as genes_with_errors_f, open('primers.gff','w') as gff_f:
    files = (pyrg_f, ptra_f)
    for f in files:
        f.write(','.join(['primer name','primer sequence','penalty', 'product size'])+'\n')
    for gene in coord[:10]:

        try:
            print(gene.name)
            if gene.min_coord - 1400  > 0 and gene.max_coord + 1400 < len(fasta_dict[gene.chrom]) and gene.max_coord - gene.min_coord > 200:
                left = bindings.designPrimers({
                    'SEQUENCE_ID': gene.name,
                    'SEQUENCE_TEMPLATE': str(fasta_dict[gene.chrom].seq),
                    'SEQUENCE_INCLUDED_REGION': [gene.min_coord-1400, 1500],
                    'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST':[-1,-1,gene.min_coord, 100],
                    },
                    {'PRIMER_PRODUCT_SIZE_RANGE': [[1100,1300]]})
                right = bindings.designPrimers({
                    'SEQUENCE_ID': gene.name,
                    'SEQUENCE_TEMPLATE': str(fasta_dict[gene.chrom].seq),
                    'SEQUENCE_INCLUDED_REGION': [gene.max_coord-100, 1500],
                    'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST':[gene.max_coord-100, 100,-1,-1],
                    },
                    {'PRIMER_PRODUCT_SIZE_RANGE': [[1100,1300]]})
                left_to_right_dist = (right['PRIMER_RIGHT_0'][0] - right['PRIMER_RIGHT_0'][1]) - (left['PRIMER_LEFT_0'][0] + left['PRIMER_LEFT_0'][1]) 
                nested_range = [left_to_right_dist-200, left_to_right_dist]
                nest = bindings.designPrimers({
                    'SEQUENCE_ID': gene.name,
                    'SEQUENCE_TEMPLATE': str(fasta_dict[gene.chrom].seq),
                    'SEQUENCE_INCLUDED_REGION': [sum(left['PRIMER_LEFT_0']), left_to_right_dist]
                    },
                    {'PRIMER_PRODUCT_SIZE_RANGE' : [nested_range]})
                (nest['PRIMER_RIGHT_0_SEQUENCE'], nest['PRIMER_LEFT_0_SEQUENCE'])
                #introns_in_deletion = [intron for intron in gene.introns if intron[0] > left['PRIMER_RIGHT_0'][0] and intron[1] < right['PRIMER_LEFT_0'][0]]
                #print(gene.name, introns_in_deletion)
                #if introns_in_deletion:
                #    pprint(introns_in_deletion)
                #    qpcr = bindings.designPrimers({
                #        'SEQUENCE_ID': gene.name,
                #        'SEQUENCE_TEMPLATE': str(fasta_dict[gene.chrom].seq),
                #        'SEQUENCE_INCLUDED_REGION': [left['PRIMER_RIGHT_0'][0], right['PRIMER_LEFT_0'][0]-left['PRIMER_RIGHT_0'][0]],
                #        'SEQUENCE_TARGET' : [[x[0],x[1]-x[0]] for x in gene.introns],
                #        'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST':[[x[0]-100, 63, x[1]+37, 63] for x in introns_in_deletion],
                #        },
                #        {'PRIMER_PRODUCT_SIZE_RANGE': [[74,max([x[1]-x[0] for x in introns_in_deletion])+75]]})
                #else:
                #    qpcr = bindings.designPrimers({
                #        'SEQUENCE_ID': gene.name,
                #        'SEQUENCE_TEMPLATE': str(fasta_dict[gene.chrom].seq),
                #        'SEQUENCE_INCLUDED_REGION': [gene.min_coord, gene.max_coord-gene.min_coord],
                #        'SEQUENCE_EXCLUDED_REGION': gene.introns
                #        },
                #        {'PRIMER_PRODUCT_SIZE_RANGE': [75, 150]})
                #(qpcr['PRIMER_RIGHT_0_SEQUENCE'], qpcr['PRIMER_LEFT_0_SEQUENCE'])
                #pprint(left)
                #pprint(right)
                #pprint(nest)
            else: continue
        except KeyError as e:
            traceback.print_exc()
            genes_with_errors.append(gene.name)
            pprint(gene.name, genes_with_errors_f)
            continue
        else:
            if gene.strand == '+':
                gff_f.write('\t'.join([gene.chrom,'primer3', forward_primer,str(left['PRIMER_LEFT_0'][0]+1),str(left['PRIMER_LEFT_0'][0]+left['PRIMER_LEFT_0'][1]),str(left['PRIMER_LEFT_0_PENALTY']),'+','.','Name='+gene.name+'-5F'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', reverse_primer,str(left['PRIMER_RIGHT_0'][0]+2-left['PRIMER_RIGHT_0'][1]),str(left['PRIMER_RIGHT_0'][0]+1),str(left['PRIMER_RIGHT_0_PENALTY']),'-','.','Name='+gene.name+'-5R'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', forward_primer,str(right['PRIMER_LEFT_0'][0]+1),str(right['PRIMER_LEFT_0'][0]+right['PRIMER_LEFT_0'][1]),str(right['PRIMER_LEFT_0_PENALTY']),'+','.','Name='+gene.name+'-3F'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', reverse_primer,str(right['PRIMER_RIGHT_0'][0]+2-right['PRIMER_RIGHT_0'][1]),str(right['PRIMER_RIGHT_0'][0]+1),str(right['PRIMER_RIGHT_0_PENALTY']),'-','.','Name='+gene.name+'-3R'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', forward_primer,str(nest['PRIMER_LEFT_0'][0]+1),str(nest['PRIMER_LEFT_0'][0]+nest['PRIMER_LEFT_0'][1]),str(nest['PRIMER_LEFT_0_PENALTY']),'+','.','Name='+gene.name+'-nestF'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', reverse_primer,str(nest['PRIMER_RIGHT_0'][0]+2-nest['PRIMER_RIGHT_0'][1]),str(nest['PRIMER_RIGHT_0'][0]+1),str(nest['PRIMER_RIGHT_0_PENALTY']),'-','.','Name='+gene.name+'-nestR'])+'\n')
            elif gene.strand == '-':
                gff_f.write('\t'.join([gene.chrom,'primer3', forward_primer,str(right['PRIMER_LEFT_0'][0]+1),str(right['PRIMER_LEFT_0'][0]+right['PRIMER_LEFT_0'][1]),str(right['PRIMER_LEFT_0_PENALTY']),'+','.','Name='+gene.name+'-5R'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', reverse_primer,str(right['PRIMER_RIGHT_0'][0]+2-right['PRIMER_RIGHT_0'][1]),str(right['PRIMER_RIGHT_0'][0]+1),str(right['PRIMER_RIGHT_0_PENALTY']),'-','.','Name='+gene.name+'-5F'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', forward_primer,str(left['PRIMER_LEFT_0'][0]+1),str(left['PRIMER_LEFT_0'][0]+left['PRIMER_LEFT_0'][1]),str(left['PRIMER_LEFT_0_PENALTY']),'+','.','Name='+gene.name+'-3F'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', reverse_primer,str(left['PRIMER_RIGHT_0'][0]+2-left['PRIMER_RIGHT_0'][1]),str(left['PRIMER_RIGHT_0'][0]+1),str(left['PRIMER_RIGHT_0_PENALTY']),'-','.','Name='+gene.name+'-3F'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', forward_primer,str(nest['PRIMER_LEFT_0'][0]+1),str(nest['PRIMER_LEFT_0'][0]+nest['PRIMER_LEFT_0'][1]),str(nest['PRIMER_LEFT_0_PENALTY']),'+','.','Name='+gene.name+'-nestR'])+'\n')
                gff_f.write('\t'.join([gene.chrom,'primer3', reverse_primer,str(nest['PRIMER_RIGHT_0'][0]+2-nest['PRIMER_RIGHT_0'][1]),str(nest['PRIMER_RIGHT_0'][0]+1),str(nest['PRIMER_RIGHT_0_PENALTY']),'-','.','Name='+gene.name+'-nestF'])+'\n')
            for f in files:
                if f.name == pyrg_f.name:
                    overlap5R = pyrg5R
                    overlap3F = pyrg3F
                    insert_size = 1500
                elif f.name == ptra_f.name:
                    overlap5R = ptra5R
                    overlap3F = ptra3F
                    insert_size = 2000
                if gene.strand == '+':
                    f.write(','.join([gene.name+'-5F', left['PRIMER_LEFT_0_SEQUENCE'], str(left['PRIMER_LEFT_0_PENALTY']), str(left['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-5R', overlap5R + left['PRIMER_RIGHT_0_SEQUENCE'], str(left['PRIMER_RIGHT_0_PENALTY']), str(left['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-3F', overlap3F + right['PRIMER_LEFT_0_SEQUENCE'], str(right['PRIMER_LEFT_0_PENALTY']), str(right['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-3R', right['PRIMER_RIGHT_0_SEQUENCE'], str(right['PRIMER_RIGHT_0_PENALTY']), str(right['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-nestF', nest['PRIMER_LEFT_0_SEQUENCE'], str(nest['PRIMER_LEFT_0_PENALTY']), str(nest['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-nestR', nest['PRIMER_RIGHT_0_SEQUENCE'], str(nest['PRIMER_RIGHT_0_PENALTY']), str(nest['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-5F-3R KO', '', '', str(left['PRIMER_PAIR_0_PRODUCT_SIZE'] + insert_size + right['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-5F-3R WT', '', '', str(right['PRIMER_RIGHT_0'][0] -  left['PRIMER_LEFT_0'][0])])+'\n')
                    #f.write(','.join([gene.name+'-qF', qpcr['PRIMER_LEFT_0_SEQUENCE'], str(qpcr['PRIMER_LEFT_0_PENALTY']), str(qpcr['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    #f.write(','.join([gene.name+'-qR', qpcr['PRIMER_RIGHT_0_SEQUENCE'], str(qpcr['PRIMER_RIGHT_0_PENALTY']), str(qpcr['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    
                if gene.strand == '-':
                    f.write(','.join([gene.name+'-5F', right['PRIMER_RIGHT_0_SEQUENCE'], str(right['PRIMER_RIGHT_0_PENALTY']), str(right['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-5R', overlap5R + right['PRIMER_LEFT_0_SEQUENCE'], str(right['PRIMER_LEFT_0_PENALTY']), str(right['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-3F', overlap3F + left['PRIMER_RIGHT_0_SEQUENCE'], str(left['PRIMER_RIGHT_0_PENALTY']), str(left['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-3R', left['PRIMER_LEFT_0_SEQUENCE'], str(left['PRIMER_LEFT_0_PENALTY']), str(left['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-nestF', nest['PRIMER_RIGHT_0_SEQUENCE'], str(nest['PRIMER_RIGHT_0_PENALTY']), str(nest['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-nestR', nest['PRIMER_LEFT_0_SEQUENCE'], str(nest['PRIMER_LEFT_0_PENALTY']), str(nest['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-5F-3R KO', '', '', str(left['PRIMER_PAIR_0_PRODUCT_SIZE'] + insert_size + right['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    f.write(','.join([gene.name+'-5F-3R WT', '', '', str(right['PRIMER_RIGHT_0'][0] -  left['PRIMER_LEFT_0'][0])])+'\n')
                    #f.write(','.join([gene.name+'-qF', qpcr['PRIMER_RIGHT_0_SEQUENCE'], str(qpcr['PRIMER_RIGHT_0_PENALTY']), str(qpcr['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')
                    #f.write(','.join([gene.name+'-qR', qpcr['PRIMER_LEFT_0_SEQUENCE'], str(qpcr['PRIMER_LEFT_0_PENALTY']), str(qpcr['PRIMER_PAIR_0_PRODUCT_SIZE'])])+'\n')

