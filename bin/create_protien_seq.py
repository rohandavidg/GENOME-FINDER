#/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
create protien sequence
"""


from Bio import SeqIO
import sys
import pyfaidx
from pyfaidx import Fasta
import argparse
import collections
import pprint
import os
import csv
from collections import defaultdict
import logging
import time
import datetime
from Bio import Seq
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
from collections import Counter
from do_logging import configure_logger
import parser_gtf

log_filename = "create_protien_seq"
test1 = "M306I"
gene1 = 'embB'

def main():
   logger = configure_logger(log_filename)
   args = parse_args()
   bed_file = parser_gtf.main(args.gtf_file, args.organism_name)
   genes = Fasta(args.fasta)
   bed_dict = parse_bed_file(bed_file, logger)
   seq_dict = faidx_the_lot(bed_dict, genes, logger)
   amino_dict = gene_level_aa(seq_dict, gene1, logger)
#   genomic_position = get_renomic_postion(amino_dict, test1, gene1)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', dest='organism_name',
                        help='name of the organism',
                        required=True)
    parser.add_argument('-g', dest='gtf_file',
                        help="ensembl gtf file",
                        required=True)
    parser.add_argument('-f', dest='fasta', required=True,
                        help="fasta sequences of the organism")
    args = parser.parse_args()
    return args



def parse_bed_file(bed_file, logger):
    bed_transcipt_dict = defaultdict(list)
    with open(bed_file) as bin:
       some_bin = bin.readlines()
       for i, raw_line in enumerate(some_bin):
          good_line = raw_line.strip().split("\t")
          check_space = [bool(s.isspace()) or " " in s for s in good_line]
          for l in check_space:
             if l == True:
                logger.debug("space found on line {0}".format(i))
             else:
                pass
             chrom = good_line[0]
             if not chrom.startswith('chr'):
                logger.debug("{0} should start with chr on line {1}".format(chrom, i))
             else:
                pass
             start = good_line[1]
             stop = good_line[2]
             assert start.isdigit()
             assert stop.isdigit()
             if int(start) < int(stop):
                pass
             else:
                logger.debug("start:{0} is greater than stop {1} on line {3}".format(start, stop, i))
             try:
                gene_transcript = good_line[3]
                gene = gene_transcript.split(":")[0]
                chrom_start_stop = [chrom, start, stop, gene]
                bed_transcipt_dict[gene].append(chrom_start_stop)
             except IndexError:
                logger.debug("{0} gene not present in the 3rd column".format(gene_transcript))
    return bed_transcipt_dict
    

def faidx_the_lot(bed_transcipt_dict, genes, logger):
    bed_seq_dict = {}
    for key, value in bed_transcipt_dict.items():
       gene = key.split(":")[0]
       for i in value:
          chrom = i[0]
          start = int(i[1])
          stop = int(i[2])
          gene = i[3]
          try:
             sequence = genes[chrom][start:stop].seq
             length = int(len(sequence)/3)
             new_value  = i + [sequence, length]
             bed_seq_dict[key] = new_value
          except ValueError:
             logger.debug('Missing fields in target.bed file. please check if target.bed is in std format')
    return bed_seq_dict


def get_amino_acid(sequence):
    seq = Seq(sequence)
    table = 1
    min_pro_len = 20
    for strand, nuc in [(1, seq)]:
        for frame in range(3):
           yield nuc[frame:].translate(table=11)


def generate_aa(gen_sequence):
    for aa in gen_sequence:
       return aa
#        if aa.startswith('M'):
#            return aa
#        else:
#            pass


def gene_level_aa(bed_seq_dict, logger, gene1):
    gene_seq_dict = defaultdict(list)
    exon_seq = ''
    for key, value in bed_seq_dict.items():
       if key == 'embB':
          print generate_aa(get_amino_acid(value[4]))
#             try:
#                exon_seq += "".join(str(j) for j in i[-2:-1])
#                print exon_seq
#             except IndexError:
#                logger.debug('Missing strand information in {0}'.format(i))
#    amino_acid = get_amino_acid(exon_seq)
#    print generate_aa(amino_acid)
#    print new_value
#       gene_seq_dict[key].append(new_value)
#    return gene_seq_dict


def get_genomic_postion(gene_seq_dict, test1, gene):
   for k, v in gene_seq_dict.items():
      print k


if __name__ == "__main__":
   main()
