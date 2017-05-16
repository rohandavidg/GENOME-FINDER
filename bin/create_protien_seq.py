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

def main():
   logger = configure_logger()
   args = parse_args()
   bed_file = parser_gtf(args.organism_name, args.gtf_file)
   print bed_file
   genes = Fasta(args.fasta)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', dest='organism_name',
                       help='name of the organism',
                       required=True)
#    parser.add_argument('-b', dest='bed_file',
#                        help="bed file of the organism cds region",
#                        required=True)
    parser.add_argument('-g', dest='gtf_file',
                        help="ensembl gtf file",
                        required=True)
    parser.add_argument('-f', dest='fasta',
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
                chrom_start_stop = [chrom, start, stop]
                bed_transcipt_dict[gene].append(chrom_start_stop)
             except IndexError:
                logger.debug("{0} gene not present in the 3rd column".format(gene_transcript))
    return bed_transcipt_dict
    

def faidx_the_lot(bed_transcipt_dict, logger):
    bed_amino_dict = defaultdict(list)
    for key, value in bed_transcipt_dict.items():
        try:
            transcript = key.split(":")[1]
            for i in value:
                chrom = i[0]
                start = int(i[1])
                stop = int(i[2])
                strand = i[3]
                position = i[4]
                region = i[5]
                try:
                    sequence = genes[chrom][start:stop].seq
                    length = int(len(sequence)/3)
                    new_value  = i + [sequence, length]
                    bed_amino_dict[key].append(new_value)
                except ValueError:
                    logger.debug('Missing fields in target.bed file. please check if target.bed is in std format')
        except IndexError:
            sys.exit('Missing transcript in bed, check validate_bed_<date>.log for errors')
    return bed_amino_dict


if __name__ == "__main__":
   main()
