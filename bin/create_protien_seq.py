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


def main():
   args = parse_args()
   genes = Fasta(args.fasta)


def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', dest='organism_name',
                       help='name of the organism',
                       required=True)
    parser.add_argument('-b', dest='bed_file',
                        help="bed file of the organism cds region",
                        required=True)
    parser.add_argument('-f', dest='fasta',
                        help="fasta sequences of the organism")
    args = parser.parse_args()
    return args


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
