#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
create only a coding bedfile from ensemble"
"""

import argparse
import csv
import gzip
import shlex


def main():
    args = parse_args()
    cds_bed = parse_gtf(args.gtf_file, args.organism_name)
    return cds_bed

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
    args = parser.parse_args()
    return args



def parse_gtf(gtf_file, organism_name):
    outfile = organism_name + "_cds.bed"
    with gzip.open(gtf_file) as csvfile, open(outfile, 'w') as fout:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            try:
                if row[2] == 'CDS':
                    chrom = row[0]
                    start = row[3]
                    stop = row[4]
                    transcript_id = row[8].split(';')[1]
                    exon_number = row[8].split(';')[2]
                    exon = shlex.split(exon_number.split(" ")[2])[0]
                    gene_name = row[8].split(';')[3]
                    gene = shlex.split(gene_name.split(" ")[2])[0]
                    transcript = shlex.split(transcript_id.split(" ")[2])[0]
                    out = (chrom, start, stop, gene, transcript, exon)
                    fout.write('\t'.join(str(i) for i in out) + '\n')
                else:
                    pass
            except IndexError:
                pass
    return outfile

                
if __name__ == "__main__":
    main()
