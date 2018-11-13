"""
Run this after prefilter_annotation_file.py
This script will remove contigs/chromosomes from a genome which contain
no transcripts in a given reference file.
It is critical that this be run after filtering filtering your reference,
or there is a chance you will encounter unnecessary runtime
or errors in crossmap.
"""

__author__='james'
import warnings
import argparse
import inspect
from Bio import SeqIO
from plastid.util.io.openers import get_short_name, argsopener
from plastid.util.scriptlib.argparsers import AnnotationParser
from plastid.util.io.filters import NameDateWriter, AbstractReader
from plastid.util.services.exceptions import FileFormatWarning

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


def main():

    ap = AnnotationParser()
    annotation_file_parser = ap.get_parser()
    parser = argparse.ArgumentParser(parents=[annotation_file_parser])
    parser.add_argument("genome", metavar='genome', type=str,
                        default = None, help="File handle for genome to filter")
    parser.add_argument("outbase", metavar="outbase", type=str,
                 help="Basename for output file")
    args = parser.parse_args()

    printer.write('Reading reference file...')
    source = ap.get_transcripts_from_args(args,printer=printer)  # Iterator w/ transcript info
    chromosomes_with_transcripts = set()

    for tx in source:
        chromosomes_with_transcripts.add(tx.spanning_segment.chrom)

    printer.write('Filtering genome...')
    unfiltered_genome_iterator = SeqIO.parse(args.genome, 'fasta')
    filtered_genome_iterator = (record for record in unfiltered_genome_iterator \
                      if record.name in chromosomes_with_transcripts)

    printer.write('Writing genome...')
    SeqIO.write(filtered_genome_iterator, args.outbase + '_filtered.fa', 'fasta')



if __name__ == '__main__':
    # main()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FileFormatWarning)
        main()