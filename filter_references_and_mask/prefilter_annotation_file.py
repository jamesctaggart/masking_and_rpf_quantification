"""
Script based on framework of cs generate (jdunn 2011)
Used to pre-process genome annotation files for use in crossmap and read alignment
Sorts features in annotation file and removes overlapping transcripts
"""

import warnings
import argparse
import subprocess
import inspect
import itertools
import gc
from plastid.util.io.openers import get_short_name, argsopener
from plastid.util.scriptlib.argparsers import AnnotationParser
from plastid.util.io.filters import NameDateWriter, AbstractReader
from plastid.util.services.exceptions import FileFormatWarning

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


def sort_reference_file(annotation_fh,annotation_format,printer):
    if annotation_format == 'GFF3':
        cmd = 'cat '+annotation_fh+' | grep -v "#" | sort -k1,1 -k4,4n >'+annotation_fh+'_sorted'
        subprocess.call(cmd,shell=True)
    elif annotation_format == 'GTF2':
        cmd = 'cat ' + annotation_fh + ' | grep -v "#" | sort -k1,1 -k4,4n >' + annotation_fh + '_sorted'
        subprocess.call(cmd,shell=True)
    elif annotation_format == 'BED':
        cmd = 'cat ' + annotation_fh + ' | grep -v "#" | sort -k1,1 -k2,2n >' + annotation_fh + '_sorted'
        subprocess.call(cmd,shell=True)
    else:
        printer.write('File format not recognized. File not sorted!')

def APPRIS_filter(transcript_dict, APPRIS_file):
    """
    This function will take dictionary of transcripts and filter them based
    on their presence in the APPRIS annotation as a primary/alternative transcript.
    This should get rid of transcripts annotated as minor.
    :param transcript_dict: keys: transcript name
    :param APPRIS_file:
    :return: Filtered transcript dictionary
    """
    # First, generate a set of known transcripts.
    APPRIS_transcripts = set()
    with open(APPRIS_file, 'r') as f:
        for line in f:
            fields = line.split()
            # gene_name = fields[0]
            tx_name = fields[1]
            # tag = fields[2]
            APPRIS_transcripts.add(tx_name)

    # Next, check each transcript name in dictionary against this set, returning those which match
    filtered_transcripts = {}
    count_removed = 0
    for tx in transcript_dict:
        if tx in APPRIS_transcripts:
            filtered_transcripts[tx] = transcript_dict[tx]
        else:
            count_removed+=1
    report_string = "%s of %s transcripts removed by APPRIS filter." % (str(count_removed),str(len(transcript_dict)))
    return filtered_transcripts, report_string


def attribute_filter(transcript_dict, attr_pair):
    """
    Remove transcripts who have an attribute with a particular value.
    :param transcript_dict:
    :param attr_tuple: Format (attribute name, attribute value)
    :return: filtered transcript dictionary
    """
    if len(attr_pair)>2:
        raise RuntimeError("Multi-attribute filtering not yet implemented.")
    count_removed = 0
    filtered_transcripts = {}
    for tx in transcript_dict:
        if transcript_dict[tx].attr[attr_pair[0]] == attr_pair[1]:
            count_removed += 1
        else:
            filtered_transcripts[tx] = transcript_dict[tx]
    report_string = "%s of %s transcripts removed by attribute filter." % (str(count_removed), str(len(transcript_dict)))
    return filtered_transcripts, report_string


def find_overlaps(transcript_dict):
    has_overlap = set()
    for tx1 in transcript_dict:
        cds1 = transcript_dict[tx1].get_cds()
        parent_gene1 = transcript_dict[tx1].get_gene()
        for tx2 in transcript_dict:
            cds2 = transcript_dict[tx2].get_cds()
            parent_gene2 = transcript_dict[tx2].get_gene()
            if tx1 != tx2:
                if cds1.overlaps(cds2) and parent_gene1 != parent_gene2:
                    has_overlap.add(parent_gene1)
                    has_overlap.add(parent_gene2)
                else:
                    continue
    return has_overlap


def main():

    ap = AnnotationParser()
    annotation_file_parser = ap.get_parser()
    parser = argparse.ArgumentParser(parents=[annotation_file_parser])
    parser.add_argument("outbase", metavar="outbase", type=str,
                 help="Basename for output file")
    parser.add_argument("--appris", metavar='APPRIS', type=str,
                        default = None, help="File handle for primary/alternative transcript reference")
    parser.add_argument("--remove_overlaps", action='store_true', default = False,
                        help = 'If true, remove overlapping transcripts from reference file.')
    parser.add_argument("--attr_filter", metavar='attrfilter', type=str, nargs = '+',
                        default = None, help="If specified, filter based on attribute specified")

    args = parser.parse_args()

    printer.write('Reading reference file...')
    source = ap.get_transcripts_from_args(args,printer=printer)  # Iterator w/ transcript info

    accepted_transcripts = {}
    transcripts = {}
    # code here adapted from cs generate
    # loop conditions

    last_chrom = None
    do_loop = True
    if args.remove_overlaps:
        removed_genes = open('%s_removed_genes' % args.outbase,'w')  # record genes removed in run

    # to save memory, we process one chromosome at a time if input file is sorted
    # knowing that at that moment all transcript parts are assembled
    counter = 0
    while do_loop == True:
        try:
            tx = next(source)
            counter+=1
            if counter%5000 == 0:
                printer.write('Parsed %s features.' % str(counter))
        except StopIteration:
            do_loop = False

        try:
            # enter here if chromosome is completely processed or EOF
            if tx.spanning_segment.chrom != last_chrom or do_loop == False:
                if last_chrom is not None or do_loop == False:
                    printer.write("Removing minor and overlapping transcripts on chromosome/contig '%s'" % last_chrom)

                    # Remove minor transcripts
                    if args.appris:
                        transcripts, report = APPRIS_filter(transcripts,args.appris)
                        printer.write(report)
                    if args.attr_filter:
                        transcripts, report = attribute_filter(transcripts, args.attr_filter)
                        printer.write(report)

                    if args.remove_overlaps:
                        overlapped_genes = find_overlaps(transcripts)  # Determine genes on chromosome which have overlapping transcripts
                        for gene in overlapped_genes:  # Record which genes are removed
                            removed_genes.write(gene+'\n')

                        to_remove = set()
                        for tx_handle in transcripts:  # Remove any transcripts which are from overlapped genes
                            if transcripts[tx_handle].get_gene() in overlapped_genes:
                                to_remove.add(tx_handle)
                        for tx_handle in to_remove:
                            transcripts.pop(tx_handle, None)

                    accepted_transcripts.update(transcripts)

                    del transcripts 
                    gc.collect()
                    del gc.garbage[:]
                    transcripts = {}


                # reset last chrom
                last_chrom = tx.spanning_segment.chrom

            # otherwise, remember transcript
            else:
                transcripts[tx.get_name()] = tx

        # exit gracefully if no transcripts found
        except UnboundLocalError:
            print "UNBOUND LOCAL ERROR"
            pass

    if args.remove_overlaps:
        removed_genes.close()

    printer.write("Writing output ...")
    printer.write("ALERT: Reference files must be manually sorted once this is done!")
    printer.write("See: http://plastid.readthedocs.io/en/latest/examples/a1_genome_setup.html")
    if args.annotation_format == 'GTF2':
        outf = argsopener('%s.gtf' % args.outbase, args)
        for tx in accepted_transcripts:
            outf.write(accepted_transcripts[tx].as_gtf())
        outf.close()
        # printer.write('Sorting output...')
        # sort_reference_file('%s.gtf' % args.outbase,args.annotation_format,printer)  # FIXME: automatic sorting is truncating some stuff?
    elif args.annotation_format == 'GFF3':
        outf = argsopener('%s.gff' % args.outbase, args)
        for tx in accepted_transcripts:
            outf.write(accepted_transcripts[tx].as_gff3())
        # printer.write('Sorting output...')
        # sort_reference_file('%s.gff' % args.outbase,args.annotation_format,printer)

    else:
        printer.write('Unable to find annotation format. Aborting...')

    # Additionally write an extended bed version of annotation, which is far more memory efficient
    printer.write("Writing extended bed...")
    outf = argsopener('%s.bed' % args.outbase, args)
    for tx in accepted_transcripts:
        outf.write(accepted_transcripts[tx].as_bed(extra_columns=['gene_id']))
    outf.close()
    # sort_reference_file('%s.bed' % args.outbase, 'BED', printer)

    printer.write("Done!")

if __name__ == '__main__':
    # main()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FileFormatWarning)
        main()