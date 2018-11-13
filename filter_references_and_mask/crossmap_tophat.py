#!/usr/bin/env python
"""
jct edits:
This script has been modified by JCT to use TopHat, rather than bowtie, to do the mapping.
It additionally will create k-mers from the transcriptome, guided by a GFF file, rather than genome.
Compare to plastid script crossmap.py, written by Joshua Dunn.

original jdunn docstring follows:
This script empirically determines which positions in a genome yield
multimapping reads under a given set of alignment parameters. These positions
are saved in a `BED`_-formatted :term:`mask file`, so that they may be excluded
as from further analyses.

To identify multimapping regions, a genome sequence is diced into :term:`k-mers <k-mer>`
and the :term:`k-mers` aligned back to the genome. Positions in the genome that
give rise to :term:`k-mers <k-mer>` that align equally well to more than one
genomic location are then marked as multimapping..

`k` is specified by the user, as are the alignment parameters.


Output files
------------
The following files are made:

    OUTBASE_READLENGTH_MISMATCHES_crossmap.bed
        Final :term:`mask file` annotation, in `BED`_ format
    
    OUTBASE_READLENGTH_MISMATCHES_CHROMOSOME_kmers.fa
        :term:`K-mers <k-mer>` derived from chromosome `CHROMOSOME`. These files can
        be reused in subsequent runs allowing a different number of mismatches,
        using the ``--have_kmers`` option

where:

  - `OUTBASE` is a name meaningful to the user
  - `READLENGTH` is the :term:`k-mer` length chosen by the user
  - `MISMATCHES` is the number of mismatches permitted during alignment,
    also set by the user.

   

Considerations for large genomes
--------------------------------

For large genomes (e.g. vertebrate, plant, or some *very* big amoebas):

  - |crossmap| can require a ton of memory if genome sequence is stored 
    in a fasta file. If |crossmap| maxes out your system's memory, it may
    be terminated by your system before it completes.
   
    Consider converting the file to a `2bit`_ file to save memory and
    avoid this potential problem

  - Enabling mismatches with short read sizes will make |crossmap| take
    a lot longer, especially on large genomes. Consider using ``--mismatches
    0`` if you run into this problem 

  - Using multiple processes (e.g. via ``-p 2``) will speed |crossmap|'s
    execution time, but will increase its memory footprint, as each process
    will need its own memory space to create and align k-mers from chromosomal
    sequence

  - By default, |crossmap| creates `BED`_ files. Consider converting these to
    `BigBed`_ files will save substantial amounts of time and memory in the future.

    This can be achieved using Jim Kent's ``bedToBigBed`` utility as follows
    (from the terminal)::

        $ bowtie-inspect --summary BOWTIE_INDEX | grep Sequence |\\
                         cut -f2,3 | sed -e "s/\([^ ]\+\).*\\t/\\1\\t/"  >OUTFILE.sizes
        $ sort -k1,1 -k2,2n OUTBASE.bed > OUTBASE_sorted.bed
        $ bedToBigBed OUTBASE_sorted.bed OUTBASE.sizes OUTBASE_sorted.bb


    See https://github.com/ENCODE-DCC/kentUtils/tree/master/src/product/scripts
    for download & documentation of Kent utilities
"""
__author__ = "joshua, modified by james"
import argparse
import sys
import os
import subprocess
import re
import inspect
import multiprocessing
import shutil
import functools
import pysam
from collections import OrderedDict
import kmers_from_genome_reference as makekmer

from plastid.util.io.filters import NameDateWriter, AbstractReader
from plastid.util.io.openers import get_short_name, argsopener
from plastid.genomics.roitools import SegmentChain, positionlist_to_segments, GenomicSegment
from plastid.util.scriptlib.help_formatters import format_module_docstring
from plastid.util.services.mini2to3 import xrange
from plastid.util.services.exceptions import MalformedFileError
from plastid.util.scriptlib.argparsers import SequenceParser, BaseParser
from plastid.genomics.genome_hash import GenomeHash
from Bio import SeqIO

namepat = re.compile(r"(.*):([0-9]+)\(([\+-])\)")  # Anything:(Any#)(+) -- as in original crossmap, can only handle read names on the + strand!
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


BigBedMessage = """Crossmap complete and saved as 'OUTFILE.bed'.

    For large (e.g. mammalian) genomes, it is highly recommended to convert
    the BED-format output to a BigBed file, using Jim Kent's bedToBigBed
    utility as follows (from the terminal):
    
        $ bowtie-inspect --summary BOWTIE_INDEX | grep Sequence |\\
                         cut -f2,3 | sed -e "s/\\([^ ]\+\\).*\\t/\\1\\t/"  >OUTFILE.sizes
        $ sort -k1,1 -k2,2n OUTFILE.bed > OUTFILE_sorted.bed
        $ bedToBigBed OUTFILE_sorted.bed OUTFILE.sizes OUTFILE_sorted.bb
    
    
    See https://github.com/ENCODE-DCC/kentUtils/tree/master/src/product/scripts
    for download & documentation of Kent utilities
    
    
"""


def simulate_reads(seq_record,fh=sys.stdout,k=30):
    """Chops a DNA sequence into :term:`k-mers <k-mer>`, mimicking a sequencing run.
    Output is delivered in fasta format. Sequences are named for position of
    origin using 0-based indices.
    
    Parameters
    ----------
    seq_record : :py:class:`Bio.SeqRecord.SeqRecord`
        DNA sequence
    
    fh : file-like
        filehandle to write output 
        
    k : int, optional
        length of k-mers to generate (Default: `30`)
    """
    seq = str(seq_record.seq)

    for x in xrange(0,len(seq_record)-k+1):
        fh.write(">%s:%s(+)\n" % (seq_record.name,x))
        fh.write("%s\n" % seq[x:x+k])

    return None


class FastaNameReader(AbstractReader):
    """Returns names of sequences in a fasta file"""

    def filter(self,line):
        """Return next sequence name in a fasta file

        Parameters
        ----------
        line : str
            Line of text

        Returns
        -------
        str
            Name of next sequence, excluding prefix `'>'` and line terminator
        """
        if line.startswith(">"):
            return line[1:].rstrip()
        else:
            return self.__next__()


def revcomp_mask_chain(seg,k,offset=0):
    """Reverse-complement a single-interval mask, correcting for `offset`.
    
    Parameters
    ----------
    seg : |SegmentChain|
        Plus-strand mask, including `offset`

    k : int
        Length of k-mers

    offset : int, optional
        Offset from 5' end of read at which to map mask (Default: `0`)

    Returns
    -------
    |SegmentChain|
        Mask on minus strand corresponding to `seg`
    """
# Algorithm note:
#
#     Let
#         FW = plus-strand coordinate
#         RC = minus-strand coordinate
#     
#     Then
#         RC = FW + k - 1 - offset
#     
#     But we are given FW + offset, so:
#     
#         RC + offset = (FW + offset) + k - 1 - offset
#         RC = (FW + offset) + k - 1 - 2*offset   
    span = seg.spanning_segment
    new_offset = k - 1 - 2*offset
    ivminus = GenomicSegment(span.chrom,
                             span.start + new_offset,
                             span.end + new_offset,
                             "-")
    return SegmentChain(ivminus)


def read_cigar_offset(cigar_tuples):
    """
    Since tophat returns the leftmost position on top strand, need to switch to 5' end of read for mask
    Use cigar string to incorporate read length and skipped genome bases if read spans intron
    Using pysam AlignedSegment.cigartuples as input
    Returns length to add to AlignedSegment.pos in mapping
    :param cigar_tuples:
    :return:
    """

    length = 0
    for tup in cigar_tuples:
        length+=tup[1]

    return length


# def set_to_segchain(chrom,position_set,strand):
#     """
#
#     :param chrom: Name of chromosome (str)
#     :param position_set: Set of positions on chromsome (each as int)
#     :param strand: Strand of chromosome of interest (str)
#     :return: SegmentChain object for that strand of that chromosome
#     """
#     sc = SegmentChain()
#     for run in runs(position_set):
#         sc.add_segments(GenomicSegment(chrom,run[0],run[1]+1,strand))
#     return sc

class ReferencePositionContainer(object):  # Object which contains genome positions which should be masked
    def __init__(self):
        self.chromosomes = dict()

    def add_position(self, position, strand, chromosome):
        # Convert strand to bool, accepting flexible inputs
        if strand == '+' or strand == 1 or strand == 'plus':  
            reverse_strand_bool = 0 
        elif strand == '-' or strand == 0 or strand == 'minus':
            reverse_strand_bool = 1

        if chromosome not in self.chromosomes:  # Populate chromsome when you add a position to it
            self.chromosomes[chromosome] = (set(), set())
        if reverse_strand_bool == 1:
            self.chromosomes[chromosome][1].add(position)
        elif reverse_strand_bool == 0:
            self.chromosomes[chromosome][0].add(position)

    def has_single_position(self):
        total_positions = 0
        for chromosome in self.chromosomes:
            total_positions += len(self.chromosomes[chromosome][0]) + len(self.chromosomes[chromosome][1])

        if total_positions > 1:
            return False
        elif total_positions == 1:
            return True
        else:
            raise RuntimeError("Position sets for read are empty.")

    # def generate_segment_chains(self):
    #     """
    #     Iterate through the chromosomes of the reference, and generate SegmentChains
    #     that correspond to the positions contained by this object.
    #
    #     :yield: SegmentChains for all chromosomes in the container object.
    #     """
    #     for chromosome in self.chromosomes:
    #         for strand_bool in range(2): # 0 is plus strand, 1 is minus strand
    #             positions = self.chromosomes[chromosome][strand_bool]
    #             if strand_bool == 0:  
    #                 strand = '+'
    #             elif strand_bool == 1:
    #                 strand = '-'
    #             yield set_to_segchain(chromosome,positions,strand)


def get_strand(read_name):
    return namepat.search(read_name).groups()[2]


def get_pos(read_name):
    return int(namepat.search(read_name).groups()[1])


def sort_read_names(read_name_list):
    plus_reads = []
    minus_reads = []

    for read in read_name_list:  # Split into plus and minus lists.
        if get_strand(read) == '+':
            plus_reads.append(read)
        elif get_strand(read) == '-':
            minus_reads.append(read)

    plus_reads.sort(key=lambda l: get_pos(l))
    minus_reads.sort(key=lambda l: get_pos(l))

    return plus_reads, minus_reads


def bam_to_bed(bam_file_name,k,offset=0):

    """
    Generates SegmentChain of all multimapping reads (toomany file), given a .bam alignment file
    This can be used in modified chrom_worker function, which saves them as bed files
    This is a replacement for fa_to_bed, as we cannot use the --max output from bowtie as our input (considering spliced alignment)

    Parameters
    ----------
    bam_file_name : bam file to be parsed for multiply mapping reads

    Yields
    ------
    Returns masked positions, but produces fasta file to be converted to bed by fa_to_bed
    :return: None
    """

    # multiply_mapped: keys are read names
    # values are two sets, containing positions the reads map to on plus and minus strands respectively
    # the idea is to establish all reads which multiply map
    multiply_mapped = OrderedDict()

    reader = pysam.AlignmentFile(bam_file_name)

    # Go through the alignments in the bam file and collect positions they map to.
    # This generates a new ReferencePositionContainer to be iterated through for each read!
    for alignment in reader:
        if alignment.mapq != 50:  # mapq of 50 means uniquely mapping, skip these.
            read_name = alignment.query_name
            cigar = alignment.cigartuples
            chromosome = alignment.reference_name
            multiply_mapped[read_name] = multiply_mapped.get(read_name,
                                                             ReferencePositionContainer())  # Add read to dict if not there

            if alignment.is_reverse:  
                p_to_add = alignment.pos + read_cigar_offset(cigar)
                multiply_mapped[read_name].add_position(p_to_add, '-', chromosome)
            else:
                p_to_add = alignment.pos
                multiply_mapped[read_name].add_position(p_to_add, '+', chromosome)


    # Iterate through multiply_mapped reads, generating SegmentChains for those which map to multiple locations.
    # These are fed to rest of chrom_worker, which saves them in bed file
    last_strand = None
    plus_reads, minus_reads = sort_read_names(multiply_mapped.keys())  

    # multif = open('multiplymap_testing','a')  # For checking positions

    for reads in (plus_reads,minus_reads):  # two sorted lists of read names to iterate through

        last_chrom = None  # last_chrom is previous read's chromosome
        last_pos = None  # last position considered
        start_pos = None  # start of this segment in bed file

        for read in reads:

            if multiply_mapped[read].has_single_position():  # Skip reads which multiply mapped to same transcript
                continue

            else:  # Will only enter here if read maps to multiple positions
                # multif.write(read+'\t')  # For validating positions are correct.
                # multif.write(str(multiply_mapped[read].chromosomes)+'\n')

                chrom, pos, strand = namepat.search(read).groups()  # Extract chromosome and position of read's origin from its name
                if strand == '+':
                    pos = int(pos) + offset  # Apply an offset if required (p-site offset, etc).
                elif strand == '-':
                    pos = int(pos) - offset
                else:
                    raise RuntimeError('Strand must be + or -')
                if chrom != last_chrom or strand != last_strand:  # Triggers this when it hits new chromosome, writing the last segment from the old one. we also want to trigger if strand changes.
                    if last_chrom is not None and last_strand is not None:
                        chain = SegmentChain(GenomicSegment(last_chrom, start_pos, last_pos + 1, strand))
                        last_chrom = chrom
                        last_strand = strand
                        start_pos = pos
                        last_pos = pos
                        yield chain  # yield region
                    else:
                        last_chrom = chrom
                        last_strand = strand
                        start_pos = pos
                        last_pos = pos

                else:  # Goes here if still on same chromosome and strand
                    delta = pos - last_pos
                    if delta > 1:  # If it skips more than 1 position, yield (and ultimately write) the segment
                        chain = SegmentChain(GenomicSegment(chrom, start_pos, last_pos + 1, strand))
                        last_pos = pos
                        start_pos = pos
                        yield chain
                    elif delta == 1:  # If the positions are contiguous, however, continue onto next read. (Preventing adjacent, len(1) segments in bed)
                        last_pos = pos
                    else:  # Triggers this if order of kmers is not correct (i.e. if positions are getting smaller)
                        msg = "k-mers are not sorted at read %s! Aborting." % read_name
                        raise MalformedFileError(bam_file_name, msg, line_num=read)

        # export final feature
        if last_pos is not None:
            chain = SegmentChain(GenomicSegment(chrom, start_pos, last_pos + 1, strand))

            yield chain
    # multif.close()


def chrom_worker(chrom_seq,args=None):
    name, seq_or_kmers = chrom_seq  # Chrom_seq is a tuple (chromosome name, SeqRecord)

    printer.write("Processing chromosome %s..." % name)
    base         = "%s_%s_%s_%s" % (args.outbase, args.read_length, args.mismatches, name)
    kmer_file    = "%s_%s_%s_kmers.fa"     % (args.outbase,args.read_length,name)  # This is changed from original crossmap - kmer files have shorter name
    toomany_file = "%s/accepted_hits.bam"  % name  # Changed from original crossmap.
    bed_file     = "%s_crossmap.bed" % base

    if args.have_kmers == False:
        if args.bed_file == None:
            # printer.write('Using annotation file to generate kmers.')
            # args.annotation_file is used to generate kmers if separate bed file is not provided
            # We can't use BED as args.annotation_files, since we need GFF/GTF for tophat.
            makekmer.main(argv=['--annotation_files', args.annotation_files, '--annotation_format', args.annotation_format,
                            '--genome_fh', args.sequence_file, '--outbase', args.outbase, '-k', str(args.read_length), name]) 
        else:
            # printer.write('Using bed file to generate kmers.')
            # If bed file is provided, we use it to make kmers rather than GFF, as it will be much more efficient.
            makekmer.main(argv=['--annotation_files', args.bed_file, '--annotation_format', 'BED',
                            '--genome_fh', args.sequence_file, '--outbase', args.outbase, '-k', str(args.read_length), name]) 

    else:
        kmer_file = seq_or_kmers

    argdict = { "mismatches" : args.mismatches,  # These are the args for use in tophat cmd
                "processors" : 1, 
                "tophat"     : args.tophat,
                "tophatoutdir"  : name,
                "kmers"      : kmer_file,
                "ebwt"       : args.ebwt,
                "annotation" : args.annotation_files
                }

    cmd  = "%(tophat)s -G %(annotation)s --bowtie1 --o %(tophatoutdir)s " \
           "--transcriptome-only --no-novel-juncs -N %(mismatches)s " \
           "-p %(processors)s %(ebwt)s %(kmers)s " % argdict

    printer.write("Aligning %s-mers for chromosome '%s' :\n\t'%s'" % (args.read_length,name,cmd))
    try:
        retcode = subprocess.call(cmd,shell=True)
        if retcode < 0 or retcode == 2:
            printer.write("Alignment for chromosome '%s' terminated with status %s" % (name,retcode))
        else:
            if os.path.exists(toomany_file):
                printer.write("Assembling multimappers from chromosome '%s' into crossmap..."% name)
                with argsopener(bed_file,args,"w") as bed_out:
                    for chain in bam_to_bed(open(toomany_file),  # iterate through generator, producing all bed segments
                                                         args.read_length,
                                                         offset=args.offset):
                        bed_out.write(chain.as_bed())

                    # for plus_chain, minus_chain in bam_to_bed(open(toomany_file),  # iterate through generator, producing all bed segments
                    #                                      args.read_length,  # THIS IS THE OLD UNSTRANDED VERSION
                    #                                      offset=args.offset):
                    #     bed_out.write(plus_chain.as_bed())
                    #     bed_out.write(minus_chain.as_bed())
                
                    bed_out.close()
            
            else:
                printer.write("Could not find multimapper source file '%s' ." % toomany_file)
    except OSError as e:
        printer.write("Alignment failed for chromosome '%s': %s" % (name,e))

    printer.write("Cleaning up chromosome '%s'..." % name)
    os.remove(toomany_file)  # Just delete accepted_hits.bam
    if args.have_kmers == False and args.save_kmers == False:
        os.remove(kmer_file)

    return bed_file


def main(argv=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
	argv : list, optional
		A list of command-line arguments, which will be processed
		as if the script were called from the command line if
		:py:func:`main` is called directly.

        Default: `sys.argv[1:]`. The command-line arguments, if the script is
        invoked from the command line
    """
    sp = SequenceParser()
    bp = BaseParser()

    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[bp.get_parser(),sp.get_parser()])
    parser.add_argument("-k",dest="read_length",metavar="READ_LENGTH",
                        type=int,default=29,
                        help="K-mer length to generate from input file. "+
                             "(Default: 29)")
    parser.add_argument("--offset",type=int,default=0,
                        help="Offset from 5' end of plus-strand read at which to attribute score (Default: 14)")
    parser.add_argument("--mismatches",metavar="N",
                        type=int,default=0,
                        help="Number of mismatches tolerated in alignment. "+
                           "(Default: 0)")
    parser.add_argument("--tophat",dest="tophat",default="/usr/local/bin/tophat2/tophat2",
                        type=str,
                        help="Location of tophat binary (Default: ``/usr/local/bin/tophat2/tophat2``)")
    parser.add_argument("--have_kmers",default=False,action="store_true",
                        help="If specified, use k-mer files from previous run. "+\
                             " In this case 'sequence_file' should be the value "+\
                             "'outbase' from the k-mer files you want to use.")
    parser.add_argument("--save_kmers",default=False,action="store_true",
                        help="Save k-mer files for reuse in a subsequent run.")
    parser.add_argument("-p","--processes",type=int,default=2,metavar="n",
                        help="Number of processes to use (should be <= number of chromosomes")
    parser.add_argument("--annotation_files",type=str,
                        help="Annotation file for genome to make k-mers")
    parser.add_argument("--annotation_format", type=str,
                        help="Annotation file format for genome to make k-mers")
    parser.add_argument("--bed_file", type = str,
                        help = "If specified, provide bed file for faster kmer generation",)
    parser.add_argument("ebwt",type=str,
                        help="Bowtie index of genome against which crossmap will be made. In most cases, should be generated from the same sequences that are in `sequence_file`.")
    parser.add_argument("outbase",default='TopXMap',type=str,
                        help="Basename for output files")

    # print parser
    args = parser.parse_args(argv)

    bp.get_base_ops_from_args(args)

    #filenames
    base         = "%s_%s_%s" % (args.outbase, args.read_length, args.mismatches)
    bed_file     = "%s_topxmap.bed" % base

    #if not os.path.exists(args.sequence_file):
    #    printer.write("Could not find source file: %s" % args.sequence_file)
    #    printer.write("Exiting.")
    #    sys.exit(1)
    printer.write('Getting kmers...')
    if args.have_kmers == True:
        import glob
        kmer_files = glob.glob(args.outbase+"*_kmers.fa")
        seq_pat = re.compile(r".*_([^_]*)_kmers.fa") #
        seqs = { seq_pat.search(X).groups()[0] : X for X in kmer_files }
    else:
        seqs = sp.get_seqdict_from_args(args,index=True) 

    worker = functools.partial(chrom_worker,args=args)
    chroms = seqs.items()  # Generator that produces the chromosomes to apply the chrom_worker to. Output is tuple: (val,key)
    pool = multiprocessing.Pool(processes=args.processes)
    bed_filenames = pool.map(worker,chroms,1)  # Apply worker to chromosomes - bedfiles generated are named based on chr of reads origin
    pool.close()
    pool.join()
   
    with open(bed_file,"w") as fout:  # Merge bed files into one.
        for f in sorted(bed_filenames):
            shutil.copyfileobj(open(f,"r"),fout)
            os.remove(f)

    fout.close()

    printer.write("Done.")
    printer.write(BigBedMessage.replace("OUTFILE",bed_file.replace(".bed","")).replace("BOWTIE_INDEX",args.ebwt))


if __name__ == "__main__":
    main()
