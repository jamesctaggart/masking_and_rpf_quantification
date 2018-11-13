import argparse
import inspect
import gc
import os
import sys
from shutil import copyfileobj
from Bio import SeqIO
from plastid.util.io.openers import get_short_name, argsopener
from plastid.util.scriptlib.argparsers import AnnotationParser
from plastid.util.io.filters import NameDateWriter, AbstractReader

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

def cat_files(fh_list,output_fh,save_partial=False):
    with open(output_fh, 'wb') as output:
        for fh in fh_list:
            with open(fh,'rb') as fd:
                copyfileobj(fd, output, 1024*1024*10)
            if not save_partial:
                os.remove(fh)


def write_kmers(gene_dict,k,genome_dict,out_file,make_fastq=False):
    """This should generate all unique reads which could be produced from \
    each genomic position encoding a transcript.
    These can then be fed into TopHat, and positions which produced \
    Multimapping reads can be excluded from analysis (via mask .bed file).
    """
    outf = open(out_file,'w')
    gene_counter = 0
    chroms_not_in_genome = []
    for gene in gene_dict:

        gene_counter += 1  # Track how many genes
        if gene_counter%2500 == 0:
            print 'Number of genes processed: ', gene_counter

        isoforms = gene_dict[gene]

        reads_included = set()  # This will prevent double-writing reads of identical position/sequence

        for tx in isoforms:

            gene_name = tx.get_gene()
            chromosome = tx.chrom
            tx_name = tx.get_name()
            if chromosome in chroms_not_in_genome:
                continue
            try:
                seq = tx.get_sequence(genome_dict)
            except(KeyError):
                printer.write('Chromosome/contig %s not found in genome. Skipping... ' % chromosome)
                chroms_not_in_genome.append(chromosome)
                continue
            genome_positions = tx.get_position_list()

            strand = tx.strand

            if make_fastq == False:
                for ii in range(len(seq) - k):
                    read = seq[ii:ii + k]
                    if strand == '-':
                        pos5pr = genome_positions[-(ii+1)]
                    elif strand == '+':
                        pos5pr = genome_positions[ii]
                    if (pos5pr,read) not in reads_included:  # check if pos,seq of read already in file
                        outf.write('\n'+'>%s:%s(%s)' % (chromosome,pos5pr,strand))
                        outf.write('\n'+read)
                        reads_included.add((pos5pr,read))

            elif make_fastq == True:  # TODO: This may give errors - only updated make_fastq=False so far
                for ii in range(len(seq) - k):
                    read = seq[ii:ii + k]
                    if strand == '-':
                        pos5pr = genome_positions[-(ii+1)]
                    elif strand == '+':
                        pos5pr = genome_positions[ii]
                    if (pos5pr,read) not in reads_included:  # check if pos,seq of read already in file
                        out_file.write('\n'+'@'+gene_name+'_'+tx_name+'_'+str(pos5pr)+':'+chromosome+'('+strand+')')
                        out_file.write('\n'+read)
                        out_file.write('\n'+'+'+gene_name+'_'+tx_name+'_'+str(pos5pr)+':'+chromosome+'('+strand+')')
                        out_file.write('\n'+'I'*len(read))
                        reads_included.add((pos5pr,read))
    outf.close()


def main(argv=sys.argv[1:]):
    """This function will generate kmers from a reference file
Need kmers to be sorted by position within each chromosome for crossmap to successfully import them
"""
    # TODO: pipe arguments into this more gracefully than manually writing out argv in crossmap_tophat.py...
    chrom_genes = {}

    ap = AnnotationParser()
    annotation_file_parser = ap.get_parser()
    parser = argparse.ArgumentParser(parents=[annotation_file_parser])
    parser.add_argument('--genome_fh',metavar='genome', type=str,
                        help="Genome fasta file to retrieve sequences")
    # parser.add_argument("outbase", metavar="outbase", type=str,
    #              help="Basename for output file")
    parser.add_argument("-k",dest="read_length",metavar="READ_LENGTH",
                        type=int,default=29,
                        help="K-mer length to generate from input file. "+
                             "(Default: 29)")
    parser.add_argument('--outbase',dest="outbase",metavar='OUTBASE',
                        type=str,default='TopXMap')
    parser.add_argument('chromosome',type=str)
    args = parser.parse_args(argv)


    source = ap.get_transcripts_from_args(args, printer=printer)  # Iterator w/ transcript info
    # code here adapted from cs generate
    # loop conditions
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome_fh, 'fasta'))

    wrote_something = False  # Bool so we can create an empty kmer file if there are no segments in refs for chromosome
    last_chrom = None
    do_loop = True

    # to save memory, we process one chromosome at a time if input file is sorted
    # knowing that at that moment all transcript parts are assembled
    while do_loop == True:
        try:
            tx = next(source)
        except StopIteration:
            do_loop = False

        try:
            # if chromosome is completely processed or EOF
            if tx.spanning_segment.chrom != last_chrom or do_loop == False:
                if last_chrom is not None or do_loop == False:
                    if last_chrom == args.chromosome:
                        printer.write("Writing kmers for chromosome/contig '%s'..." % last_chrom)
                        write_kmers(chrom_genes,int(args.read_length),genome_dict,'%s_%s_%s_kmers.fa' % (args.outbase,args.read_length,last_chrom))
                        wrote_something = True
                        do_loop = False  # If it has written the kmers for the chromosome of interest, break while loop.
                    del chrom_genes
                    gc.collect()
                    del gc.garbage[:]
                    chrom_genes = {}

                # reset last chrom
                last_chrom = tx.spanning_segment.chrom

            # otherwise, remember transcript
            else:
                if tx.get_gene() not in chrom_genes:
                    chrom_genes[tx.get_gene()] = [tx]
                else:
                    chrom_genes[tx.get_gene()].append(tx)

        # exit gracefully if no transcripts found
        except UnboundLocalError:
            print "UNBOUND LOCAL ERROR"
            pass

    if not wrote_something:
        #  If there were no annotated transcripts on chromosome/contig (happens a lot with APPRIS filter)
        #  use this to just generate an empty file, so tophat won't throw an error.
        fh = '%s_%s_%s_kmers.fa' % (args.outbase, args.read_length, args.chromosome)
        print "No transcripts found on this chromosome! Generating empty kmer file..."
        print "Warning: expect an IndexError from Tophat call."
        open(fh, 'a').close()


if __name__ == '__main__':
    main()