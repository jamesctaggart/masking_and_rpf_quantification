from plastid.util.scriptlib.argparsers import AnnotationParser,\
                                              AlignmentParser,\
                                              MaskParser,\
                                              PlottingParser, \
                                              BaseParser
import os
import inspect
import pickle
import argparse
import numpy as np
from plastid.genomics.genome_hash import GenomeHash, BigBedGenomeHash
from plastid.util.io.filters import NameDateWriter
from plastid.util.io.openers import get_short_name, argsopener
from plastid.genomics.genome_array import BAMGenomeArray
from plastid.genomics.map_factories import FivePrimeMapFactory, SizeFilterFactory
from plastid.readers.bed import BED_Reader
from corrections import calculate_metagene_correction, apply_metagene_correction
from utils import get_common_exons, has_different_coding_regions, bootstrap_error_bars

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


def main():
    ap = AnnotationParser()
    annotation_file_parser = ap.get_parser()
    parser = argparse.ArgumentParser(parents=[annotation_file_parser])
    parser.add_argument('--trim5', metavar='trim5', type = int, default = 0,
                        help = "Number of nt to trim from 5' end of profiles")
    parser.add_argument('--trim3', metavar='trim3', type = int, default = 0,
                        help = "Number of nt to trim from 3' end of profiles")
    parser.add_argument("--mask_file", metavar='mask', type = str,
                        help = "Mask file, calculated with crossmap_tophat.py")
    parser.add_argument("--mask_bb", metavar = 'maskbb', type = bool,
                        default = False, help = "Specify true if mask formatted as BigBed file")
    parser.add_argument("input_bam", type = str,
                        help = "bam file containing ribosome profiling alignments.")
    parser.add_argument("outbase", metavar="outbase", type=str,
                        help = "Basename for output file")
    args = parser.parse_args()

    OUTDIR = 'synthesis_tables/'+args.outbase
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)

    # Read in .bam file, extract profiles as a SegmentChain
    # Sequencing alignments imported from tophat output (bam file).
    printer.write('Reading alignments...')
    alignments = BAMGenomeArray([args.input_bam])  # Read in alignment data.

    # Correct for psite offset
    alignments.set_mapping(FivePrimeMapFactory(offset=13))  # Set p-site offsets.

    # Filter read lengths
    alignments.add_filter('size',SizeFilterFactory(min=27,max=33))

    printer.write('Initializing genes...')

    source = ap.get_transcripts_from_args(args, printer=printer)  # Iterator w/ transcript info

    # Group all transcripts into lists by their parent gene
    gene_dict = {}  # Key = gene name, val = list of transcript objects
    for transcript in source:
        current_gene = transcript.get_gene()
        if current_gene in gene_dict:
            gene_dict[current_gene].append(transcript)
        else:
            gene_dict[current_gene] = [transcript]

    printer.write('Importing mask...')

    if not args.mask_bb:  # Go here if mask .bed formatted (i.e. NOT BigBed)
        mask_features = list(BED_Reader(args.mask_file))  # Extract masks from mask file
        mask_hash = GenomeHash(mask_features)  # Organize features into indexed GenomeHash
    else:  # If specified that mask is BigBed, go here.
        mask_hash = BigBedGenomeHash(args.mask_file)

    printer.write('Calculating metagene...')
    metagene_correction = calculate_metagene_correction(alignments, gene_dict, args.trim5,args.trim3, 3000,
                                                        OUTDIR, mask = mask_hash)

    synthesis_file = argsopener('%s/%s' % (OUTDIR,args.outbase),args,'w')
    synthesis_file.write('Gene Name' + '\t' + 'Synthesis Rate' + '\t' + 'Bootstrap25' + '\t' + 'Bootstrap75' + '\t' + 'Has Multiple CDS' + '\t' + 'Counts < 128')

    # troubleshooting_file = argsopener('%s/%s_extra' % (OUTDIR,args.outbase),args,'w')
    # troubleshooting_file.write('Gene Name' + '\t' + 'Synthesis Rate' '\t' + 'Bootstrap25' + '\t' + 'Bootstrap75' + '\t' + 'Unmasked_corrected' + '\t'
    #                            + 'Masked_uncorrected' + '\t' + 'Unmasked_uncorrected'+ '\t' + 'Has Multiple CDS'
                               # + '\t' + 'Counts < 128')

    printer.write('Writing synthesis rates...')

    counter = 0
    bootstrap_distributions = {}  # Store bootstrap distributions if you need to save them. Can remove to save memory.
    measured_distributions = {}  # Store measured distributions if you need to save them. Can remove to save memory.
    for gene in gene_dict:
        if counter % 500 == 0:
            print 'Writing gene number %s' % str(counter)
        isoforms = gene_dict[gene]
        shared_exons = get_common_exons(isoforms)

        # Identify masked positions for shared_exons
        masks = mask_hash[shared_exons]
        for mask_chain in masks:
            shared_exons.add_masks(*mask_chain)

        mask_count_array = shared_exons.get_masked_counts(alignments)
        mask_count_array = mask_count_array[args.trim5:-args.trim3]  # Trim profile - do this now, since metagene calculations use trimmed.

        uncorrected_profile = mask_count_array.data  # Unmasked data
        mask_bools = mask_count_array.mask  # Array of mask bools (True = exclude)
        masked_profile = np.zeros(len(mask_bools) - np.sum(mask_bools))  # initialize masked_profile
        uncorrected_masked_profile = np.zeros(len(mask_bools) - np.sum(mask_bools))  # To test effects of metagene correction, generating profile w/o
        
        # Apply metagene correction
        corrected_profile = apply_metagene_correction(uncorrected_profile, metagene_correction)

        # Build masked profile based on bool list in masked_count_array
        out_index = 0
        for ii in range(len(uncorrected_profile)):
            if mask_bools[ii]==False:
                masked_profile[out_index] = corrected_profile[ii]
                uncorrected_masked_profile[out_index] = uncorrected_profile[ii]
                out_index += 1
        if out_index!=len(masked_profile):
            raise RuntimeError('Error: mask profile not filled.')

        multiple_coding = has_different_coding_regions(isoforms)  # True if gene has multiple coding isoforms
        if sum(uncorrected_masked_profile) < 128:
            sub128 = True
        else:
            sub128 = False

        # Calculate synthesis rate, as well as those used for troubleshooting
        synthesis = np.mean(masked_profile)
        err25, err75, bootstrap_distributions[gene] = bootstrap_error_bars(masked_profile)
        measured_distributions[gene] = masked_profile
        unmasked_corrected_synthesis = np.mean(corrected_profile)
        masked_uncorrected_synthesis = np.mean(uncorrected_masked_profile)
        unmasked_uncorrected_synthesis = np.mean(uncorrected_profile)

        # Write to tables
        synthesis_file.write('\n'+str(gene)+'\t'+str(synthesis)+'\t'+str(err25)+'\t'+str(err75)+'\t'+str(multiple_coding)+'\t'+str(sub128))

        # troubleshooting_file.write('\n'+gene+'\t'+str(synthesis)+'\t'+str(err25)+'\t'+str(err75)+'\t'+str(unmasked_corrected_synthesis)+'\t'+
        #                            str(masked_uncorrected_synthesis)+'\t'+str(unmasked_uncorrected_synthesis)
        #                            + '\t' + str(multiple_coding) + '\t' + str(sub128))
        counter+=1

    synthesis_file.close()
    # troubleshooting_file.close()


if __name__ == '__main__':
    main()
