import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mstats
from utils import has_different_coding_regions, calculate_masked_fraction


def smooth_profile(profile,averaging_window):
    """
    Given an averaging window radius, average the profile at each position over
    the values within this radius.
    :param profile: Profile (reads per position) to smooth
    :param averaging_window: Averages each position by the mean of range defined by +/- this value
    :return: Smoothed profile
    """
    smoothed = []

    for ii in range(len(profile)): # Iterate over positions of profile

        if ii<averaging_window:  # Enter if position considered is less than window radius

            val = np.mean([profile[0]]*(averaging_window-ii)+profile[0:ii+averaging_window])
            smoothed.append(val)

        elif (len(profile)-ii)<averaging_window:  # Enter if curve is within averaging window of profile end
            val = np.mean(profile[ii-averaging_window:]+[profile[-1]]*(ii+averaging_window-len(profile)))
            smoothed.append(val)

        else:  # Enter here for all positions in middle of profile
            smoothed.append(np.median(profile[ii-averaging_window:ii+averaging_window]))

    return smoothed


def calculate_metagene_correction(ga, gene_dictionary, nt_trimmed_5, nt_trimmed_3, L, outdir, mask = None):    # rawdata formatted as dictionary of profiles
    """

    :param ga: genome alignments (imported from bam file)
    :param gene_dictionary: Dictionary -- {gene name : [isoforms]}
    :param nt_trimmed: nucleotides trimmed from 5' and 3' end of profile
    :param L: Length of calculated metagene region
    :param outdir: Directory to save metagene plot to.
    :param mask: GenomeHash containing mask.
    :return:
    """
    reads_each_position = [[] for x in range(L)]  # generate a list of positions, which will contain a list storing the value at that position for each gene

    # normalize each profile to the average read count per position
    counter = 0
    for gene in gene_dictionary:
        counter +=1
        if counter%100 == 0:
            print "Genes processed: %s" % str(counter)

        primary_isoform = None
        isoforms = gene_dictionary[gene]
        if len(isoforms)==0:
            continue
        if len(isoforms)==1:
            primary_isoform = isoforms[0]
        if len(isoforms)>1:
            if not has_different_coding_regions(isoforms):  # If all isoforms have identical CDS, can use any, so choose first
                primary_isoform = isoforms[0]
            else:
                continue  # Skip anything with multiple distinct CDS

        profile = ga[primary_isoform.get_cds()]  # np array of values in CDS of primary isoform
        if sum(profile) < 128:  # Depth cutoff
            continue

        if mask != None:  # If we have a mask, exclude gene w/ primary isoform >50% masked.
            if calculate_masked_fraction(mask,primary_isoform)>0.5: 
                continue

        profile = profile[nt_trimmed_5:-nt_trimmed_3] # Calculate using trimmed profiles

        grouped_profile = [profile[ii:ii + 3] for ii in range(0, len(profile), 3)]  # Smooth over codons
        codon_signals = []
        for codon in grouped_profile:
            if len(codon) != 3:
                print 'Warning: length of gene %s not multiple of 3' % gene
            codon_signals.append(np.mean(codon))
        smoothed_profile = []
        for val in codon_signals:
            for i in range(3):
                smoothed_profile.append(val)

        average_reads = np.mean(smoothed_profile[:150])  # Normalize to first 50 codons
        # average_reads = np.mean(smoothed_profile)  # Use this if you want to average each over whole gene
        if average_reads == 0:
            continue
        normalized_profile = smoothed_profile / average_reads
        # add the reads at each position to our metagene storage list
        for ii in range(len(reads_each_position)):
            if ii < len(normalized_profile):
                reads_each_position[ii].append(normalized_profile[ii])

    # create and populate metagene plot
    metagene = []
    print "Calculating average per position..."
    for ii in range(len(reads_each_position)):
        for_mean = mstats.winsorize(reads_each_position[ii], limits=(0.05,0.05))
        metagene.append(np.mean(for_mean))  
        
    print 'Number of genes at positon 0: '+str(len(reads_each_position[0]))

    print "Smoothing metagene..."
    smoothed_metagene = smooth_profile(metagene,25)  # Smooth over 51 nt window (pos plus 25 each side)

    print "Plotting metagene..."
    plt.plot(range(len(metagene)),metagene)
    plt.xlabel('Position relative to start (nt)')
    plt.ylabel('Mean density')
    plt.plot(range(L),smoothed_metagene)
    plt.savefig(outdir+'/metagene.pdf',format='pdf',dpi=250)

    return smoothed_metagene


def apply_metagene_correction(values,metagene_profile):
    # Apply metagene correction to read counts at each position in a gene's cds

    # Values in form of np.array
    # gene_cds in form of plastid SegmentChain (transcript.get_cds())
    # metagene profile is np.array or list of correction values

    cds_length = len(values)
    if cds_length > len(metagene_profile):
        end_value = np.median(metagene_profile[-100:])
        metagene_correction = metagene_profile+[end_value]*(cds_length-len(metagene_profile))
    else:
        metagene_correction = metagene_profile[:cds_length]

    return values / metagene_correction
