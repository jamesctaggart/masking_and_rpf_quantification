from plastid import SegmentChain
import numpy as np


def bootstrap_error_bars(profile, number_drawn = 1000, masked=False):
    '''
    Generate error bars on read count density by resampling read count distributions for each gene. 
    Return the 25th/75th percentiles and the full distribution of bootstrapped means
    '''
    drawn_means = []
    if len(profile) == 0:
        return None, None, None
    for ii in range(number_drawn):
        if masked:
            sample = np.random.choice(profile[~profile.mask], len(profile[~profile.mask]))
        else:
            sample = np.random.choice(profile, len(profile))
        drawn_means.append(np.mean(sample))

    return np.percentile(drawn_means, 25), np.percentile(drawn_means, 75), drawn_means


def get_common_exons(transcript_list):  
    # Input list of plastid Transcript objects
    if len(transcript_list)==1:
        return transcript_list[0].get_cds()
    out = SegmentChain()
    for seg in transcript_list[0].get_cds():
        add = True
        for tx in transcript_list:
            if seg not in tx.get_cds():
                add = False
        if add:
            out.add_segments(seg)
    return out  # Returns SegmentChain containing just the exons common to transcripts


def has_different_coding_regions(transcript_list):
    """
    Given a list of plastid transcript objects, return True if
    there are at least two non-identical CDS regions among them.
    :param transcript_list:
    :return:
    """
    previous_cds = transcript_list[0].get_cds()
    for tx in transcript_list[1:]:
        if not identical_cds(previous_cds,tx.get_cds()):
            return True
    return False


def calculate_masked_fraction(mask_hash, gene_segment_chain):
    """

    :param mask_hash: GenomeHash containing mask
    :param gene_segment_chain: SegmentChain for gene's CDS
    :return:
    """
    mask_set = set()

    for segment in mask_hash[gene_segment_chain]:
        mask_set = mask_set.union(segment.get_position_set())

    masked_positions = mask_set & gene_segment_chain.get_position_set()

    return float(len(masked_positions)) / gene_segment_chain.get_length()
