# Masking and quantification from Taggart et al. 2018, Cell Systems

Scripts used for quantification of protein synthesis from a ribosome profiling dataset.
Included Python scripts allow for 1. the generation of a mask file against a transcriptome to avoid counting multimapping reads and 2. quantification of protein synthesis from each gene using this mask

## Prefiltering annotations and mask generation

The included scripts will take a fasta-formatted genome and GTF2-formatted annotation from [Ensembl](https://useast.ensembl.org/info/data/ftp/index.html), and an optional list of transcript isoforms from [APPRIS](http://apprisws.bioinfo.cnio.es/pub/), filter out overlapping genes and those not included in the transcript isoform annotation, and generate a BED-formatted transcriptome mask reference file from them. Scripts to generate a mask file are based off of the [plastid package](https://plastid.readthedocs.io/en/latest/index.html) (version 0.4.8), published by Joshua Dunn and Jonathan Weissman. 

We mask against the transcriptome rather than the genome to both accelerate the runtime of mask generation and ensure correct masking near splice junctions. The simpler but slower workflow which masks the genome instead of transcriptome is available through the [crossmap function](https://plastid.readthedocs.io/en/latest/examples/using_masks.html) of plastid.

The workflow used here to generate a mask is as follows:
First, filter the GTF2-formatted annotation file using `prefilter_annotation_file.py`. This script can be run without the `--remove_overlaps` flag to preserve genes which share overlapping coding sequence. An example call might look like: 
```
python prefilter_annotation_file.py\
 --appris $APPRIS_ANNOTATION\
 --annotation_files $ANNOTATION_FILE\
  --annotation_format GTF2 $OUTFILE_BASE
```

The output of this script must be sorted. This can be done as follows:
```
cat annotation_file.gtf | grep -v "#" | sort -k1,1 -k4,4n >annotation_file_sorted.gtf &
cat annotation_file.bed | grep -v "#" | sort -k1,1 -k2,2n >annotation_file_sorted.bed &
```

Using the sorted, filtered GTF2, the reference genome can now be filtered using `prefilter_genome.py`, removing any contigs which do not encode a transcript. A typical call to do this might look like:
```
python prefilter_genome.py --annotation_file $FILTERED_ANNOTATION\
 --annotation_format 'GTF2'\
  $REFERENCE_GENOME $OUTFILE_BASE
```

From this filtered genome, build a bowtie1 index for use in mask generation. You can then run `crossmap_tophat.py` to generate a mask file from your filtered genome and annotation. To reduce run time, annotation file can be provided as both an extended BED file (used in making k-mers) and a GTF2 file (required for tophat mapping of k-mers). This mask will be specific to the transcriptome assembly you provide. The number of allowed mismatches is specified with `--mismatches`, and `--offset` must match the offset applied to mapped read positions used in downstream quantification (13 nt by default for included scripts) An example call:
```
python crossmap_tophat.py -k 27 --mismatches 1 --offset 13 --processes 2 --sequence_format fasta\
 --sequence_file $FILTERED_GENOME_FASTA\
  --bed_file $ANNOTATION_BED_FILE\
   --annotation_files $ANNOTATION_GTF_FILE\
    --annotation_format 'GTF2' $BOWTIE1_INDEX\
     $OUTFILE_BASE
```

The way k-mers are generated in this process means the ends of a transcript (size defined by k-mer length) will not be included in the mask, but this does not substantively impact protein synthesis quantification due to the end-trimming used in this process.


## Protein synthesis quantification from ribosome profiling
Included scripts allow for quantification of protein synthesis given an indexed BAM file of aligned ribosome protected fragments, genome and annotation, and mask file. Given these inputs, the scripts will produce a table of synthesis rates, as well as bootstrapped estimates of error on these synthesis rates. Error is reported as 25th and 75th percentiles of a bootstrapped distribution of synthesis rates given each gene's read count distribution. The output table additionally contains columns indicating if a particular gene has multiple possible coding regions and if a gene has fewer than 128 reads mapping to its quantified positions.

To quantify synthesis, run `quantify_bam.py`. By default, this script considers the 5' end of each read, with a +13 nt offset, and considers only reads 27-33 nt in length. The number of nucleotides to remove from the ends of each coding region can be specified with the `--trim5` and `--trim3` flags.

This script will import the mask and data from the BAM file (with mask file specified by the `--mask_file` flag), extract the common exons for all isoforms from each gene, calculate a metagene correction from the data, and subsequently calculate masked, metagene-corrected read density for each gene.

An a quantification call might look as follows:
```
python quantification/quantify_bam.py --trim5 15 --trim3 15 \
--annotation_files $ANNOTATION_BED\
 --annotation_format "BED" --bed_extra_columns "gene_id"\
 --mask_file $MASK_FILE\
  $RPF_BAM_FILE $OUTPUT_BASE
```

