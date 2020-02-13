#!/usr/bin/env python
import argparse
import os

from variant_calling_modules_repo \
import collectGCBiasMetrics, \
collectInsertSizeMetrics, \
qcheck


def main():
    # init cmd-line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_path', required=True, help="Path to sample bam file to run QC metrics")
    parser.add_argument('--sampleid', required=True, help="sampleid assigned to the bamfile")
    parser.add_argument('--reference', required=True, help="Path to the reference genome fasta file")
    parser.add_argument('--path_to_picard', required=True, help="Path to picard cmd line tools")
    parser.add_argument('--path_to_dir', required=True, help="Path to intermediate files")
    parser.add_argument('--output_dir', required=True, help="Path to the output directory")
    args = parser.parse_args()

    bamfile = args.sample_path
    sampleid = args.sampleid
    ref = args.reference
    picard = args.path_to_picard
    tmp_dir = args.path_to_dir
    outdir = args.output_dir

    qcheck(bamfile, sampleid, tmp_dir, ref, picard, outdir)

    return()
    
if __name__ == "__main__":
    main()
