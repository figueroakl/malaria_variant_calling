#!/usr/bin/env python
import argparse

from variant_calling_modules_repo \
    import collectGCBiasMetrics, \
    collectInsertSizeMetrics, \
    qcheck


def main():
    # init cmd-line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--metafile_path', required=True, help="Path to metafile excel or any other tsv")
    parser.add_argument('--reference', required=True, help="Path to the reference genome fasta file")
    parser.add_argument('--path_to_picard', required=True, help="Path to picard cmd line tools")
    #parser.add_argument('--upload_to_GC', default='gs://bucket/', type=str,
    #                    help="Directory to upload ouput on GC.")
    parser.add_argument('--verbose', action="store_true") # verbose replace verbose in if with dry run
    parser.add_argument('--dry_run', action="store_true")
    args = parser.parse_args()

    metafile = args.metafile_path
    ref = args.reference
    picard = args.path_to_picard
    
    qcheck(metafile, ref, picard)
    

