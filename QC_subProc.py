#!/usr/bin/env python
import argparse

from mal095_modules \
    import collectGCBiasMetrics, \
    collectInsertSizeMetrics


def main():
    # init cmd-line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--metafile_path', required=True, help="Path to metafile excel")
    parser.add_argument('--upload_to_GC', default='gs://bucket/', type=str,
                        help="Directory to upload ouput on GC.")
    parser.add_argument('--verbose', action="store_true") # verbose replace verbose in if with dry run
    parser.add_argument('--dry_run', action="store_true")
    args = parser.parse_args()



