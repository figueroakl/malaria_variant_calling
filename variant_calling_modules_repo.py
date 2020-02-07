#!/usr/bin/env python

import re
import os
import sys
import json
import warnings
import argparse
import numpy as np
import pandas as pd
import datetime as dt
import openpyxl as excel
#from google.cloud import storage

#example module: to align
def bwaAlign():
    return()

# QC metrics wrapper
def qcheck(metaPath, refPath, path_to_picard, output_dir):
	if os.path.isfile(metaPath):
        metadf = pd.read_csv(metaPath, sep = '\t')
        if metadf.empty:
            raise Exception('Input Meta File Found empty')
        else:
            fqc = open("qc_summary.tsv", "a+")
        	for num in range(len(metadf)):
        		sampleid = metadf.iloc[num,0]
				ipath = metadf.iloc[num,1]
                qc_summary = pd.DataFrame(data = {'sampleid' : [sampleid], 'MaxCoverage_GC' : [np.nan], 'MaxCoverage_INS' : [np.nan]})
                gc_met = collectGCBiasMetrics(path_to_picard, sampleid, ipath, refPath, output_dir)
                ins_met = collectInsertSizeMetrics(path_to_picard, sampleid, ipath, refPath, output_dir)
				if os.path.isfile(gc_met):
                    gc_df = pd.read_csv(gc_met, sep = '\t', skiprows = 6)
                    maxgc = gc_df['NORMALIZED_COVERAGE'].idxmax()
                    qc_summary['NORMALIZED_COVERAGE'] = maxgc
                else:
                    warnings.warn('Path to GcMetrics output File does not exist: Path given'.format(gc_met))
                if os.path.isfile(ins_met):
                    ins_df = pd.read_csv(ins_met, sep = '\t')
                else:
                    warnings.warn('Path to InsertsizeMetrics File does not exist: Path given'.format(ins_met))
                fqc.write(qc_summary)
            fqc.close()
    else:
    	raise Exception('Input Metafile not found: The path mentioned is'.format(metaPath))
    return()

# to compute GC Bias Metrics
def collectGCBiasMetrics(path_to_picard, sampleid, ipath, refPath, output_dir):
    opath = output_dir + sampleid + "_gc_bias_metrics.txt"
    chart = " CHART_OUTPUT=" + output_dir + sampleid + "_chart_output.pdf"
    summary = " s=" output_dir + sampleid + "_summary_output.txt"
    cmd = path_to_picard + "CollectGcBiasMetrics I=" + ipath + " O=" + opath + " R=" + refPath + chart + summary
    err = os.system(cmd)
    if err:
        print 'Execution of "%s" failed!\n' % cmd
        sys.exit(1)
    return(opath)

# to compute Insert Size Metrics
def collectInsertSizeMetrics(path_to_picard, sampleid, ipath, output_dir):
    opath = output_dir + sampleid + "insert_size_metrics.txt"
    cmd = path_to_picard + "CollectInsertSizeMetrics I=" + ipath + " O=" + opath
    err = os.system(cmd)
    if err:
        print 'Execution of "%s" failed!\n' % cmd
        sys.exit(1)
    return(opath)

def depthOfCoverage():
    return()

