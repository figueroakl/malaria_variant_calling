#!/usr/bin/env python

import re
import os
import sys
import json
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
def qcheck(metaPath, refPath, path_to_picard):
	if os.path.isfile(metaPath):
        metadf = pd.read_csv(metaPath, sep = '\t')
        if metadf.empty:
            raise Exception('Input Meta File Found empty')
        else:
        	for num in range(len(metadf)):
        		sampleid = metadf.iloc[num,0]
				ipath = metadf.iloc[num,1]
                gc_met = collectGCBiasMetrics(path_to_picard, sampleid, ipath, refPath)
                ins_met = collectInsertSizeMetrics(path_to_picard, sampleid, ipath, refPath)
				

    else:
    	raise Exception('Input Metafile not found: The path mentioned is'.format(metaPath))
	return()

# to compute GC Bias Metrics
def collectGCBiasMetrics(path_to_picard, sampleid, ipath, refPath):
    opath = sampleid + "_gc_bias_metrics.txt"
    cmd = path_to_picard + "CollectGcBiasMetrics I=" + ipath + " O=" + opath + " R=" + refPath
    err = os.system(cmd)
    if err:
        print 'Execution of "%s" failed!\n' % cmd
        sys.exit(1)
    return(opath)

# to compute Insert Size Metrics
def collectInsertSizeMetrics(path_to_picard, sampleid, ipath):
    opath = sampleid + "insert_size_metrics.txt"
    cmd = path_to_picard + "CollectInsertSizeMetrics I=" + ipath + " O=" + opath
    err = os.system(cmd)
    if err:
        print 'Execution of "%s" failed!\n' % cmd
        sys.exit(1)
    return(opath)

