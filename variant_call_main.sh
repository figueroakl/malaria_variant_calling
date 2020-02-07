#!/bin/bash

use Java-1.8
use Picard-Tools

meta_path="/path/to/metafile/"
ref_path="/path/to/ref/"
path_to_fq="/path/to/fq"

fastqarray=(`find $path_to_fq -maxdepth 1 -name "*.fastq"`)
fqarray=(`find $path_to_fq -maxdepth 1 -name "*.fq"`)

if [ ${#fastqarray[@]} -gt 0 ]; then 
    
elif [ ${#fqarray[@]} -gt 0  ]; then 
    
else
	echo false
fi