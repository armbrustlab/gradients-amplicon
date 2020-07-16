#!/bin/bash

# run-qiime.sh START_ARTIFACT EARLY_STOP

# processing pipeline for 18S analysis in qiime2
# read each dataset into qiime2 artifact format using a
# manifest.

file=$1
filename=$(basename -- "$file")
filename_no_ext="${filename%.*}"

es=$2

mkdir int_files
mkdir run_stats
mkdir output_tables

# trim adapter sequences
qiime cutadapt trim-single --i-demultiplexed-sequences "${file}" \
--p-cores 8 --p-front ^CAGCAGCCGCGGTAATTCC \
--p-discard-untrimmed \
--p-no-indels --p-error-rate 0.2 \
--o-trimmed-sequences "int_files/${filename_no_ext}-trimmed.qza"

# visualize read quality following trimming.
qiime demux summarize \
  --i-data "int_files/${filename_no_ext}-trimmed.qza" \
  --o-visualization "run_stats/${filename_no_ext}-demux.qzv"

# the first time we run this we want to stop to look at read quality, 
# which we need for later.
if [ es == true ]
  then
    exit(1)

# begin the denoising process (this may take a while)
qiime dada2 denoise-single \
  --i-demultiplexed-seqs "int_files/${filename_no_ext}-trimmed.qza" \
  --p-trunc-len 220 \
  --p-n-threads 8 \
  --o-representative-sequences "int_files/${filename_no_ext}-rep-seqs.qza" \
  --o-table "int_files/${filename_no_ext}-table.qza" \
  --o-denoising-stats "run_stats/${filename_no_ext}-stats.qza"

# create a visual for the denoising statistics
qiime  metadata  tabulate \
--m-input-file "run_stats/${filename_no_ext}-stats.qza" \ 
--o-visualization "run_stats/${filename_no_ext}-stats.qzv"

# run vsearch to get annotations
qiime feature-classifier classify-consensus-vsearch \ 
--i-reference-reads db/pr2-mothur.qza 
--i-reference-taxonomy db/pr2-tax.qza \
--i-query "int_files/${filename_no_ext}-rep-seqs.qza" \ 
--p-threads 8 \
--o-classification "int_files/${filename_no_ext}-taxonomy-pr2.qza"
