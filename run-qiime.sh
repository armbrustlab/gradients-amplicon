#!/bin/bash

# run-qiime.sh START_ARTIFACT QC PROC

# processing pipeline for 18S analysis in qiime2
# read each dataset into qiime2 artifact format using a
# manifest.

file=$1
filename=$(basename -- "$file")
filename_no_ext="${filename%.*}"

qc=$2
process=$3

if [qc == true]
then
qiime cutadapt trim-single --i-demultiplexed-sequences "./int_files/${filename}"  --p-cores 8 \
--p-adapter AGACCAAGTCTCTGC CAAGCAGAAGACGGCATACGAGAT CCCGTGTTGAGTCAAATTAAGC  --p-front CAGCAGCCGCGGTAATTCC ACACTGACGACATGGTTCTACA AATGATACGGCGACCACCGAGATCT \
--p-discard-untrimmed  --p-no-indels --p-error-rate 0.2  --o-trimmed-sequences "./int_files/${filename_no_ext}-trimmed.qza"

qiime demux summarize --i-data "./int_files/${filename_no_ext}-trimmed.qza" --o-visualization "int_files/viz/${filename_no_ext}-demux.qzv"
fi

if [ process == true]
then
nohup qiime dada2 denoise-single --i-demultiplexed-seqs "./int_files/${filename_no_ext}-trimmed.qza" --p-trunc-len 200 --p-n-threads 8 \ 
--o-representative-sequences "./int_files/${filename_no_ext}-rep-seqs.qza" \
--o-table ./int_files/g1-18S-base-artifact-table.qza --o-denoising-stats "./int_files/${filename_no_ext}-denoising-stats.qza" &

qiime metadata tabulate --m-input-file "./int_files/${filename_no_ext}-denoising-stats.qza" --o-visualization "int_files/viz/${filename_no_ext}-denoising-stats.qzv"

nohup qiime feature-classifier classify-consensus-vsearch --i-reference-reads ../db/pr2/pr2-mothur.qza \ 
--i-reference-taxonomy ../db/pr2/pr2-tax.qza --i-query "./int_files/${filename_no_ext}-rep-seqs.qza" --p-threads 8 \ 
--o-classification "./int_files/${filename_no_ext}-taxonomy.qza" &
fi
