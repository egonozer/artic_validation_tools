#!/bin/bash

if [ $# -lt 2 ]; then
    echo "artic_remask.sh <prefix> <minimum depth>"
    exit
fi

pref=$1
depth=$2

if [ ! -e "$pref.primertrimmed.rg.sorted.bam" ]; then
    echo "ERROR: $pref.primertrimmed.rg.sorted.bam not found!"
    exit
fi

source ~/miniconda3/etc/profile.d/conda.sh    
conda activate artic-ncov2019

artic_make_depth_mask \
    --depth $depth \
    /home/eoz639/Applications/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta \
    $pref.primertrimmed.rg.sorted.bam \
    $pref.coverage_mask.txt
artic_mask \
    /home/eoz639/Applications/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta \
    $pref.coverage_mask.txt \
    $pref.fail.vcf \
    $pref.preconsensus.fasta
bcftools consensus \
    -f $pref.preconsensus.fasta \
    $pref.pass.vcf.gz \
    -m $pref.coverage_mask.txt \
    -o $pref.consensus.fasta
artic_fasta_header $pref.consensus.fasta "$pref/ARTIC/nanopolish/min$depth"

conda deactivate
