#!/bin/bash
# 
# Bosco Gracia Alvira, 2025

### VARIABLES
# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/metaphlan"
LOCATION_COLD="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"
DB="/Users/bgracia/PhD_local/db/metaphlan"
VISUALS="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/metaphlan"

# Create the working directory and visuals folders
if [[ ! -f "$WORKDIR" ]]
then
    mkdir -p "$WORKDIR"
fi

if [[ ! -f "$VISUALS" ]]
then
    mkdir -p "$VISUALS"
fi

conda activate metaphlan-4.1

for i in "$LOCATION_HOT"/F*
do
    pool=$(basename "$i")
    metaphlan \
        "${i}/Pooled_${pool}_noCont_1.fq.gz","${i}/Pooled_${pool}_noCont_2.fq.gz" \
        --bowtie2db "$DB" \
        -t rel_ab_w_read_stats \
        --nproc 12 \
        --input_type fastq \
        --bowtie2out "$WORKDIR/HOT_${pool}.bowtie2.bz2" \
        -o "$WORKDIR/HOT_${pool}_profile.txt"
        #-s "$WORKDIR/HOT_${pool}.sam.bz2" \
done

for i in "$LOCATION_COLD"/F*
do
    pool=$(basename "$i")
    metaphlan \
        "${i}/Pooled_${pool}_noCont_1.fq.gz","${i}/Pooled_${pool}_noCont_2.fq.gz" \
        --bowtie2db "$DB" \
        -t rel_ab_w_read_stats \
        --nproc 12 \
        --input_type fastq \
        --bowtie2out "$WORKDIR/COLD_${pool}.bowtie2.bz2" \
        -o "$WORKDIR/COLD_${pool}_profile.txt"
        #-s "$WORKDIR/COLD_${pool}.sam.bz2" \
done

# Merge the tables
merge_metaphlan_tables.py "$WORKDIR"/*_profile.txt > "$WORKDIR/merged_abundance_table.txt"
