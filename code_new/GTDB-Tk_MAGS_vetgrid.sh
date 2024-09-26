#!/bin/bash
# This script assesses the taxonomy and completeness of the MAGs binned from the metagenomes.
# It requires connexion to the VetLinux server since my local computer does not support CheckM2 nor GTDB-TK2.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data"

### COMMANDS
IFS="
"

ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

if [[ ! -d ~/Bosco/Ancestral_microbiome/GTDB-Tk/MAGs ]]
then  
        mkdir -p ~/Bosco/Ancestral_microbiome/GTDB-Tk/MAGs
fi

FOO

# Copy the genomes to the server
rsync -av \
    "$WORKDIR"/DeNovo_Assembly_*/06_MERGED/SUMMARY_FINAL/bin_by_bin/*/*.fa \
    vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Ancestral_microbiome/GTDB-Tk/MAGs

# Bins cleaned from contamination are in fa format. I change the extension to fasta
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO
POOL=$POOL
for i in ~/Bosco/Ancestral_microbiome/GTDB-Tk/MAGs/*.fa
do
    mv "\${i}" "\${i%.fa}.fasta"
done
FOO

# We connect to the server again to run GTDB-TK2
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

cd ~/Bosco/Ancestral_microbiome/GTDB-Tk

eval \$(conda shell.bash hook)
conda activate gtdbtk-2.1.1

export GTDBTK_DATA_PATH="/home/vetlinux05/Bosco/db/gtdbtk_r214_database"

gtdbtk classify_wf \
        --genome_dir MAGs \
        -x fasta \
        --out_dir output

FOO

# I copy the results and the genomes back to my computer
rsync -av \
    vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Ancestral_microbiome/GTDB-Tk \
    "$WORKDIR"
