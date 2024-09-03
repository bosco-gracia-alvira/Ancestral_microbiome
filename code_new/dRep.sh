#!/bin/bash
# This script collects all the metagenomic and isolate genomes and de-replicates them, keeping only one representative for each taxon.
# https://drep.readthedocs.io/en/latest/overview.html#a-drep-based-metagenomic-workflow
# Bosco Gracia Alvira, 2023

### VARIABLES
WORKDIR="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox/Ancestral_microbiome/data"
ASSEMBLY="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox/Isolates_assembly"
TEMP="/Users/bgracia/PhD_local/temp"

### COMMANDS
if [ ! -d "$TEMP"/Genomes ]
then
    mkdir -p "$TEMP"/Genomes
fi

cp "$WORKDIR"/DeNovo_Assembly_*/06_MERGED/SUMMARY_FINAL/bin_by_bin/*/*.fa "$TEMP"/Genomes/
cp "$ASSEMBLY"/Pool_???/07.GTDB-Tk/Genomes/*.fasta "$TEMP"/Genomes/

for i in $(basename -a "$TEMP"/Genomes/*.fa)
do
    mv "$TEMP"/Genomes/${i} "$TEMP"/Genomes/${i}sta
done

eval "$(conda shell.bash hook)"
conda activate drep-3.4
export CHECKM_DATA_PATH=/Users/bgracia/PhD_local/db/CheckM

# Note that I use a secondary alignment threshold of 98% because I want strain resolution.
dRep dereplicate \
    "$TEMP"/dRep \
    -g "$TEMP"/Genomes/*.fasta \
    --S_algorithm fastANI \
    --multiround_primary_clustering \
    -ms 10000 \
    -pa 0.9 \
    -sa 0.98 \
    -nc 0.30 \
    -cm larger \
    -p 10

mv "$TEMP"/dRep "$WORKDIR"

# Now we want to run GTDB-Tk2 on the de-replicated genomes
# We connect to the server to create a directory and copy the genomes
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

if [[ ! -d ~/Bosco/Ancestral_microbiome/dRep/GTDB-Tk/Genomes ]]
then  
    mkdir -p ~/Bosco/Ancestral_microbiome/dRep/GTDB-Tk/Genomes
fi

FOO

scp -r "$WORKDIR"/dRep/dereplicated_genomes/*.fasta \
        vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Ancestral_microbiome/dRep/GTDB-Tk/Genomes

# We connect to the server again to run GTDB-TK2
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

cd ~/Bosco/Ancestral_microbiome/dRep/GTDB-Tk/

eval \$(conda shell.bash hook)
conda activate gtdbtk-2.1.1

export GTDBTK_DATA_PATH="/home/vetlinux05/Bosco/db/gtdbtk_r214_database"

gtdbtk classify_wf \
        --genome_dir ~/Bosco/Ancestral_microbiome/dRep/GTDB-Tk/Genomes \
        -x fasta \
        --cpus 12 \
        --pplacer_cpus 12 \
        --out_dir ~/Bosco/Ancestral_microbiome/dRep/GTDB-Tk/output

FOO

# Finally we copy the results and the genomes back to my computer
scp -r vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Ancestral_microbiome/dRep/GTDB-Tk/output/* "$WORKDIR"/dRep/GTDB-Tk/

# Now create a name2taxon file with the lowest taxonomic level available for each genome
cut -f1 "$WORKDIR"/dRep/GTDB-Tk/gtdbtk.bac120.summary.tsv > "$WORKDIR"/dRep/genome.tmp
cut -f2 "$WORKDIR"/dRep/GTDB-Tk/gtdbtk.bac120.summary.tsv > "$WORKDIR"/dRep/taxon.tmp

# This loop removes the unassigned taxonomic levels in each row
for i in s g f o c p
do
    sed -i "" "s/;${i}__\$//g" "$WORKDIR"/dRep/taxon.tmp
done

# We only leave the last rank for each genome
rev "$WORKDIR"/dRep/taxon.tmp | cut -d ";" -f1 | rev > "$WORKDIR"/dRep/taxon2.tmp

# Some genomes have the same taxon, so we will add a number to differentiate them
awk '
{
    count[$0]++
    if (count[$0] > 1) {
        print $0 "_" (count[$0])
    } else {
        print $0
    }
}' "$WORKDIR"/dRep/taxon2.tmp | sed "s/ /\_/g" > "$WORKDIR"/dRep/taxon3.tmp

# We paste the genome and taxon files to create the name2taxon file and remove temporary files
paste "$WORKDIR"/dRep/genome.tmp "$WORKDIR"/dRep/taxon3.tmp > "$WORKDIR"/dRep/name2taxon.tsv
rm "$WORKDIR"/dRep/*.tmp
