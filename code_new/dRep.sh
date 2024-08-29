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

conda deactivate

eval "$(conda shell.bash hook)"
conda activate gtdbtk-2.0
export GTDBTK_DATA_PATH=/Users/bgracia/PhD_local/db/GTDB/release207

mkdir "$TEMP"/GTDB-Tk

for i in $(basename "$TEMP"/dRep/dereplicated_genomes/*.fasta);
do  mkdir "$TEMP"/GTDB-Tk/${i%.fasta};
    cp "$TEMP"/dRep/dereplicated_genomes/${i} "$TEMP"/GTDB-Tk/${i%.fasta}/${i}
    gtdbtk classify_wf \
                --genome_dir "$TEMP"/GTDB-Tk/${i%.fasta} \
                -x fasta \
                --cpus 8 \
                --out_dir "$TEMP"/GTDB-Tk/${i%.fasta};
    rm "$TEMP"/GTDB-Tk/${i%.fasta}/${i};
done

cat "$TEMP"/GTDB-Tk/*/gtdbtk.bac120.summary.tsv | sort -ur > "$TEMP"/GTDB-Tk/summary.tsv

mv "$TEMP"/dRep "$WORKDIR"
mv "$TEMP"/GTDB-Tk/summary.tsv "$WORKDIR"/dRep
rm -r "$TEMP"

# We will create a name2taxon file with the lowest taxonomic level available
cut -f1 "$WORKDIR"/dRep/summary.tsv > "$WORKDIR"/dRep/genome.tmp
cut -f2 "$WORKDIR"/dRep/summary.tsv > "$WORKDIR"/dRep/taxon.tmp

# This removes the unassigned taxonomic levels in each row
for i in s g f o c p
do
    sed -i "" "s/;${i}__\$//g" "$WORKDIR"/dRep/taxon.tmp
done

rev "$WORKDIR"/dRep/taxon.tmp | cut -d ";" -f1 | rev > "$WORKDIR"/dRep/taxon2.tmp

awk '
{
    count[$0]++
    if (count[$0] > 1) {
        print $0 "_" (count[$0])
    } else {
        print $0
    }
}' "$WORKDIR"/dRep/taxon2.tmp | sed "s/ /\_/g" > "$WORKDIR"/dRep/taxon3.tmp
paste "$WORKDIR"/dRep/genome.tmp "$WORKDIR"/dRep/taxon3.tmp > "$WORKDIR"/dRep/name2taxon.tsv
rm "$WORKDIR"/dRep/*.tmp
