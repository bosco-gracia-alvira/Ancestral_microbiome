#!/bin/bash
# This script creates a graph pangenome from all the isolates assigned to one species
# Bosco Gracia Alvira, 2024

### VARIABLES
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data"
ISOLATES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Isolates"
MAGS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/GTDB-Tk/MAGs"

### COMMANDS
IFS=$'\n'

# If the Graph_pangenome does not exist, create it
if [[ ! -d "$WORKDIR"/Graph_pangenome/Genomes ]]
then  
    mkdir -p "$WORKDIR"/Graph_pangenome/Genomes
fi

# Link the genomes to the Genomes directory
cp "$ISOLATES"/*.fasta "$WORKDIR"/Graph_pangenome/Genomes
cp "$MAGS"/*MAG*.fasta "$WORKDIR"/Graph_pangenome/Genomes

# Compile the taxonomy of all the available good quality taxa (not the metagenomic bins, only the MAGs)
cat "/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/taxonomy.tsv" > "$WORKDIR"/Graph_pangenome/taxonomy.tsv
tail -n +2 "$WORKDIR"/GTDB-Tk/output/gtdbtk.bac120.summary.tsv |\
  grep -v "_Bin_" >> "$WORKDIR"/Graph_pangenome/taxonomy.tsv
TAXONOMY="$WORKDIR/Graph_pangenome/taxonomy.tsv"

# List the number of taxa available
TAXON_LIST=$(cut -f2 "$TAXONOMY" | awk -F 's__' '{print $2}' | grep " " | sed 's/ /_/' | sort | uniq -c | sort -r | sed 's/^ *//')

echo "$TAXON_LIST" | while read -r i
do
  IFS=$'\n'
  # Number of genomes available for the given species
  count=$(echo "$i" | cut -d " " -f 1)
  # Species name
  sp=$(echo "$i" | cut -d " " -f 2)
  # Species name with space between genus and species
  SPECIES=$(echo "$sp" | sed 's/_/ /')
  # Name of the samples belonging to that species
  SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' "$TAXONOMY" | grep -v "user")

  if [ "$count" -eq 1 ]
  then

    echo "There is only one genome available for the species $sp. Using this lonely genome as species reference."

    mkdir -p "$WORKDIR/Graph_pangenome/$sp"
    rsync -av \
          "$WORKDIR/Graph_pangenome/Genomes/$SAMPLES.fasta" \
          "$WORKDIR/Graph_pangenome/$sp/$sp.fasta"

  else

    echo "There are $count genomes available for the species $sp. Thus, we will make a pangenome."

    # Create the folder for the species
    ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

    if [[ ! -d ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp"/Genomes ]]
    then  
        mkdir -p ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp"/Genomes
    fi

FOO

    # Copy the genomes to the VetLinux server
    for j in $SAMPLES
    do
            rsync -av \
              "$WORKDIR/Graph_pangenome/Genomes/$j.fasta" \
              vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp"/Genomes
    done

    # In this long EOF I run the whole SuperPang pipeline remotely in VetLinux
    ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

    eval "\$(conda shell.bash hook)"
    conda activate SuperPang-0.9

    cd ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp"

    echo "Making pangenome for $sp"

    SuperPang.py \
            --fasta ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp"/Genomes/*.fasta \
            --assume-complete \
            --output-dir ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp" \
            --force-overwrite \
            --threads 24

    mv \
      ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp"/assembly.fasta \
      ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp"/"$sp".fasta

    rm -r ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp"/Genomes

FOO

  # Copy back the results to the Graph_pangenome folder
  rsync -av \
    vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Ancestral_microbiome/Graph_pangenome/"$sp" \
    "$WORKDIR"/Graph_pangenome

  fi

done
