#!/bin/bash
# This script annotates the vcf files obtained by mapping each poolseq and isolate to the graph pangenomes.
# Bosco Gracia Alvira, 2024

### VARIABLES
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome"

VISUALS="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/Competitive_mapping_microbiome"

LOCATION_REFERENCES="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Graph_pangenome"

# List the available taxa
TAXON_LIST=$(cut -f2 "$LOCATION_REFERENCES"/name2taxon.tsv | grep "__")


### COMMANDS

# Annotate each graph pangenome with Bakta
conda activate bakta
export BAKTA_DB="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/db/Bakta/db-light"
i="s__Lactiplantibacillus_plantarum"
echo "$TAXON_LIST" | while read -r i
do
    echo "Annotating $i"
    bakta \
        -p "${i}.annotated" \
        -o "${LOCATION_REFERENCES}/${i}"/bakta/ \
        --species "${i}" \
        --threads 16 \
        "$LOCATION_REFERENCES"/"$i"/"$i".fasta
done

# Integrate the annotations into the vcf files mapped to the same reference pangenomes using snpeff

