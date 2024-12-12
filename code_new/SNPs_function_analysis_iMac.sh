#!/bin/bash
# This script annotates the vcf files obtained by mapping each poolseq and isolate to the graph pangenomes.
# Bosco Gracia Alvira, 2024

### VARIABLES
LOCATION_REFERENCES="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome/genomes"

LOCATION_SNPS="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/SNPs_analysis"

### COMMANDS

# Annotate each graph pangenome with Bakta
conda activate bakta
export BAKTA_DB="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/db/Bakta/db-light"

for i in $(basename "$LOCATION_REFERENCES"/*__*.fasta)
do
    genome=$(basename "${i%.fasta}")

    echo "Annotating $genome with Bakta"

    bakta \
        -p "${genome}.annotated" \
        -o "${LOCATION_SNPS}/${genome}/bakta/" \
        --species "${genome}" \
        --threads 8 \
        --keep-contig-headers \
        "${LOCATION_REFERENCES}/${genome}.fasta"

done

# Integrate the annotations into the vcf files mapped to the same reference pangenomes using snpeff
for i in $(basename "$LOCATION_REFERENCES"/*__*.fasta)
do
    genome=$(basename "${i%.fasta}")

    vcf-annotator \
        --output "${LOCATION_SNPS}/${genome}/${genome}_annotated.vcf" \
        "${LOCATION_SNPS}/${genome}/${genome}.vcf" \
        "${LOCATION_SNPS}/${genome}/bakta/${genome}.annotated.gbff"
done
