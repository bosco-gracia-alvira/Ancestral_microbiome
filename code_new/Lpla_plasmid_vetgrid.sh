#!/bin/bash
# This script maps all the Lpla genomes' readsets against the colonisation plasmid
# Bosco Gracia Alvira, 2024

### VARIABLES
# Set the paths
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Lpla_plasmid"
LOCATION_COLD="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"
VISUALS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/Lpla_plasmid"
LOCATION_ISOLATES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"
LOCATION_SNPS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Lactiplantibacillus_plantarum/Phylogeny/roary"
TAXON2SAMPLE="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/taxonomy.tsv"
RAW_READS="$WORKDIR/reads"
REFERENCE="$WORKDIR/reference"
LOGS="$WORKDIR/logs"
MAPPED="$WORKDIR/mapped"
TREE="$WORKDIR/tree"

# Create the folder where plots will go
if [[ ! -f "$VISUALS" ]]
then
    mkdir -p "$VISUALS"
fi

# Create subfolders in the working directory
if [[ ! -d "$RAW_READS" ]]
then
  mkdir -p "$RAW_READS"
fi

if [[ ! -d "$REFERENCE" ]]
then
  mkdir -p "$REFERENCE"
fi

if [[ ! -d "$MAPPED" ]]
then
  mkdir -p "$MAPPED"
fi

if [[ ! -d "$TREE" ]]
then
  mkdir -p "$TREE"
fi

if [[ ! -d "$LOGS" ]]
then
  mkdir "$LOGS"
fi

# Link the reads from the pools to the raw reads folder
for i in $(basename "$LOCATION_HOT"/F*)
do
  for j in $(seq 1 10)
  do
    # Only link the pool if the file exists
    if [[ -f "$LOCATION_HOT"/${i}/${i}_${j}_noCont_1.fq.gz ]]
    then
      ln -fs "$LOCATION_HOT"/${i}/${i}_${j}_noCont_1.fq.gz "$RAW_READS"/h${i}_${j}_1.fq.gz
      ln -fs "$LOCATION_HOT"/${i}/${i}_${j}_noCont_2.fq.gz "$RAW_READS"/h${i}_${j}_2.fq.gz
    fi
  done
done

for i in $(basename "$LOCATION_COLD"/F*)
do
  for j in $(seq 11 20)
  do
    # Only link the pool if the file exists
    if [[ -f "$LOCATION_COLD"/${i}/${i}_${j}_noCont_1.fq.gz ]]
    then
      ln -fs "$LOCATION_COLD"/${i}/${i}_${j}_noCont_1.fq.gz "$RAW_READS"/c${i}_$(($j-10))_1.fq.gz
      ln -fs "$LOCATION_COLD"/${i}/${i}_${j}_noCont_2.fq.gz "$RAW_READS"/c${i}_$(($j-10))_2.fq.gz
    fi
  done
done

# Select the L. plantarum isolates
cut -f1,2 "$TAXON2SAMPLE" | grep "s__Lactiplantibacillus plantarum" | cut -f1 > "$WORKDIR/isolates.tsv"

# Link the reads from each isolate to the the raw reads folder
while IFS=$'\t' read -r sample
do
    isolate=${sample%_?-contigs}
    r1=$(ls "$LOCATION_ISOLATES"/Pool_???/02.Rm_adapters/fastq_clean/*.clean_1.fq.gz | grep -v '/Pool_589/' | grep "${isolate}.clean_1.fq.gz")
    r2="${r1%1.fq.gz}2.fq.gz"

    ln -fs "${r1}" "${RAW_READS}/${isolate}_1.fq.gz"
    ln -fs "${r2}" "${RAW_READS}/${isolate}_2.fq.gz"

done < "$WORKDIR/isolates.tsv"

# Link the short reads from LpWF to the raw reads folder
ln -fs "${REFERENCE}/SRR28557241_1.fastq.gz" "${RAW_READS}/LpWF_1.fq.gz"
ln -fs "${REFERENCE}/SRR28557241_2.fastq.gz" "${RAW_READS}/LpWF_2.fq.gz"

# Prepare the reference
bowtie2-build --threads 16 "$REFERENCE"/LplaWF.fa "$REFERENCE"/LplaWF

# Map each isolate against the reference
numsamples=$(basename -a "$RAW_READS"/*_1.fq.gz | wc -l | bc)
processed=1
echo -e "Starting mapping of ${numsamples} samples"

# Remove the temporary files
rm "$WORKDIR"/*.t*mp

for i in $(basename -a "$RAW_READS"/*_1.fq.gz)
do

  # If the sample is an isolates genome (doesn't start with ?F), keep only the first field, else keep the two first fields
  if [[ "$i" != ?F* ]]
  then
    name=$(echo "$i" | cut -d "_" -f1)
  else
    name=$(echo "$i" | cut -d "_" -f1,2)
  fi

  echo -e "Mapping sample ${name} (${processed}/${numsamples})"

  # Map paired end reads using bowtie with stringent settings and output the result to a sam file
  bowtie2 \
    -x "$REFERENCE"/LplaWF \
    -q --very-sensitive \
    --no-mixed \
    --no-discordant \
    -1 "$RAW_READS"/${name}_1.fq.gz \
    -2 "$RAW_READS"/${name}_2.fq.gz \
    -S "$MAPPED"/${name}.sam \
    --threads 16 \
    --rg-id "${name}" \
    --rg "SM:${name}" > "$LOGS"/bowtie2_${name}.log 2>&1

  # Turn the sam into bam to save memory, and sort it
  samtools view \
    -bS \
    -@ 16 \
    "$MAPPED"/${name}.sam |\
  samtools sort \
    -@ 16 \
    -o "$MAPPED/${name}_sorted.bam"

  # Delete the sam
  rm "$MAPPED"/${name}.sam

  # Index the bam
  samtools index "$MAPPED/${name}_sorted.bam"

  # Extract the coverage of the plasmid and the housekeeping gene recA
  samtools depth -r contig_2 "$MAPPED/${name}_sorted.bam" | awk -v name="$name" 'BEGIN {OFS="\t"} {print $0, name, "plasmid"}' > "$WORKDIR/${name}_cov_plasmid.tmp"
  samtools depth -r contig_3:2097241-2098383 "$MAPPED/${name}_sorted.bam" | awk -v name="$name" 'BEGIN {OFS="\t"} {print $0, name, "recA"}' > "$WORKDIR/${name}_cov_recA.tmp"

  # Next iteration
  processed=$((processed+1))
done

# Create the header and the body of the coverage table
echo -e "contig\tposition\tcoverage\tsample\tgene" > "$WORKDIR/header.temp"
cat "$WORKDIR"/*.tmp > "$WORKDIR/body.temp"
cat "$WORKDIR/header.temp" "$WORKDIR/body.temp" > "$WORKDIR/coverage.tsv"

# Remove the temporary files
rm "$WORKDIR"/*.t*mp

# Build a phylogenetic tree of the isolates
conda activate

rm "$TREE/accessions.txt"
touch "$TREE/accessions.txt"

# Create a list of all the accessions
for i in $(basename -a "$RAW_READS"/*_1.fq.gz)
do
  name=${i%_1.fq.gz}
  if [[ "$name" != ?F* ]]
  then
    echo "$name" >> "$TREE/accessions.txt"
  fi
done

# Extract the accessions of interest from the MSA file from the phylogeny
seqkit grep -f "$TREE/accessions.txt" \
            -p \
            -i \
            -o "$TREE/Lpla.aln" \
            "$LOCATION_SNPS/core_gene_alignment.aln"

# Run snp-sites to extract only the SNPs
snp-sites \
    -mvp \
    -o "$TREE/Lpla" \
    "$TREE/Lpla.aln"

# Run the tree
iqtree \
    -s "$TREE/Lpla.phylip" \
    --boot 1000 \
    -m GTR+ASC \
    -T 16
