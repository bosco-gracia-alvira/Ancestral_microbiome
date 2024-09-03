#!/bin/bash
# This script maps competitively the reads from all the poolseqs and isolate genomes against the collection of de-replicated genomes.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
LOCATION_COLD="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"
LOCATION_REFERENCES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/dRep/dereplicated_genomes"
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome"
RAW_READS="$WORKDIR/reads"
GENOMES="$WORKDIR/genomes"
LOGS="$WORKDIR/logs"
MAPPED="$WORKDIR/mapped"

### COMMANDS
IFS="
"

if [[ ! -d "$MAPPED" ]]
then
  mkdir -p "$MAPPED"
fi

if [[ ! -d "$LOGS" ]]
then
  mkdir "$LOGS"
fi

# Create the genomes folder and set the path to it
if [[ ! -d "$GENOMES" ]]
then
  mkdir "$GENOMES"
fi

# Link the genomes to the genomes folder using the name2taxon table
IFS=$'\t'
while read -r name taxon
do
  cat "$LOCATION_REFERENCES/${name}.fasta" |\
    seqkit seq -m 2000 |\
    seqkit replace -p .+ -r "${taxon}_{nr}" --nr-width 3 > "$GENOMES/${taxon}.fasta"
done < <(tail -n +2 "$LOCATION_REFERENCES/../name2taxon.tsv")
unset IFS

# Combine all the genomes
cat "$GENOMES"/*.fasta > "$GENOMES"/combined.fa

# We create the raw reads folder and link the poolseqs to it
if [[ ! -d "$RAW_READS" ]]
then
  mkdir "$RAW_READS"
fi

for i in $(basename "$LOCATION_HOT"/F*)
do
  ln -sf "$LOCATION_HOT"/F*/Pooled_${i}_Clean_noCont_1.fq.gz "$RAW_READS"/h${i}_1.fq.gz
  ln -sf "$LOCATION_HOT"/F*/Pooled_${i}_Clean_noCont_2.fq.gz "$RAW_READS"/h${i}_2.fq.gz
done

for i in $(basename "$LOCATION_COLD"/F*)
do
  ln -sf "$LOCATION_COLD"/F*/Pooled_${i}_Clean_noCont_1.fq.gz "$RAW_READS"/c${i}_1.fq.gz
  ln -sf "$LOCATION_COLD"/F*/Pooled_${i}_Clean_noCont_2.fq.gz "$RAW_READS"/c${i}_2.fq.gz
done

# Create a file linking all the isolate genomes and the location of their reads
for i in $(ls /Volumes/Data/PopGen\ Dropbox/Martin\ McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_???/02.Rm_adapters/fastq_clean/*.clean_1.fq.gz)
do
  sample=$(echo "${i}" | cut -d "/" -f12 | cut -d "." -f1)
  ln -sf "${i}" "$RAW_READS"/i${sample}_1.fq.gz
  ln -sf "${i%_1.fq.gz}_2.fq.gz" "$RAW_READS"/i${sample}_2.fq.gz
done

### Competetive reads mapping to the representative microbiome and extraction of the unmapped reads

# Create an index for the reference combined genome
bowtie2-build --threads 16 "$GENOMES"/combined.fa "$GENOMES"/combined

# Competitive mapping against each of the reads sets
for i in $(basename -a "$RAW_READS"/*_1.fq.gz)
do

  name=$(echo "$i" | cut -d "_" -f1)

  # Create the results folder
  mkdir -p "$MAPPED"/${name}

  # Map paired end reads using bowtie with stringent settings and output the result to a sam file
  bowtie2 \
    -x "$GENOMES"/combined \
    -q --very-sensitive \
    --no-mixed \
    --no-discordant \
    -1 "$RAW_READS"/${name}_1.fq.gz \
    -2 "$RAW_READS"/${name}_2.fq.gz \
    -S "$MAPPED"/${name}/combined.sam \
    --threads 16 \
    --rg-id "${name}" \
    --rg "SM:${name}" > "$LOGS"/bowtie2_${name}.log 2>&1

  # Turn the sam into bam to save memory
  samtools view \
    -bS \
    -@ 16 \
    "$MAPPED"/${name}/combined.sam > "$MAPPED"/${name}/combined.bam

  # Delete the sam
  rm "$MAPPED"/${name}/combined.sam

  # Sort the bam
  samtools sort \
    -@ 16 \
    -O bam \
    -o "$MAPPED"/${name}/combined_sorted.bam \
    "$MAPPED"/${name}/combined.bam

  # Delete the unsorted bam
  rm "$MAPPED"/${name}/combined.bam

  # Extract the header of the bam
  samtools view -H "$MAPPED"/${name}/combined_sorted.bam > "$MAPPED"/${name}/header.sam

  # Index the combined bam
  samtools index "$MAPPED"/${name}/combined_sorted.bam

  # For each of the original genomes...
  for j in $(basename -a "$GENOMES"/*.fasta | cut -d "." -f1)
  do  
    # Make a list with all the contigs belonging to that genome
    contigs=$(samtools idxstats "$MAPPED"/${name}/combined_sorted.bam | cut -f 1 | grep "^${j}")

    # From the combined bam, extract the reads that map the contigs of the specific genome "j"
    # Then, add the header and create a new file with the reads mapping to "j"
    samtools view -@ 16 -b "$MAPPED"/${name}/combined_sorted.bam $contigs |\
      samtools reheader "$MAPPED"/${name}/header.sam - > "$MAPPED"/${name}/${j}.bam
  done

  # Remove the combined bam file, the index and the header
  rm "$MAPPED"/${name}/combined_sorted.bam
  rm "$MAPPED"/${name}/combined_sorted.bam.bai
done

# This chunk calculates the statistics that we are intersted in:
# Number of reads mapped to each genome
# Number of reads mapped UNIQUELY to each genome

# Remove any previous temporary file
rm "$WORKDIR"/*.col

# For each of the original genomes we create a temporary file to store the number of reads mapped to it, as well as the sample and genome names and sizes

echo -e "Sample" > "$WORKDIR"/sample_name.col
echo -e "Genome" > "$WORKDIR"/genome_name.col
echo -e "Size" > "$WORKDIR"/genome_size.col

for j in $(basename -a "$GENOMES"/*.fasta | cut -d "." -f1)
do
  echo -e "${j}" > "$WORKDIR"/${j}_reads.col
  echo -e "${j}" > "$WORKDIR"/${j}_uniq.col
done

for j in $(basename -a "$GENOMES"/*.fasta | cut -d "." -f1)
do
  echo -e "${j}" >> "$WORKDIR"/genome_name.col
  genome="$GENOMES/${j}.fasta"
  seqkit stats "$genome" | awk 'NR==2 {print $7}' >> "$WORKDIR"/genome_size.col
done

# For each sample...
for i in $(basename -a "$RAW_READS"/*_1.fq.gz | cut -d "_" -f1)
do

  # Add the sample name to the file "sample_name.tmp"
  echo ${i} >> "$WORKDIR"/sample_name.col

  # For each of the original genomes...
  for j in $(basename -a "$GENOMES"/*.fasta | cut -d "." -f1)
  do

    # Extract the number of reads mapped to the genome and add it to ""$WORKDIR"/${j}_reads.tmp"
    samtools view -c -F 4 "$MAPPED"/${i}/${j}.bam >> "$WORKDIR"/${j}_reads.col

    # Extract the number of reads mapped uniquely to the genome (with MAPQ>3) and add it to ""$WORKDIR"/${j}_uniq.tmp"
    samtools view -c -F 4 -q 4 "$MAPPED"/${i}/${j}.bam >> "$WORKDIR"/${j}_uniq.col

  done 

done

paste "$WORKDIR"/genome_name.col "$WORKDIR"/genome_size.col > "$WORKDIR"/genome_size.tsv
paste "$WORKDIR"/sample_name.col "$WORKDIR"/*_reads.col > "$WORKDIR"/reads_mapped.tsv
paste "$WORKDIR"/sample_name.col "$WORKDIR"/*_uniq.col > "$WORKDIR"/uniq_mapped.tsv

rm "$WORKDIR"/*.col
