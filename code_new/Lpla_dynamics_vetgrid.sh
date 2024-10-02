#!/bin/bash
# This script maps competitively the reads from each replicate and generation against the L_plantarum genomes S103 (cold), S239 (hot) and B89 (base).
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
LOCATION_COLD="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"
LOCATION_ISOLATES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"

WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_Lpla_3strains_replicates"
RAW_READS="$WORKDIR/reads"
GENOMES="$WORKDIR/genomes"
LOGS="$WORKDIR/logs"
MAPPED="$WORKDIR/mapped"
METADATA="$WORKDIR/Metadata.csv"

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

# Copy the three competting genomes. Entry names must be in the format S103_0001, S103_0002...
S103="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Isolates/S103.fasta"
ln -s "$S103" "$GENOMES"
S239="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Isolates/S239.fasta"
ln -s "$S239" "$GENOMES"
B89="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Isolates/B89.fasta"
ln -s "$B89" "$GENOMES"

# Combine the three genomes
cat "$GENOMES"/S103.fasta "$GENOMES"/S239.fasta "$GENOMES"/B89.fasta > "$GENOMES"/combined.fa

# We create the raw reads folder and link the poolseqs and isolates reads to it
if [[ ! -d "$RAW_READS" ]]
then
  mkdir "$RAW_READS"
fi

for i in $(basename "$LOCATION_HOT"/F*)
do
  for j in $(seq 1 10)
  do
    ln -s "$LOCATION_HOT"/${i}/${i}_${j}_noCont_1.fq.gz "$RAW_READS"/h${i}_${j}_1.fq.gz
    ln -s "$LOCATION_HOT"/${i}/${i}_${j}_noCont_2.fq.gz "$RAW_READS"/h${i}_${j}_2.fq.gz
  done
done

for i in $(basename "$LOCATION_COLD"/F*)
do
  for j in $(seq 11 20)
  do
    ln -s "$LOCATION_COLD"/${i}/${i}_${j}_noCont_1.fq.gz "$RAW_READS"/c${i}_$(($j-10))_1.fq.gz
    ln -s "$LOCATION_COLD"/${i}/${i}_${j}_noCont_2.fq.gz "$RAW_READS"/c${i}_$(($j-10))_2.fq.gz
  done
done

ln -s "$LOCATION_ISOLATES"/Pool_503/02.Rm_adapters/fastq_clean/S103.clean_1.fq.gz "$RAW_READS"/iS103_1.fq.gz
ln -s "$LOCATION_ISOLATES"/Pool_503/02.Rm_adapters/fastq_clean/S103.clean_2.fq.gz "$RAW_READS"/iS103_2.fq.gz

ln -s "$LOCATION_ISOLATES"/Pool_503/02.Rm_adapters/fastq_clean/S239.clean_1.fq.gz "$RAW_READS"/iS239_1.fq.gz
ln -s "$LOCATION_ISOLATES"/Pool_503/02.Rm_adapters/fastq_clean/S239.clean_2.fq.gz "$RAW_READS"/iS239_2.fq.gz

ln -s "$LOCATION_ISOLATES"/Pool_643/02.Rm_adapters/fastq_clean/B89.clean_1.fq.gz "$RAW_READS"/iB89_1.fq.gz
ln -s "$LOCATION_ISOLATES"/Pool_643/02.Rm_adapters/fastq_clean/B89.clean_2.fq.gz "$RAW_READS"/iB89_2.fq.gz

### Competetive reads mapping to the representative microbiome and extraction of the unmapped reads

# Create an index for the reference combined genome
bowtie2-build --threads 16 "$GENOMES"/combined.fa "$GENOMES"/combined

# Competitive mapping mapping against each of the reads sets
for i in $(basename -a "$RAW_READS"/*_1.fq.gz)
do

  # If the sample is an isolates genome (starts with "i"), keep only the first field, else keep the two first fields
  if [[ "$i" == i* ]]
  then
  name=$(echo "$i" | cut -d "_" -f1)
  else
  name=$(echo "$i" | cut -d "_" -f1,2)
  fi
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
    --threads 16 > "$LOGS"/bowtie2_${name}.log 2>&1

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
rm "$WORKDIR"/*.tmp

# For each sample...
for i in $(basename -a "$RAW_READS"/*_1.fq.gz | cut -d "_" -f1,2)
do

  # Add the sample name to the file "sample_name.tmp"
  echo ${i} >> "$WORKDIR"/sample_name.tmp

  # For each of the original genomes...
  for j in $(basename -a "$GENOMES"/*.fasta | cut -d "." -f1)
  do

    # Extract the number of reads mapped to the genome and add it to ""$WORKDIR"/${j}_reads.tmp"
    samtools view -c -F 4 "$MAPPED"/${i}/${j}.bam >> "$WORKDIR"/${j}_reads.tmp

    # Extract the number of reads mapped uniquely to the genome (with MAPQ>3) and add it to ""$WORKDIR"/${j}_uniq.tmp"
    samtools view -c -F 4 -q 4 "$MAPPED"/${i}/${j}.bam >> "$WORKDIR"/${j}_uniq.tmp

  done 

done

# Paste the columns generated in the previous chunk
paste \
  "$WORKDIR"/sample_name.tmp \
  "$WORKDIR"/S103_reads.tmp \
  "$WORKDIR"/S239_reads.tmp \
  "$WORKDIR"/B89_reads.tmp \
  "$WORKDIR"/S103_uniq.tmp \
  "$WORKDIR"/S239_uniq.tmp \
  "$WORKDIR"/B89_uniq.tmp > "$WORKDIR"/body.tmp

# Create the header with the column names
echo -e "Sample\tS103_reads\tS239_reads\tB89_reads\tS103_uniq\tS239_uniq\tB89_uniq" > "$WORKDIR"/head.tmp

# Bind the header to the body and create a file with the depth of each genome in each sample and the ratio between genomes
#cat "$WORKDIR"/head.tmp "$WORKDIR"/body.tmp > "$WORKDIR"/coverage.tsv

# Remove any previous temporary file
#rm "$WORKDIR"/*.tmp
