#!/bin/bash
# This script maps competitively the reads from all the poolseqs against the collection of de-replicated genomes.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
LOCATION_COLD="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"
LOCATION_ISOLATES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/dRep/dereplicated_genomes"

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
while read -r name taxon
do
  cat "$LOCATION_ISOLATES/${name}.fasta" |\
    seqkit seq -m 2000 |\
    seqkit replace -p .+ -r "${taxon}_{nr}" --nr-width 3 > "$GENOMES/${taxon}.fasta"
done < <(tail -n +2 "$LOCATION_ISOLATES/../name2taxon.tsv")

# Combine all the genomes
cat "$GENOMES"/*.fasta > "$GENOMES"/combined.fa

# We create the raw reads folder and link the poolseqs to it
if [[ ! -d "$RAW_READS" ]]
then
  mkdir "$RAW_READS"
fi

for i in $(basename "$LOCATION_HOT"/F*)
do
  ln -s "$LOCATION_HOT"/F*/Pooled_${i}_Clean_noCont_1.fq.gz "$RAW_READS"/h${i}_1.fq.gz
  ln -s "$LOCATION_HOT"/F*/Pooled_${i}_Clean_noCont_2.fq.gz "$RAW_READS"/h${i}_2.fq.gz
done

for i in $(basename "$LOCATION_COLD"/F*)
do
  ln -s "$LOCATION_COLD"/F*/Pooled_${i}_Clean_noCont_1.fq.gz "$RAW_READS"/c${i}_1.fq.gz
  ln -s "$LOCATION_COLD"/F*/Pooled_${i}_Clean_noCont_2.fq.gz "$RAW_READS"/c${i}_2.fq.gz
done


### Competetive reads mapping to the representative microbiome and extraction of the unmapped reads

# Create an index for the reference combined genome
bowtie2-build --threads 16 "$GENOMES"/combined.fa "$GENOMES"/combined

# Competitive mapping mapping against each of the reads sets
for i in $(basename -a "$RAW_READS"/*_1.fq.gz)
do

  name=$(echo "$i" | cut -d "_" -f1,2)

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
rm "$WORKDIR"/*.col

# For each of the original genomes we create a temporary file to store the number of reads mapped to it, as well as the sample names
echo -e "Sample" > "$WORKDIR"/sample_name.col
for j in $(basename -a "$GENOMES"/*.fasta | cut -d "." -f1)
do
  echo -e "${j}" > "$WORKDIR"/${j}_reads.col
  echo -e "${j}" > "$WORKDIR"/${j}_uniq.col
done

# For each sample...
for i in $(basename -a "$RAW_READS"/*_1.fq.gz | cut -d "_" -f1,2)
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

paste "$WORKDIR"/sample_name.col "$WORKDIR"/*_reads.col > "$WORKDIR"/reads_mapped.tsv
paste "$WORKDIR"/sample_name.col "$WORKDIR"/*_uniq.col > "$WORKDIR"/uniq_mapped.tsv

rm "$WORKDIR"/*.col


 








### Make the compositional table per generation
declare x=1
for i in $(basename "$WORKDIR"/logs/*_depth.txt | cut -d "_" -f1)
do
    sed "s/jTAD/${i}/" "$WORKDIR"/logs/${i}_depth.txt | sort -rk 1 > "$WORKDIR"/depth_${x}.tmp
    ((x++))
done

output_file="DepthGen.tsv"
cp "$WORKDIR"/depth_1.tmp "$WORKDIR"/"$output_file"

declare z=$(basename "$WORKDIR"/logs/*_depth.txt | wc -l)

for ((i=2; i<=${z}; i++))
do
  join -t $'\t' -1 1 -2 1 "$WORKDIR"/"$output_file" "$WORKDIR"/"depth_${i}.tmp" > "$WORKDIR"/"tmp_file.txt"
  mv "$WORKDIR"/"tmp_file.txt" "$WORKDIR"/"$output_file"
done

# Make the ANIr table
for i in $(basename "$WORKDIR"/logs/*_depth.txt | cut -d "_" -f1)
do
  for z in $(basename "$GENOMES"/*.fasta | cut -d "." -f1)
  do
    anir.rb --m-format bam -m "$WORKDIR"/"${i}"/"${z}".sorted.bam -a fix --identity 95 -l ${i}.tmp
  done
  grep "ANIr:" ${i}.tmp | rev | cut -d " " -f1 | rev > "${i}".txt
  echo "${i}" > "$WORKDIR"/0000.tmp
  cat *.tmp > "$WORKDIR"/"${i}".col
  rm *.tmp
done
for i in $(basename "$WORKDIR"/logs/*_depth.txt | cut -d "_" -f1)
do

for z in $(basename "$GENOMES"/*.fasta | cut -d "." -f1)
do
  echo "${z}" >> "$WORKDIR"/0000.txt
done

paste *.txt > "$WORKDIR"/ANIr.tsv


### Record the number of reads per generation
for i in $(basename "$WORKDIR"/logs/*_depth.txt | cut -d "_" -f1)
do
    echo "${i}" >> "$WORKDIR"/total_header.tmp
    grep "reads; of these:" "$WORKDIR"/logs/bowtie2_${i}.txt | cut -d " " -f1 >> "$WORKDIR"/total_reads.tmp
    grep "overall alignment rate" "$WORKDIR"/logs/bowtie2_${i}.txt | cut -d "%" -f1 >> "$WORKDIR"/aligned_reads.tmp
done

paste "$WORKDIR"/total_header.tmp "$WORKDIR"/total_reads.tmp "$WORKDIR"/aligned_reads.tmp > "$WORKDIR"/total_reads.tsv

rm "$WORKDIR"/*.tmp
