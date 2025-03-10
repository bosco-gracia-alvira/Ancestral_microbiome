#!/bin/bash
# This script maps competitively the reads from all the poolseq replicates and isolate genomes against the collection of de-replicated genomes.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
LOCATION_COLD="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"
LOCATION_ISOLATES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"
LOCATION_REFERENCES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Graph_pangenome"
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome_replicates"
VISUALS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/Competitive_mapping_microbiome_replicates"
RAW_READS="$WORKDIR/reads"
GENOMES="$WORKDIR/genomes"
LOGS="$WORKDIR/logs"
MAPPED="$WORKDIR/mapped"

### COMMANDS

# Print the current shell
echo "Current shell: $SHELL"

# Change IFS to just new line
IFS=$'\n'

# Create subfolders in the working directory
if [[ ! -d "$MAPPED" ]]
then
  mkdir -p "$MAPPED"
fi

if [[ ! -d "$LOGS" ]]
then
  mkdir "$LOGS"
fi

# Create the folder where plots will go
if [[ ! -f "$VISUALS" ]]
then
    mkdir -p "$VISUALS"
fi


# Create the genomes folder and set the path to it
if [[ ! -d "$GENOMES" ]]
then
  mkdir "$GENOMES"
fi

# Copying the genomes to the genomes folder using the name2taxon table
echo -e "Copying the reference pangenomes to the genomes folder"
for i in $(basename "$LOCATION_REFERENCES"/*__*)
do
  echo "${i}"
  cat "$LOCATION_REFERENCES/${i}/${i}.fasta" |\
    seqkit seq -m 2000 |\
    seqkit replace -p .+ -r "${i}_{nr}" --nr-width 3 > "$GENOMES/${i}.fasta"
done

# Combine all the genomes
echo -e "Combining all the pangenomes"
cat "$GENOMES"/*.fasta > "$GENOMES"/combined.fa

# We create the raw reads folder and link the poolseqs to it
echo -e "Linking the reads to the raw reads folder"
if [[ ! -d "$RAW_READS" ]]
then
  mkdir "$RAW_READS"
fi

for i in $(basename "$LOCATION_HOT"/F*)
do
  for j in $(seq 1 10)
  do
    # Only link the pool if the file exists
    if [[ -f "$LOCATION_HOT"/${i}/${i}_${j}_noCont_1.fq.gz ]]
    then
      ln -s "$LOCATION_HOT"/${i}/${i}_${j}_noCont_1.fq.gz "$RAW_READS"/h${i}_${j}_1.fq.gz
      ln -s "$LOCATION_HOT"/${i}/${i}_${j}_noCont_2.fq.gz "$RAW_READS"/h${i}_${j}_2.fq.gz
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
      ln -s "$LOCATION_COLD"/${i}/${i}_${j}_noCont_1.fq.gz "$RAW_READS"/c${i}_$(($j-10))_1.fq.gz
      ln -s "$LOCATION_COLD"/${i}/${i}_${j}_noCont_2.fq.gz "$RAW_READS"/c${i}_$(($j-10))_2.fq.gz
    fi
  done
done

# Create a file linking all the isolate genomes and the location of their reads
# THIS DOES NOT WORK IN THE SCRIPT, JUST IF YOU RUN IT IN THE TERMINAL
for i in "$LOCATION_ISOLATES"/Pool_???/02.Rm_adapters/fastq_clean/*.clean_1.fq.gz
do
  pool=$(echo "${i}" | cut -d "/" -f9 | cut -d "_" -f2)
  if [[ $pool == 589 ]] # We exclude pool 589
  then
    continue
  fi
  sample=$(echo "${i}" | cut -d "/" -f12 | cut -d "." -f1)
  ln -sf "${i}" "$RAW_READS"/i${sample}_1.fq.gz
  ln -sf "${i%_1.fq.gz}_2.fq.gz" "$RAW_READS"/i${sample}_2.fq.gz
done

### Competitive reads mapping to the representative microbiome and extraction of the unmapped reads

# Create an index for the reference combined genome
bowtie2-build --threads 16 "$GENOMES"/combined.fa "$GENOMES"/combined

# Competitive mapping against each of the reads sets
numsamples=$(basename -a "$RAW_READS"/*_1.fq.gz | wc -l | sed "s/ //g")
processed=1
echo -e "Starting competitive mapping of ${numsamples} samples"

for i in $(basename -a "$RAW_READS"/*_1.fq.gz)
do
  
  # If the sample is an isolates genome (starts with "i"), keep only the first field, else keep the two first fields
  if [[ "$i" == i* ]]
  then
    name=$(echo "$i" | cut -d "_" -f1)
  else
    name=$(echo "$i" | cut -d "_" -f1,2)
  fi

  echo -e "Mapping sample ${name} (${processed}/${numsamples})"

  # Create the results folder
  mkdir -p "$MAPPED"/${name}

  # Map paired end reads using bowtie with stringent settings and output the result to a sam file
  bowtie2 \
    -x "$GENOMES"/combined \
    -q \
    -D 500 -R 40 -N 0 -L 20 -i S,1,0.50 \
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
  rm "$MAPPED"/${name}/header.sam

  processed=$((processed+1))
done

# Compute the mappability of the genomes using GenMap (installed in base conda)
conda activate
genmap index -FD "$GENOMES" -I "$GENOMES"/genmap
genmap map -I "$GENOMES"/genmap -O "$GENOMES"/genmap -K 100 -E 0 -bg -t

# Filter the bed files with "low mappability". I set the threshold in 0, but can be adjusted
for beds in $(basename "$GENOMES"/genmap/*.bedgraph)
do
  awk '{if ($4+0 <= 0) print $0}' "$GENOMES"/genmap/${beds} > "$GENOMES"/genmap/${beds%genmap.bedgraph}bed
done

# This chunk calculates the statistics that we are intersted in:
# Number of reads mapped to each genome
# Number of reads mapped UNIQUELY to each genome

echo -e "Calculating reads mapped to each genome"

# Remove any previous temporary file
rm -r "$WORKDIR"/*.col

# For each of the original genomes we create a temporary file to store the number of reads mapped to it, as well as the sample and genome names and sizes
echo -e "Sample" > "$WORKDIR"/sample_name.col
echo -e "Genome" > "$WORKDIR"/genome_name.col
echo -e "Size" > "$WORKDIR"/genome_size.col

for j in "$GENOMES"/*.fasta
do

  # Species variable is the species name
  species=$(basename -a "${j}" | cut -d "." -f1)
  # Genome variable is the fasta file for each species
  genome=$(basename -a "${j}")
  # Create columns for reads mapping to each species: reads mapped, uniq reads mapped, median reads mapped
  echo -e "${species}" > "$WORKDIR"/${species}_reads.col
  echo -e "${species}" > "$WORKDIR"/${species}_uniq.col
  echo -e "${species}" > "$WORKDIR"/${j}_MedCov.col
  # Create a column with species names
  echo -e "${species}" >> "$WORKDIR"/genome_name.col
  # Create a column with species genome size
  seqkit stats "$GENOMES"/"${genome}" | awk 'NR==2 {print $7}' >> "$WORKDIR"/genome_size.col

done

# Reset the number of samples to process
processed=1

# For each sample...
for i in $(basename -a "$MAPPED"/*)
do

  # Define the variables
  sample=${i}
  echo -e "Extracting reads from sample ${sample} (${processed}/${numsamples})"
  processed=$((processed+1))

  # Add the sample name to the file "sample_name.tmp"
  echo ${sample} >> "$WORKDIR"/sample_name.col

  # For each of the original genomes...
  for j in $(basename -a "$GENOMES"/*.fasta | cut -d "." -f1)
  do
    # Species variable is the species name
    species=$(basename -a "${j}" | cut -d "." -f1)
    # Filter the bam file based on the mappability bed file
    bedtools intersect -v -abam "$MAPPED"/${sample}/${species}.bam -b "$GENOMES"/genmap/combined.bed > "$MAPPED"/${sample}/${species}_filt.bam
    # Extract the number of reads mapped to the genome and add it to ""$WORKDIR"/${j}_reads.tmp"
    samtools view -c -F 4 "$MAPPED"/${sample}/${species}_filt.bam >> "$WORKDIR"/${species}_reads.col
    # Extract the number of reads mapped uniquely to the genome (with MAPQ>3) and add it to ""$WORKDIR"/${j}_uniq.tmp"
    samtools view -c -F 4 -q 4 "$MAPPED"/${sample}/${species}_filt.bam >> "$WORKDIR"/${species}_uniq.col
    # Extract the median coverage of the species in that sample
    samtools depth -a "$MAPPED"/${sample}/${species}_filt.bam | \
      awk '{print $3}' | sort -n | awk '{
        count[NR] = $1;
        }
        END {
            if (NR % 2) {
                median = count[(NR + 1) / 2];
            } else {
                median = (count[NR / 2] + count[NR / 2 + 1]) / 2;
            }
            print median;
        }' >> "$WORKDIR"/${species}_MedCov.col
  done 
  
done

paste "$WORKDIR"/genome_name.col "$WORKDIR"/genome_size.col > "$WORKDIR"/genome_size.tsv
paste "$WORKDIR"/sample_name.col "$WORKDIR"/*_reads.col > "$WORKDIR"/reads_mapped.tsv
paste "$WORKDIR"/sample_name.col "$WORKDIR"/*_uniq.col > "$WORKDIR"/uniq_mapped.tsv
paste "$WORKDIR"/sample_name.col "$WORKDIR"/*_MedCov.col > "$WORKDIR"/MedCov.tsv

rm -r "$WORKDIR"/*.col

echo -e "Mapping finished"
