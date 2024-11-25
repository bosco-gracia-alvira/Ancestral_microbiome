#!/bin/bash
# This script maps competitively the reads from each replicate and generation against the L_plantarum genomes S103 (cold), S239 (hot) and B89 (base).
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Lpla_dynamics_replicates"
VISUALS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/Lpla_dynamics_replicates"
BAMS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome_replicates/mapped"
NAME2TAXON="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Graph_pangenome/name2taxon.tsv"
RAW_READS="$WORKDIR/reads"
GENOMES="$WORKDIR/genomes"
LOGS="$WORKDIR/logs"
MAPPED="$WORKDIR/mapped"

### COMMANDS
IFS="
"

# Create the folder where plots will go
if [[ ! -f "$VISUALS" ]]
then
    mkdir -p "$VISUALS"
fi

# Create subfolders in the working directory
if [[ ! -d "$MAPPED" ]]
then
  mkdir -p "$MAPPED"
fi

if [[ ! -d "$LOGS" ]]
then
  mkdir "$LOGS"
fi

# Create the isolates metadata used later for plotting
cat "/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_503/metadata.tsv" > "$WORKDIR"/metadata.tmp
echo "" >> "$WORKDIR"/metadata.tmp
tail -n +2 "/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_591/metadata.tsv" >> "$WORKDIR"/metadata.tmp
echo "" >> "$WORKDIR"/metadata.tmp
tail -n +2 "/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_643/metadata.tsv" >> "$WORKDIR"/metadata.tmp
echo "" >> "$WORKDIR"/metadata.tmp
tail -n +2 "/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_644/metadata.tsv" >> "$WORKDIR"/metadata.tmp
awk -F'\t' 'NR==1 {header=$0; next} {data[$3]=$0} END {print header; for (name in data) print data[name]}' "$WORKDIR"/metadata.tmp > "$WORKDIR"/metadata.tsv
rm "$WORKDIR"/*.tmp

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

# We create the raw reads folder and extract the reads from the bam files
if [[ ! -d "$RAW_READS" ]]
then
  mkdir "$RAW_READS"
fi

# Extract the reads from the bam files...
# From the hot pools
for j in "$BAMS"/h*
do
    sample=$(basename "${j}")
    echo "Extracting reads from ${sample}"
    samtools fastq \
        -F 4 \
        -1 "$RAW_READS"/"${sample}_1.fq.gz" \
        -2 "$RAW_READS"/"${sample}_2.fq.gz" \
        "${BAMS}/${sample}/s__Lactiplantibacillus_plantarum.bam"
done

# From the cold pools
for j in "$BAMS"/c*
do
    sample=$(basename "${j}")
    echo "Extracting reads from ${sample}"
    samtools fastq \
        -F 4 \
        -1 "$RAW_READS"/"${sample}_1.fq.gz" \
        -2 "$RAW_READS"/"${sample}_2.fq.gz" \
        "${BAMS}/${sample}/s__Lactiplantibacillus_plantarum.bam"
done

# Select the L. plantarum isolates
grep "s__Lactiplantibacillus_plantarum" "$NAME2TAXON" |\
    grep -v "MAG"|\
    cut -f1 |\
    cut -d "_" -f1 > "$WORKDIR"/isolates.tsv

while IFS=$'\t' read -r sample
do
    echo "Extracting reads from ${sample}"
    isolate="i${sample}"
    samtools fastq \
        -F 4 \
        -1 "$RAW_READS/${isolate}_1.fq.gz" \
        -2 "$RAW_READS/${isolate}_2.fq.gz" \
        "${BAMS}/${isolate}/s__Lactiplantibacillus_plantarum.bam"
done < "$WORKDIR"/isolates.tsv

### Competetive reads mapping to the representative microbiome and extraction of the unmapped reads

# Create an index for the reference combined genome
bowtie2-build --threads 16 "$GENOMES"/combined.fa "$GENOMES"/combined

# Competitive mapping mapping against each of the reads sets
for i in "$RAW_READS"/*_1.fq.gz
do
    # Declare the variable name
    # If the sample is an isolates genome (starts with "i"), keep only the first field, else keep the two first fields
    if [[ "$i" == i* ]]
    then
        name=$(echo "$i" | cut -d "_" -f1)
    else
        name=$(echo "$i" | cut -d "_" -f1,2)
    fi

    # Create the results folder
    mkdir -p "$MAPPED"/${name}

    echo "Mapping sample ${name}"

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
    rm "$MAPPED"/${name}/header.sam
done

# This chunk calculates the statistics that we are intersted in:
# Number of reads mapped to each genome
# Number of reads mapped UNIQUELY to each genome

echo -e "Calculating reads mapped to each L. plantarum strain"

# Remove any previous temporary file
rm -r "$WORKDIR"/*.col

# For each of the original genomes we create a temporary file to store the number of reads mapped to it, as well as the sample and genome names and sizes
echo -e "Sample" > "$WORKDIR"/sample_name.col
echo -e "Genome" > "$WORKDIR"/genome_name.col
echo -e "Size" > "$WORKDIR"/genome_size.col

for j in "$GENOMES"/*.fasta
do
    # Strain variable is the strain name
    strain=$(basename -a "${j}" | cut -d "." -f1)
    # Genome variable is the fasta file for each strain
    genome=$(basename -a "${j}")
    # Create columns for reads mapping to each strain
    echo -e "${strain}" > "$WORKDIR"/${strain}_reads.col
    echo -e "${strain}" > "$WORKDIR"/${strain}_uniq.col
    # Create a column with strain names
    echo -e "${strain}" >> "$WORKDIR"/genome_name.col
    # Create a column with strain genome size
    seqkit stats "$GENOMES"/"${genome}" | awk 'NR==2 {print $7}' >> "$WORKDIR"/genome_size.col
done

# For each sample...
for i in "$RAW_READS"/*_1.fq.gz
do
    # If the sample is an isolates genome (starts with "i"), keep only the first field, else keep the two first fields
    if [[ "$i" == i* ]]
    then
        name=$(echo "$i" | cut -d "_" -f1)
    else
        name=$(echo "$i" | cut -d "_" -f1,2)
    fi

    # Add the sample name to the file "sample_name.tmp"
    echo ${sample} >> "$WORKDIR"/sample_name.col

    # For each of the original genomes...
    for j in "$GENOMES"/*.fasta
    do
        # Strain variable is the strain name
        strain=$(basename -a "${j}" | cut -d "." -f1)
        # Extract the number of reads mapped to the genome and add it to ""$WORKDIR"/${j}_reads.tmp"
        samtools view -c -F 4 "$MAPPED"/${sample}/${strain}.bam >> "$WORKDIR"/${strain}_reads.col
        # Extract the number of reads mapped uniquely to the genome (with MAPQ>3) and add it to ""$WORKDIR"/${j}_uniq.tmp"
        samtools view -c -F 4 -q 4 "$MAPPED"/${sample}/${strain}.bam >> "$WORKDIR"/${strain}_uniq.col
    done 

done

paste "$WORKDIR"/genome_name.col "$WORKDIR"/genome_size.col > "$WORKDIR"/genome_size.tsv
paste "$WORKDIR"/sample_name.col "$WORKDIR"/*_reads.col > "$WORKDIR"/reads_mapped.tsv
paste "$WORKDIR"/sample_name.col "$WORKDIR"/*_uniq.col > "$WORKDIR"/uniq_mapped.tsv

#rm -r "$WORKDIR"/*.col

echo -e "Mapping finished"
