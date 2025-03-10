#!/bin/bash
# This script maps competitively the reads from each sample against the L_plantarum genomes S103 and S239.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
METADATA="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Lactiplantibacillus_plantarum/Pangenomic_Florida/metadata.tsv"
LOCATION_COLD="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"
LOCATION_ISOLATES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"

WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Lpla_dynamics_replicates"
RAW_READS="$WORKDIR/reads_val"
GENOMES="$WORKDIR/genomes"
LOGS="$WORKDIR/logs_val"
MAPPED="$WORKDIR/mapped_val"

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

# We create the graph pangenome for each strain using SuperPang in Vetlinux5, if the reference doesn't exist
if [[ ! -f "$GENOMES"/combined.fa ]]
then

    for strain in "S103" "S239" "B89"
    do

    # Create the folder for each strain in vetlinux
    ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

        if [[ ! -d ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain"/Genomes ]]
        then  
            mkdir -p ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain"/Genomes
        fi

FOO

    # Grep the genomes that belong to each strain and copy them to the strain folder
    grep "$strain" "$METADATA" |\
    while IFS=$'\t' read -r Strain Fasta Pool Temperature Generation Population Replicate Genotype rest
    do
        if [[ ! "$Strain" =~ X595_|X79_|X81_ ]]
        then
            rsync -av \
                        "/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Isolates/$Fasta" \
                        vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain"/Genomes/$Fasta
        fi
    done

    # Construct the graph pangenome in vetlinux
    ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

        eval "\$(conda shell.bash hook)"
        conda activate SuperPang-0.9

        cd ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain"

        echo "
            Making pangenome for $strain.
            "

        SuperPang.py \
                --fasta ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain"/Genomes/*.fasta \
                --assume-complete \
                --output-dir ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain" \
                --force-overwrite \
                --threads 24

        cat ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain"/assembly.fasta |\
            seqkit replace -p .+ -r "${strain}_{nr}" --nr-width 3 \
            > ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain"/"$strain".fasta

        rm -r ~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain"/Genomes

FOO

    # Copy back the results to the genomes folder
    rsync -av \
        vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Ancestral_microbiome/Graph_pangenome/"$strain"/"$strain".fasta \
        "$GENOMES"

    done

    # Combine the three genomes into a single reference
    cat "$GENOMES"/S103.fasta "$GENOMES"/S239.fasta "$GENOMES"/B89.fasta > "$GENOMES"/combined.fa
fi

# We create the raw reads folder and simulate reads sets in it
if [[ ! -d "$RAW_READS" ]]
then
  mkdir "$RAW_READS"
fi

# We make three mock communities with different proportions of strains
# Set the seed
seed=$(date +%s)

# Simulate 3 reads sets with 200000 reads from each reference pangenome
for strain in S103 S239 B89
do
  iss generate \
        -g "$GENOMES/${strain}.fasta" \
        -n 200000 \
        --model hiseq \
        --seed $seed \
        --cpus 16 \
        --compress \
        --quiet \
        --output "$RAW_READS/i${strain}"

  for i in 1 2
  do
    mv "$RAW_READS/i${strain}"_R${i}.fastq.gz "$RAW_READS/i${strain}"_R${i}.fq.gz 
  done
done

# Remove the temp files
rm "$RAW_READS"/*.vcf
rm "$RAW_READS"/*abundance.txt

# Now we build the different reads mixes

# Mix217
for i in 1 2
do
  # Subsample i103 to the proportion 0.2
  gzcat "$RAW_READS/iS103_R${i}.fq.gz" |\
      seqkit sample \
        -s $seed \
        -p 0.2 \
        -o "$RAW_READS/iS103_${i}_temp.fq.gz"

  # Subsample i239_1 to the proportion 0.1
  gzcat "$RAW_READS/iS239_R${i}.fq.gz" |\
      seqkit sample \
        -s $seed \
        -p 0.1 \
        -o "$RAW_READS/iS239_${i}_temp.fq.gz"

  # Subsample iB89_1 to the proportion 0.7
  gzcat "$RAW_READS/iB89_R${i}.fq.gz" |\
      seqkit sample \
        -s $seed \
        -p 0.7 \
        -o "$RAW_READS/iB89_${i}_temp.fq.gz"
  
  # Concatenate the forward reads from both strains and shuffle them, prior to further subsampling
  gzcat "$RAW_READS/iS103_${i}_temp.fq.gz" "$RAW_READS/iS239_${i}_temp.fq.gz" "$RAW_READS/iB89_${i}_temp.fq.gz" | seqkit shuffle -s $seed -o "$RAW_READS/mix217_${i}.fq.gz"

done  

# Mix333
for i in 1 2
do
  # Subsample i103 to the proportion 0.2
  gzcat "$RAW_READS/iS103_R${i}.fq.gz" |\
      seqkit sample \
        -s $seed \
        -p 0.33 \
        -o "$RAW_READS/iS103_${i}_temp.fq.gz"

  # Subsample i239_1 to the proportion 0.1
  gzcat "$RAW_READS/iS239_R${i}.fq.gz" |\
      seqkit sample \
        -s $seed \
        -p 0.33 \
        -o "$RAW_READS/iS239_${i}_temp.fq.gz"

  # Subsample iB89_1 to the proportion 0.7
  gzcat "$RAW_READS/iB89_R${i}.fq.gz" |\
      seqkit sample \
        -s $seed \
        -p 0.33 \
        -o "$RAW_READS/iB89_${i}_temp.fq.gz"
  
  # Concatenate the forward reads from both strains and shuffle them, prior to further subsampling
  gzcat "$RAW_READS/iS103_${i}_temp.fq.gz" "$RAW_READS/iS239_${i}_temp.fq.gz" "$RAW_READS/iB89_${i}_temp.fq.gz" | seqkit shuffle -s $seed -o "$RAW_READS/mix333_${i}.fq.gz"

done

# Mix811
for i in 1 2
do
  # Subsample i103 to the proportion 0.2
  gzcat "$RAW_READS/iS103_R${i}.fq.gz" |\
      seqkit sample \
        -s $seed \
        -p 0.8 \
        -o "$RAW_READS/iS103_${i}_temp.fq.gz"

  # Subsample i239_1 to the proportion 0.1
  gzcat "$RAW_READS/iS239_R${i}.fq.gz" |\
      seqkit sample \
        -s $seed \
        -p 0.1 \
        -o "$RAW_READS/iS239_${i}_temp.fq.gz"

  # Subsample iB89_1 to the proportion 0.7
  gzcat "$RAW_READS/iB89_R${i}.fq.gz" |\
      seqkit sample \
        -s $seed \
        -p 0.1 \
        -o "$RAW_READS/iB89_${i}_temp.fq.gz"
  
  # Concatenate the forward reads from both strains and shuffle them, prior to further subsampling
  gzcat "$RAW_READS/iS103_${i}_temp.fq.gz" "$RAW_READS/iS239_${i}_temp.fq.gz" "$RAW_READS/iB89_${i}_temp.fq.gz" | seqkit shuffle -s $seed -o "$RAW_READS/mix811_${i}.fq.gz"

done

rm "$RAW_READS"/*_temp.*

echo -e "sample\tproportion\tseed" > "$LOGS"/validation.tsv

# Subsample the reads by decreasing proportions
for i in $(basename -a "$RAW_READS"/*_1.fq.gz | cut -d "_" -f1)
do
  for j in 10 20 40 60 80 100 200 300 400 500 600 700 800 900 1000 2000 5000 10000 150000
  do
    echo -e "${i}\t${j}\t$seed" >> "$LOGS"/validation.tsv

    gzcat "$RAW_READS"/${i}_1.fq.gz |\
    seqkit sample \
      -s $seed \
      -n ${j} \
      -o "$RAW_READS"/${i}-${j}_1.fastq.gz
    
    gzcat "$RAW_READS"/${i}_2.fq.gz |\
    seqkit sample \
      -s $seed \
      -n ${j} \
      -o "$RAW_READS"/${i}-${j}_2.fastq.gz
  done
done

# Create an index for the reference combined genome
bowtie2-build --threads 16 "$GENOMES"/combined.fa "$GENOMES"/combined

# Competitive mapping mapping against each of the reads sets
for i in $(basename -a "$RAW_READS"/*_1.fastq.gz | cut -d "_" -f1)
do
  # Create the results folder
  mkdir -p "$MAPPED"/${i}

  # Map paired end reads using bowtie with stringent settings and output the result to a sam file
  bowtie2 \
    -x "$GENOMES"/combined \
    -q \
    -D 500 -R 40 -N 0 -L 20 -i S,1,0.50 \
    --no-mixed \
    --no-discordant \
    -1 "$RAW_READS"/${i}_1.fastq.gz \
    -2 "$RAW_READS"/${i}_2.fastq.gz \
    -S "$MAPPED"/${i}/combined.sam \
    --threads 20 > "$LOGS"/validation_${i}.log 2>&1

  # Turn the sam into bam to save memory
  samtools view \
    -bS \
    -@ 16 \
    "$MAPPED"/${i}/combined.sam > "$MAPPED"/${i}/combined.bam

  # Delete the sam
  rm "$MAPPED"/${i}/combined.sam

  # Sort the bam
  samtools sort \
    -@ 16 \
    -O bam \
    -o "$MAPPED"/${i}/combined_sorted.bam \
    "$MAPPED"/${i}/combined.bam

  # Delete the unsorted bam
  rm "$MAPPED"/${i}/combined.bam

  # Extract the header of the bam
  samtools view -H "$MAPPED"/${i}/combined_sorted.bam > "$MAPPED"/${i}/header.sam

  # Index the combined bam
  samtools index "$MAPPED"/${i}/combined_sorted.bam

  # For each of the original genomes...
  for j in $(basename -a "$GENOMES"/*.fasta | cut -d "." -f1)
  do  
    # Make a list with all the contigs belonging to that genome
    contigs=$(samtools idxstats "$MAPPED"/${i}/combined_sorted.bam | cut -f 1 | grep "^${j}")

    # From the combined bam, extract the reads that map the contigs of the specific genome "j"
    # Then, add the header and create a new file with the reads mapping to "j"
    samtools view -@ 16 -b "$MAPPED"/${i}/combined_sorted.bam $contigs |\
      samtools reheader "$MAPPED"/${i}/header.sam - > "$MAPPED"/${i}/${j}.bam
  done

  # Remove the combined bam file, the index and the header
  rm "$MAPPED"/${i}/combined_sorted.bam
  rm "$MAPPED"/${i}/combined_sorted.bam.bai
  rm "$MAPPED"/${i}/header.sam
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

echo -e "Calculating reads mapped to each L. plantarum strain"

# Remove any previous temporary file
rm -r "$WORKDIR"/*.col

# For each of the original genomes we create a temporary file to store the number of reads mapped to it, as well as the sample and genome names and sizes
echo -e "Sample" > "$WORKDIR"/sample_name_val.col
echo -e "Genome" > "$WORKDIR"/genome_name_val.col
echo -e "Size" > "$WORKDIR"/genome_size_val.col

for j in "$GENOMES"/*.fasta
do

    # Strain variable is the strain name
    strain=$(basename -a "${j}" | cut -d "." -f1)
    # Genome variable is the fasta file for each strain
    genome=$(basename -a "${j}")
    # Create columns for reads mapping to each strain
    echo -e "${strain}" > "$WORKDIR"/${strain}_reads_val.col
    echo -e "${strain}" > "$WORKDIR"/${strain}_uniq_val.col
    # Create a column with strain names
    echo -e "${strain}" >> "$WORKDIR"/genome_name_val.col
    # Create a column with strain genome size
    seqkit stats "$GENOMES"/"${genome}" | awk 'NR==2 {print $7}' >> "$WORKDIR"/genome_size_val.col
    
done

# Set the number of samples to process
numsamples=$(basename -a "$MAPPED"/* | wc -l)
processed=1

# For each sample...
for i in $(basename -a "$MAPPED"/*)
do
    
    # Define the variables
    sample=${i}
    echo -e "Extracting reads from sample ${sample} (${processed}/${numsamples})"
    processed=$((processed+1))

    # Add the sample name to the file "sample_name.tmp"
    echo ${sample} >> "$WORKDIR"/sample_name_val.col

    # For each of the original genomes...
    for j in "$GENOMES"/*.fasta
    do
        # Strain variable is the strain name
        strain=$(basename -a "${j}" | cut -d "." -f1)
        # Filter the bam file based on the mappability bed files
        bedtools intersect -v -abam "$MAPPED"/${sample}/${strain}.bam -b "$GENOMES"/genmap/S103.bed > "$MAPPED"/${sample}/${strain}_filt.bam
        bedtools intersect -v -abam "$MAPPED"/${sample}/${strain}.bam -b "$GENOMES"/genmap/S239.bed > "$MAPPED"/${sample}/${strain}_filt.bam
        bedtools intersect -v -abam "$MAPPED"/${sample}/${strain}.bam -b "$GENOMES"/genmap/B89.bed > "$MAPPED"/${sample}/${strain}_filt.bam
        # Extract the number of reads mapped to the genome and add it to ""$WORKDIR"/${j}_reads.tmp"
        samtools view -@ 16 -c -F 4 "$MAPPED"/${sample}/${strain}_filt.bam >> "$WORKDIR"/${strain}_reads_val.col
        # Extract the number of reads mapped uniquely to the genome (with MAPQ>3) and add it to ""$WORKDIR"/${j}_uniq.tmp"
        samtools view -@ 16 -c -F 4 -q 20 "$MAPPED"/${sample}/${strain}_filt.bam >> "$WORKDIR"/${strain}_uniq_val.col
    done

done

paste "$WORKDIR"/genome_name_val.col "$WORKDIR"/genome_size_val.col > "$WORKDIR"/genome_size_val.tsv
paste "$WORKDIR"/sample_name_val.col "$WORKDIR"/*_reads_val.col > "$WORKDIR"/reads_mapped_val.tsv
paste "$WORKDIR"/sample_name_val.col "$WORKDIR"/*_uniq_val.col > "$WORKDIR"/uniq_mapped_val.tsv

rm -r "$WORKDIR"/*.col

echo -e "Mapping finished"
