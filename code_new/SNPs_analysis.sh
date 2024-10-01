#!/bin/bash
# This script maps the reads to the species reference pangenome and calculates the SNPs
# It uses reads that were previously mapped competitively to the whole microbiome
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/SNPs_analysis"
GENOMES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome/genomes"
READS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome/mapped"
PLOT="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/code_new/"

# For each taxon, do the whole analysis...
for i in $(basename "$GENOMES/*.fasta")
do
    # Extract the taxon name
    taxon="${i%.fasta}"
    
    echo "
    Analysing taxon ${taxon}
    "

    # Set folders structure
    SNPS="$WORKDIR/${taxon}"
    RAW_READS="$SNPS/reads"
    LOGS="$SNPS/logs"
    BAMS="$SNPS/bams"

    ### COMMANDS
    IFS="
    "

    if [[ ! -f "$RAW_READS" ]]
    then
        mkdir -p "$RAW_READS"
    fi

    if [[ ! -f "$BAMS" ]]
    then
        mkdir -p "$BAMS"
    fi

    if [[ ! -f "$LOGS" ]]
    then
        mkdir -p "$LOGS"
    fi

    # Link the genome to the folder and index it
    ln -s "$GENOMES/${i}" "$SNPS"/ref.fa
    samtools faidx "$SNPS"/ref.fa

    # Extract the reads from the bam files...
    # From the hot pools
    for j in "$READS"/h*
    do
        sample=$(basename "${j}")
        echo "Extracting reads from ${sample}"
        samtools fastq \
            -F 4 \
            -1 "$RAW_READS"/"${sample}_1.fq.gz" \
            -2 "$RAW_READS"/"${sample}_2.fq.gz" \
            "${READS}/${sample}/${taxon}.bam"
    done

    # From the cold pools
    for j in "$READS"/c*
    do
        sample=$(basename "${j}")
        echo "Extracting reads from ${sample}"
        samtools fastq \
            -F 4 \
            -1 "$RAW_READS"/"${sample}_1.fq.gz" \
            -2 "$RAW_READS"/"${sample}_2.fq.gz" \
            "${READS}/${sample}/${taxon}.bam"
    done

    # From the isolates
    for j in $SAMPLES
    do
        sample="i${j}"
        echo "Extracting reads from ${sample}"
        samtools fastq \
                -F 4 \
                -1 "$RAW_READS/${j}_1.fq.gz" \
                -2 "$RAW_READS/${j}_2.fq.gz" \
                "${READS}/${sample}/${taxon}.bam"
    done

    # Map the reads to the reference pangenome

    # Create an index for the reference combined genome
    bowtie2-build --threads 16 "$SNPS"/ref.fa "$SNPS"/ref

    echo -e "sample\tcoverage" > "$SNPS"/coverage.txt

    # Competitive mapping against each of the reads sets
    for k in $(basename -a "$RAW_READS"/*_1.fq.gz)
    do
        sample="${k%_1.fq.gz}"
        r1="${k}"
        r2="${sample}_2.fq.gz"

        echo "Mapping ${sample}"
        # Map paired end reads using bowtie with stringent settings and output the result to a sam file
        bowtie2 \
            -x "$SNPS"/ref.fa \
            -q --very-sensitive \
            --no-mixed \
            --no-discordant \
            -1 "$RAW_READS/${r1}" \
            -2 "$RAW_READS/${r2}" \
            -S "$BAMS/${sample}.sam" \
            --threads 16 \
            --rg-id "${sample}" \
            --rg "SM:${sample}" > "$LOGS"/bowtie2_${sample}.log 2>&1

        # Turn the sam into bam and sort it
        samtools view \
            -bS \
            -@ 16 \
            "$BAMS/${sample}.sam" |\
        samtools sort \
            -@ 16 \
            -O bam \
            -o "$BAMS/${sample}_sorted.bam" \
            -

        # Delete the sam to save space
        rm "$BAMS/${sample}.sam"

        # Calculate the mean coverage of the sample (Truncated coverage 80%)
        mean_coverage=$(
            samtools depth "$BAMS/${sample}_sorted.bam" | \
            sort -k3,3nr | \
            awk '{ all[NR] = $3; sum+=$3 } END {if (NR==0) print 0; else { for (i=int(NR*0.1)+1; i<=int(NR*0.9); i++) s+=all[i]; print s/(int(NR*0.9)-int(NR*0.1)) } }')

        # Append the mean coverage to the coverage file
        echo -e "${sample}\t${mean_coverage}" >> "$SNPS/coverage.txt"
    done

    cd "$SNPS" || exit

    # List the samples with coverage above 10 and 5
    awk '$2 >= 10 {print "./bams/"$1"_sorted.bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_10.txt"
    awk '$2 >= 5 {print "./bams/"$1"_sorted.bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_5.txt"

    # This chunk counts the reference and alternative allele frequency in each position and in each sample (mpileup), then calls the SNPs (call) and filters the SNPs (no indels) with a quality above 20 and a depth above 5 
    bcftools mpileup -f "$SNPS"/ref.fa -b "$SNPS/coverage_5.txt" -Q 20 -D -d 50 -a DP,AD,QS,SCR -Ou --threads 16 | \
        bcftools call  --ploidy 1 -Ou -cv --threads 16 | \
        bcftools view -i 'QUAL > 20' -v snps -m2 -M2 -Ov - > "$SNPS/temp_${taxon}.vcf"
        # bcftools view -e 'FORMAT/DP[:0] < 3' -Ov - Filtering can be done in the R analysis step

    # We reformat the headers of the VCF file, that by default include the relative path to the bams
    bcftools view -h "$SNPS/temp_${taxon}.vcf" > "$SNPS/headers.txt"
    sed -i '' 's|./bams/||g; s|_sorted.bam||g' "$SNPS/headers.txt"
    bcftools reheader -h "$SNPS/headers.txt" -o "$SNPS/${taxon}.vcf" "$SNPS/temp_${taxon}.vcf"

    # We extract the frequency of the SNPs in each sample
    bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$SNPS/${taxon}.vcf" > "$SNPS/${taxon}.freq"

    # Delete the temporary file
    rm "$SNPS/temp_${taxon}.vcf"

    # This script plots the PCA of the samples based on the SNPs frequency
    Rscript "$WORKDIR"/../code/SNPs_plotting.Rmd "${taxon}"
done