#!/bin/bash
# This script maps the reads to the species reference pangenome and calculates the SNPs
# It uses reads that were previously mapped competitively to the whole microbiome
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/SNPs_analysis"
NAME2TAXON="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Graph_pangenome/name2taxon.tsv"
GENOMES="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome/genomes"
VISUALS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/SNPs_analysis"
BAMS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome/mapped"
CODE="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/code_new"

### COMMANDS
IFS=$'\n'

# Create the working directory
if [[ ! -f "$WORKDIR" ]]
then
    mkdir -p "$WORKDIR"
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

# Create the folder where plots will go
if [[ ! -f "$VISUALS" ]]
then
    mkdir -p "$VISUALS"
fi

# For each taxon, do the whole analysis...
for i in "$GENOMES"/*.fasta
do
    # Extract the taxon name
    taxon=$(basename "${i%.fasta}")

    echo "
    Analysing taxon ${taxon}
    "

    # Set folders structure
    SNPS="$WORKDIR/${taxon}"
    if [[ ! -f "$SNPS"/bams ]]
    then
        mkdir -p "$SNPS"/bams
    fi

    echo -e "sample\tcoverage" > "$SNPS"/coverage.txt

    # Link the pools to the bams folder
    for k in "$BAMS"/{c,h}*
    do
        # Set the sample
        sample=$(basename -a "${k}")

        # Link the bams to the SNPs analysis folder
        ln -s "$BAMS/${sample}/${taxon}.bam" "$SNPS/bams/${sample}.bam"

        # Inform the user of the current state
        echo "Calculating coverage for sample ${sample} in taxon ${taxon}"
        
        # Calculate the mean coverage of the sample (Truncated coverage 80%)
        mean_coverage=$(
            samtools depth "$SNPS/bams/${sample}.bam" | \
            sort -k3,3nr | \
            awk '{ all[NR] = $3; sum+=$3 } END {if (NR==0) print 0; else { for (i=int(NR*0.1)+1; i<=int(NR*0.9); i++) s+=all[i]; print s/(int(NR*0.9)-int(NR*0.1)) } }')

        # Append the mean coverage to the coverage file
        echo -e "${sample}\t${mean_coverage}" >> "$SNPS/coverage.txt"
    done

    # List the isolates that belong to the iterated taxon
    grep "${taxon}" "$NAME2TAXON" |\
    grep -v "MAG"|\
    cut -f1 |\
    cut -d "_" -f 1 > "$SNPS/isolates.tsv"

    # Link only the isolates belonging to the correct taxonomy to the bams folder
    while IFS=$'\t' read -r isolate
    do
        # Set the isolate
        sample="i${isolate}"

        # Link the bams to the SNPs analysis folder
        ln -s "$BAMS/${sample}/${taxon}.bam" "$SNPS/bams/${sample}.bam"

        # Inform the user of the current state
        echo "Calculating coverage for sample ${sample} in taxon ${taxon}"
        
        # Calculate the mean coverage of the sample (Truncated coverage 80%)
        mean_coverage=$(
            samtools depth "$SNPS/bams/${sample}.bam" | \
            sort -k3,3nr | \
            awk '{ all[NR] = $3; sum+=$3 } END {if (NR==0) print 0; else { for (i=int(NR*0.1)+1; i<=int(NR*0.9); i++) s+=all[i]; print s/(int(NR*0.9)-int(NR*0.1)) } }')

        # Append the mean coverage to the coverage file
        echo -e "${sample}\t${mean_coverage}" >> "$SNPS/coverage.txt"
    done < "$SNPS/isolates.tsv"

    # Change to the local directory because the vcf uses local paths from coverage_5.txt
    cd "$SNPS" || exit

    # List the samples with coverage above 10 and 5
    awk '$2 >= 10 {print "./bams/"$1".bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_10.txt"
    awk '$2 >= 5 {print "./bams/"$1".bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_5.txt"

    # This chunk counts the reference and alternative allele frequency in each position and in each sample (mpileup), then calls the SNPs (call) and filters the SNPs (no indels) with a quality above 20 and a depth above 5 
    bcftools mpileup -f "$GENOMES/${taxon}.fasta" -b "$SNPS/coverage_5.txt" -Q 20 -D -d 50 -a DP,AD,QS,SCR -Ou --threads 16 | \
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
    #Rscript "$CODE"/SNPs_plotting.Rmd "${taxon}"
done
