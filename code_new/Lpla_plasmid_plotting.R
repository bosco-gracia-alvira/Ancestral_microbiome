library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(rtracklayer)

# Paths to the files that we want to import (and export)
cov_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Lpla_plasmid/"
gff_path <- paste0(cov_path,"/reference/LplaWF.gff")
metadata_isolates_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/SNPs_analysis/metadata.tsv"
visuals_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/Lpla_plasmid"

# Import the coverage data and set the position in the beginning of each genomic element
cov <- fread(paste0(cov_path,"coverage.tsv"), header=T)
cov <- cov %>% 
        mutate(
          position = case_when(
            gene == "recA" ~ position-102018,
            gene == "plasmid" ~ position)) %>% 
          mutate(
            source = if_else(grepl("^.F", sample), "pool", "isolate"))

# Extract the coverage of recA
recA_mean <- cov %>%
  filter(gene == "recA") %>% 
  select(sample,coverage) %>% 
  group_by(sample) %>% 
  summarise(recA_mean = sum(coverage)/(1142)) # The value 1142 is the length of the recA gene

### Isolates analysis

# In this part of the analysis we want to check for the presence of the genomic island in each isolate
cov_isolates <- cov %>% filter(source == "isolate")

# Import the metadata of all the isolates and merge it with the coverage
metadata_isolates <- fread(metadata_isolates_path, header=T)
cov_merged <- merge(cov_isolates, metadata_isolates, by="sample")
cov_merged <- merge(cov_merged, recA_mean, by="sample")
cov_merged <- cov_merged %>% 
                mutate(norm_coverage = coverage/recA_mean)

# New genotype column, that matches the temperature, with three exceptions
cov_merged <- cov_merged %>% 
                mutate(Genotype = case_when(
                    sample %in% c("B25", "X76", "X371") ~ "S239",
                    Temperature == "HOT" ~ "S239",
                    Temperature == "BASE" ~ "B89",
                    Temperature == "COLD" ~ "S103",
                    Temperature == "COLD-CONSTANT" ~ "S103"))

# Subset the coverages from the plasmid and the housekeeping gene RecA
plasmid <- cov_merged %>% filter(gene=="plasmid")
recA <- cov_merged %>% filter(gene=="recA")

# Plot and save the data :)
# Coverage of the whole plasmid per isolate
plot_cov_plasmid <- ggplot(plasmid, aes(x = position, y = coverage, colour = Genotype)) +
                        annotate("rect", xmin = 102018, xmax = 184830, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightgreen") +
                        geom_point(size=0.03) +
                        #scale_y_continuous(trans='log10') +
                        scale_colour_manual(values = c("S239" = "red", "B89" = "grey", "S103" = "lightblue")) +
                        geom_hline(aes(yintercept = recA_mean), linetype = "dashed", color = "black") +
                        scale_x_continuous(labels = scales::scientific) +
                        facet_wrap(~ sample, scales = "free_y") +
                        labs(title = "Coverage of the plasmid in each the isolate", x = "Position", y = "Coverage") +
                        theme_minimal() +
                        theme(
                          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                          axis.text.y = element_text(size = 10),
                          axis.title.x = element_text(size = 14),
                          axis.title.y = element_text(size = 14))
plot_cov_plasmid

ggsave(filename = "plot_cov_plasmid.png", plot = plot_cov_plasmid, path = visuals_path, width = 24, height = 12, dpi = 300)

# Coverage of the whole plasmid for all the isolates
plot_cov_plasmid_all <- ggplot(plasmid, aes(x = position, y = coverage, colour = Genotype)) +
                            annotate("rect", xmin = 102018, xmax = 184830, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightgreen") +
                            geom_point(size=0.03) +
                            scale_colour_manual(values = c("S239" = "red", "B89" = "grey", "S103" = "lightblue")) +
                            scale_x_continuous(labels = scales::scientific) +
                            labs(title = "Overlapped coverages of the plasmid in all the isolates", x = "Position", y = "Coverage") +
                            theme_minimal() +
                            theme(
                              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                              axis.text.y = element_text(size = 10),
                              axis.title.x = element_text(size = 14),
                              axis.title.y = element_text(size = 14))

plot_cov_plasmid_all

ggsave(filename = "plot_cov_plasmid_all.png", plot = plot_cov_plasmid_all, path = visuals_path, width = 24, height = 12, dpi = 300)

# Coverage of recA gene
plot_cov_recA <- ggplot(recA, aes(x = position, y = coverage, colour = Genotype)) +
                    geom_point(size=0.03) +
                    #scale_y_continuous(trans='log10') +
                    scale_colour_manual(values = c("S239" = "red", "B89" = "grey", "S103" = "lightblue")) +
                    geom_hline(aes(yintercept = recA_mean), linetype = "dashed", color = "black") +
                    facet_wrap(~ sample) +
                    labs(title = "RecA coverage per isolate", x = "Position", y = "Coverage") +
                    theme_minimal() +
                    theme(
                    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                    axis.text.y = element_text(size = 10),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14)
                  )
plot_cov_recA

ggsave(filename = "plot_cov_recA.png", plot = plot_cov_recA, path = visuals_path, width = 24, height = 12, dpi = 300)

plot_cov_norm <- ggplot(plasmid, aes(x = position, y = norm_coverage, colour = Genotype)) +
                    annotate("rect", xmin = 102018, xmax = 184830, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightgreen") +
                    geom_point(size=0.03) +
                    scale_colour_manual(values = c("S239" = "red", "B89" = "grey", "S103" = "lightblue")) +
                    scale_x_continuous(labels = scales::scientific) +
                    facet_wrap(~ sample, scales = "free_y") +
                    labs(title = "Normalised coverage of the plasmid in each the isolate", x = "Position", y = "Normalised coverage") +
                    theme_minimal() +
                    theme(
                      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                      axis.text.y = element_text(size = 10),
                      axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14))
plot_cov_norm

ggsave(filename = "plot_cov_norm.png", plot = plot_cov_norm, path = visuals_path, width = 24, height = 12, dpi = 300)


# Import GFF of the reference genome
gff <- as.data.frame(import(gff_path))

# Extract the gene annotations from the plasmid (contig_2)
gene_annotations <- gff %>%
  filter(seqnames == "contig_2") %>% 
  select(seqnames, start, end, type, gene_id = ID, Name)

# Create a GRanges object for the plasmid's coverage data
coverage_gr <- GRanges(seqnames = plasmid$contig,
                       ranges = IRanges(start = plasmid$position, end = plasmid$position),
                       coverage = plasmid$coverage,
                       sample = plasmid$sample)

# Create a GRanges object for the gene annotations
genes_gr <- GRanges(seqnames = gene_annotations$seqnames,
                    ranges = IRanges(start = gene_annotations$start, end = gene_annotations$end),
                    gene_id = gene_annotations$gene_id,
                    Name = gene_annotations$Name)

# Find overlaps between coverage data and gene annotations
overlaps <- findOverlaps(coverage_gr, genes_gr)

# Extract the coverage values, gene IDs, sample IDs, and gene functions for the overlapping regions
coverage_values <- coverage_gr[queryHits(overlaps)]$coverage
gene_ids <- genes_gr[subjectHits(overlaps)]$gene_id
sample_ids <- coverage_gr[queryHits(overlaps)]$sample
gene_functions <- genes_gr[subjectHits(overlaps)]$Name


# Create a data frame with the coverage values, gene IDs, sample IDs, and gene functions
coverage_gene_df <- data.frame(gene_id = gene_ids, coverage = coverage_values, sample = sample_ids, gene_function = gene_functions)

# Calculate the average coverage for each gene for each sample and associate a genotype per sample
average_coverage <- coverage_gene_df %>%
  group_by(gene_id, sample, gene_function) %>%
  summarize(average_coverage = mean(coverage), .groups = 'drop') %>% 
  left_join(unique(cov_merged[,c(1,17)]), by = "sample") %>%
  mutate(transposable = if_else(grepl("Transposase|transposase|Mobile", gene_function), "red", "black"))

average_B25 <- average_coverage %>% filter(sample == "B25")

ggplot(average_B25, aes(x = gene_id, y = average_coverage, colour = gene_function)) +
  geom_point()


ggplot(average_coverage, aes(x = gene_id, y = average_coverage, colour = transposable)) +
  geom_point() +
  scale_colour_identity() +
  theme_minimal() +
  facet_wrap(~ sample, scales = "free_y") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = NULL)

  scale_colour_manual(values = c("S239" = "red", "B89" = "grey", "S103" = "lightblue")) +
  scale_x_continuous(labels = scales::scientific) +
  labs(title = "Normalised coverage of the plasmid in each the isolate", x = "Position", y = "Normalised coverage") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14))



### Pools analysis

# In this part we want to study the plasmid dynamics over the generations
cov_pools <- cov %>% filter(source == "pool")

cov_pools_meta <- cov_pools %>% 
                      mutate(
                        generation = as.numeric(str_extract(sample, "(?<=F)\\d+")),
                        replicate = str_extract(sample, "(?<=_)\\w+$"),
                        temperature = case_when(
                          str_detect(sample, "h") ~ "hot",
                          str_detect(sample, "c") ~ "cold"))

cov_summed_pools <- cov_pools_meta %>%
                        mutate(sample = str_extract(sample, "^.F\\d+")) %>% 
                        group_by(position,sample,gene,generation,temperature) %>% 
                        summarize(total_coverage = sum(coverage), .groups = 'drop')

summed_recA <- cov_summed_pools %>%
  filter(gene == "recA") %>% 
  select(sample,total_coverage) %>% 
  group_by(sample) %>% 
  summarise(recA_mean = sum(total_coverage)/(1142)) # The value 1142 is the length of the recA gene

cov_summed_pools_plasmid <- cov_summed_pools %>%  
                                    filter(gene=="plasmid") %>% 
                                    left_join(summed_recA, by="sample") %>% 
                                    mutate(recA_mean = replace_na(recA_mean, 0))

plot_cov_pools <- ggplot(cov_summed_pools_plasmid, aes(x = position, y = total_coverage, colour = temperature)) +
                      annotate("rect", xmin = 102018, xmax = 184830, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightgreen") +
                      geom_point(size=0.03) +
                      geom_hline(aes(yintercept = recA_mean), linetype = "dashed", color = "black") +
                      scale_colour_manual(values = c("hot" = "red", "cold" = "lightblue")) +
                      scale_x_continuous(labels = scales::scientific) +
                      facet_wrap(~ sample, scales = "free_y") +
                      labs(title = "Plasmid coverage per generation", x = "Position", y = "Normalised coverage") +
                      theme_minimal() +
                      theme(
                        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                        axis.text.y = element_text(size = 10),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14))
plot_cov_pools

ggsave(filename = "plot_cov_pools.png", plot = plot_cov_pools, path = visuals_path, width = 24, height = 12, dpi = 300)



