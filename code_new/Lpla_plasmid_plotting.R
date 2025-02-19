library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(gggenomes)
library(gtools)

# Paths to the files that we want to import (and export)
cov_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Lpla_plasmid"
gff_path <- normalizePath(paste0(cov_path,"/reference/LplaWF.gff"))
fa_path <- normalizePath(paste0(cov_path,"/reference/LplaWF.fa"))
metadata_isolates_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/SNPs_analysis/metadata.tsv"
visuals_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/Lpla_plasmid"

####################################################################
# DRAW THE GENETIC MAP OF THE PLASMID AND THE COVERAGE PER ISOLATE #
####################################################################

# Read the plasmid fasta file in order to get the length
fa <- read_seq_len(fa_path) %>% 
  filter(seq_id=="contig_2")

# Read the gff file of the plasmid and modify it to fit gggenomes
gff <- read_gff3(gff_path) %>%
  dplyr::filter(seq_id == "contig_2") %>% 
  dplyr::select(seq_id, start, end, strand, type, name) %>%
  dplyr::filter(!is.na(name)) %>%
  dplyr::mutate(colour = case_when(
    grepl("Transposase", name, ignore.case = TRUE) ~ "purple",
    grepl("Mobile", name, ignore.case = TRUE) ~ "yellow",
    grepl("tRNA-", name) ~ "brown",
    TRUE ~ "grey"
  ))

# Create a new feature with the genetic island range
island <- data.frame(
  seq_id = "contig_2",
  start = 102018,
  end = 184830,
  strand = "+",
  type = "Island",
  Name = "Island",
  colour = "green")

# Import the isolates' metadata and coverage data and set the position in the beginning of each genomic element
metadata_isolates <- fread(metadata_isolates_path, header=T) %>% 
  select(c(sample,Temperature,Generation,Population,Replicate))

cov_feats <- fread(paste0(cov_path,"/coverage.tsv"), header=T) %>% 
  filter(contig == "contig_2") %>%
  filter(coverage > 5) %>% 
  left_join(metadata_isolates, by=join_by(sample)) %>% 
  filter(!is.na(Temperature)) %>% 
  mutate(start = position,
         end = position,
         seq_id = contig,
         identity = 1) %>%
  select(-c(contig,position)) %>% 
  mutate(colour = case_when(
    sample %in% c("B25", "X76", "X371") ~ "red",
    Temperature == "HOT" ~ "red",
    Temperature == "BASE" ~ "grey",
    Temperature == "COLD" ~ "lightblue",
    Temperature == "COLD-CONSTANT" ~ "lightblue")) %>% 
  arrange(colour, levels(c("red", "lightblue", "grey")))

# Create y_pos, which will map each sample to the y axis of the plot
sample_mapping <- setNames(seq_along(unique(cov_feats$sample)), unique(cov_feats$sample))

cov_feats <- cov_feats %>%
  mutate(y_pos = sample_mapping[sample])

# Create a gggenomes object that contains all the genetic features
p1 <- gggenomes(seqs= fa,feats = gff) %>% 
  add_feats(island) %>% 
  add_feats(cov_feats)

# Add the genetic map and the island range on the background
p2 <- p1 +
  geom_feat(data = feats(island), aes(colour=colour), alpha = 0.2, position="pile", linewidth = 9) +
  scale_colour_identity() +
  labs(title = "Genetic map of the plasmid with Genetic Island", x = "Position") +
  geom_seq() +
  geom_gene(aes(fill = colour), shape = c(7, 4), intron_types = c("CDS","tRNA"), size = 5) +
  scale_fill_identity() +
  labs(title = "Genetic map of the plasmid", x = "Position") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20),
  )

p2

# Add the coverage data from all the sequences
p3 <- p2 + 
  geom_segment(data = feats(cov_feats), aes(x = start, xend = start, y = 0.02*(y_pos), yend = 0.02*(y_pos + identity), colour = colour), linewidth = 1) +
  geom_hline(data = feats(cov_feats), aes(yintercept = 0.02*(y_pos+1), colour = "black"), linewidth = 0.05) +
  geom_text(data = feats(cov_feats), aes(x = -2000, y = 0.02*(y_pos+0.5), label = sample, colour = colour), hjust = 1, vjust = 0.5, size = 2)
  
p3

ggsave(filename = "map_isolates_all.png", plot = p3, path = visuals_path, width = 20, height = 12, dpi = 300)


# We subset one representative per plasmid class
p4 <- p1
p4$data$feats[[3]] <- p4$data$feats[[3]] %>% 
          filter(sample %in% c("B89", "B246", "S89", "X596")) 

p5 <- p4 + 
  geom_feat(data = feats(island), aes(colour=colour), alpha = 0.2, position="pile", linewidth = 9) +
  scale_colour_identity() +
  labs(title = "Genetic map of the plasmid with Genetic Island", x = "Position") +
  geom_seq() +
  geom_gene(aes(fill = colour), shape = c(7, 4), intron_types = c("CDS","tRNA"), size = 5) +
  scale_fill_identity() +
  labs(title = "Genetic map of the plasmid", x = "Position") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20)) +
  geom_segment(data = feats(cov_feats), aes(x = start, xend = start, y = 0.7, yend = 0.8 * identity, colour = colour), linewidth = 1) +
  geom_coverage(data = feats(cov_feats), aes(z = coverage, fill = colour), height = 0.8, offset = 0.3) +
  scale_fill_identity() +
  facet_wrap(~sample)

p5

ggsave(filename = "map_isolates_rep.png", plot = p5, path = visuals_path, width = 20, height = 12, dpi = 300)


### Pools analysis

# In this part we want to study the plasmid dynamics over the generations

# Create a feature with the alignment to the pools
cov_pools <- fread(paste0(cov_path,"/coverage.tsv"), header=T) %>% 
                filter(contig == "contig_2") %>%
                filter(grepl("^.F", sample)) %>%
                mutate(sample = str_extract(sample, "^.F\\d+")) %>% 
                group_by(position,sample,contig) %>% 
                summarize(total_coverage = sum(coverage), .groups = 'drop') %>% 
                #filter(total_coverage > 5) %>%
                mutate(start = position,
                       end = position,
                       seq_id = contig,
                       identity = 1) %>%
                select(-contig) %>% 
                mutate(
                  Generation = as.numeric(str_extract(sample, "(?<=F)\\d+")),
                  Temperature = case_when(
                    str_detect(sample, "h") ~ "hot",
                    str_detect(sample, "c") ~ "cold"),
                  colour = case_when(
                    str_detect(sample, "F0") ~ "grey",
                    str_detect(sample, "h") ~ "red",
                    str_detect(sample, "c") ~ "lightblue")) %>%
                mutate(sample = factor(sample , levels = mixedsort(unique(sample)))) %>%
                arrange(sample)

# Create y_pos, which will map each sample to the y axis of the plot
sample_mapping <- setNames(seq_along(unique(cov_pools$sample)), unique(cov_pools$sample))

cov_pools <- cov_pools %>%
  mutate(y_pos = sample_mapping[sample])


t1 <- gggenomes(seqs= fa,feats = gff) %>% 
                          add_feats(island) %>% 
                          add_feats(cov_pools)

# Add the genetic map and the island range on the background
t2 <- t1 +
  geom_feat(data = feats(island), aes(colour=colour), alpha = 0.2, position="pile", linewidth = 9) +
  scale_colour_identity() +
  labs(title = "Genetic map of the plasmid with Genetic Island", x = "Position") +
  geom_seq() +
  geom_gene(aes(fill = colour), shape = c(7, 4), intron_types = c("CDS","tRNA"), size = 5) +
  scale_fill_identity() +
  labs(title = "Genetic map of the plasmid", x = "Position") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20))

t2

# Add the coverage data from all the sequences
t3 <- t2 + 
  geom_segment(data = feats(cov_pools), aes(x = start, xend = start, y = 0.02*(y_pos), yend = 0.02*(y_pos + identity), colour = colour), linewidth = 1) +
  geom_hline(data = feats(cov_pools), aes(yintercept = 0.02*(y_pos+1), colour = "black"), linewidth = 0.05) +
  geom_text(data = feats(cov_pools), aes(x = -2000, y = 0.02*(y_pos+0.5), label = sample, colour = colour), hjust = 1, vjust = 0.5, size = 2)

t3

ggsave(filename = "map_pools.png", plot = t3, path = visuals_path, width = 20, height = 12, dpi = 300)
















cov_pools <- fread(paste0(cov_path,"/coverage.tsv"), header=T) %>%
                                              filter(source == "pool") %>% 
                                              mutate(
                                                generation = as.numeric(str_extract(sample, "(?<=F)\\d+")),
                                                replicate = str_extract(sample, "(?<=_)\\w+$"),
                                                temperature = case_when(
                                                  str_detect(sample, "h") ~ "hot",
                                                  str_detect(sample, "c") ~ "cold"))

cov_pools_meta <- cov_pools %>% 
  

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


















# Import the coverage data and set the position in the beginning of each genomic element
cov <- fread(paste0(cov_path,"/coverage.tsv"), header=T)
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



