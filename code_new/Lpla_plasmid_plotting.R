library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)

# Paths to the files that we want to import (and export)
cov_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Lpla_plasmid/"
metadata_isolates_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/SNPs_analysis/metadata.tsv"
visuals_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/Lpla_plasmid"

# Import the coverage data and set the position in the beginning of each genomic element
cov <- fread(paste0(cov_path,"coverage.tsv"), header=T)
cov <- cov %>% 
        mutate(
          position = case_when(
            gene == "island" ~ position-102018,
            gene == "recA" ~ position-102018,
            gene == "plasmid" ~ position,
          )
        )

# Extract the coverage of recA
recA_mean <- cov %>%
  filter(gene == "recA") %>% 
  select(isolate,coverage) %>% 
  group_by(isolate) %>% 
  summarise(recA_mean = mean(coverage))

# Import the metadata of all the isolates and merge it with the coverage
metadata_isolates <- fread(metadata_isolates_path, header=T)
cov_merged <- merge(cov, metadata_isolates, by.x="isolate", by.y="sample")
cov_merged <- merge(cov_merged, recA_mean, by="isolate")
cov_merged <- cov_merged %>% 
                mutate(norm_coverage = coverage/recA_mean)

# New genotype column, that matches the temperature, with three exceptions
cov_merged <- cov_merged %>% 
                mutate(Genotype = case_when(
                    isolate %in% c("B25", "X76", "X371") ~ "S239",
                    Temperature == "HOT" ~ "S239",
                    Temperature == "BASE" ~ "B89",
                    Temperature == "COLD" ~ "S103",
                    Temperature == "COLD-CONSTANT" ~ "S103"))

# Subset the coverages from the plasmid and the housekeeping gene RecA
plasmid <- cov_merged %>% filter(gene=="plasmid")
island <- cov_merged %>% filter(gene=="island")
recA <- cov_merged %>% filter(gene=="recA")

# Plot and save the data :)
# Coverage of the whole plasmid per isolate
plot_cov_plasmid <- ggplot(plasmid, aes(x = position, y = coverage, colour = Genotype)) +
                        annotate("rect", xmin = 102018, xmax = 184830, ymin = 1, ymax = Inf, alpha = 0.2, fill = "lightgreen") +
                        geom_point(size=0.03) +
                        #scale_y_continuous(trans='log10') +
                        scale_colour_manual(values = c("S239" = "red", "B89" = "grey", "S103" = "lightblue")) +
                        geom_hline(aes(yintercept = recA_mean), linetype = "dashed", color = "black") +
                        scale_x_continuous(labels = scales::scientific) +
                        facet_wrap(~ isolate, scales = "free_y") +
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
  annotate("rect", xmin = 102018, xmax = 184830, ymin = 1, ymax = Inf, alpha = 0.2, fill = "lightgreen") +
  geom_point(size=0.03) +
  #scale_y_continuous(trans='log10') +
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
                    facet_wrap(~ isolate) +
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
                    annotate("rect", xmin = 102018, xmax = 184830, ymin = 1, ymax = Inf, alpha = 0.2, fill = "lightgreen") +
                    geom_point(size=0.03) +
                    scale_colour_manual(values = c("S239" = "red", "B89" = "grey", "S103" = "lightblue")) +
                    scale_x_continuous(labels = scales::scientific) +
                    facet_wrap(~ isolate, scales = "free_y") +
                    labs(title = "Normalised coverage of the plasmid in each the isolate", x = "Position", y = "Normalised coverage") +
                    theme_minimal() +
                    theme(
                      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                      axis.text.y = element_text(size = 10),
                      axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14))
plot_cov_norm

ggsave(filename = "plot_cov_norm.png", plot = plot_cov_norm, path = visuals_path, width = 24, height = 12, dpi = 300)




