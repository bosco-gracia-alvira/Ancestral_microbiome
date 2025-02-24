
# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# The first argument will be the value of $REPLY
reply <- args[1]
#reply <- "s__Lactiplantibacillus_plantarum"


# Print the value of reply (for debugging purposes)
print(paste("The following taxon will be used:", reply))

knitr::opts_knit$set(root.dir = paste0("/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/",reply,"/SNPs_analysis"))

###r Loading packages
library(lattice)
library(data.table)
library(ggplot2)
library(tidyverse)
library(stringr)
library(dplyr)
library(VariantAnnotation)

#######################################

###{r Loading the data}

# Paths to the files that we want to import
data_path <- paste0("/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/SNPs_analysis/",reply)
metadata_isolates_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/SNPs_analysis/metadata.tsv"
visuals_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/visuals/SNPs_analysis"

# freq is the table that comes from the vcf query that contains allele frequencies for each sample
freq <- fread(paste0(data_path,"/",reply,".freq"),sep="\t",header=T, colClasses = "character")

# Make the headers pretty :)
colnames(freq) <- gsub("^#", "", colnames(freq))          # Remove "#"
colnames(freq) <- gsub("^.*\\[.*?\\]", "", colnames(freq))  # Remove anything between brackets
colnames(freq) <- gsub(":AD$", "", colnames(freq))        # Remove ":AD"

# Remove four isolates that I know that are contaminated and have intermediate allele frequencies
if (reply == "s__Lactiplantibacillus_plantarum") {
  freq <- freq %>%
    dplyr::select(-iX79, -iX81, -iX595, -iX371)
}

#######################################


###{r Process the frequency table}

# Melt the dataframe to long format
freq_long <- data.table::melt(freq, id.vars = colnames(freq)[1:4], variable.name = "Sample", value.name = "AF")

# Split the AF column into Reference and alternative using "," as splitting element
freq_long <- freq_long[, c("RO", "AO") := tstrsplit(AF, ",")]

# Turn the new columns numeric
freq_long <- freq_long[, `:=`(RO = as.numeric(RO), AO = as.numeric(AO))]

# Calculate the summed allele depth (SAD) and relative frequencies from absolute number of counts
freq_long <- freq_long[, "SAD" := (RO + AO)]
freq_long <- freq_long[, `:=`(RO = RO / (RO + AO), AO = AO / (RO + AO))]

# Create a new column with unique position information
freq_long <- freq_long %>%
  mutate(CHROM_POS = paste(CHROM, POS, sep = "_"))

# Make the table wide again
freq_wide <- freq_long %>%
  dplyr::select(Sample, CHROM_POS, RO) %>%  # Select relevant columns
  pivot_wider(names_from = CHROM_POS, values_from = RO) %>%  # Pivot to wide format
  column_to_rownames(var = "Sample")  # Set Sample as row names

# Remove the columns (SNPs) that are NaN in at least 1 sample
rows_freq <- row.names(freq_wide)
freq_wide <- as.data.table(freq_wide)
freq_wide <- freq_wide[, lapply(.SD, function(x) if (any(is.nan(x))) NULL else x), .SDcols = names(freq_wide)]
freq_wide <- as.data.frame(freq_wide)
row.names(freq_wide) <- rows_freq
snps <- as.numeric(ncol(freq_wide))

# Transform data to arcosin
freq_asin <- 2*asin(sqrt(freq_wide))

#######################################

###{r SFS calculation}

# We calculate the site frequency spectrum for each sample, as well as the median allele frequency
SFS <- freq_long[,c("Sample","RO")] %>%
          mutate(RO = ifelse(RO > 0.5, 1 - RO, RO)) %>%
          na.omit()

SFS <- SFS %>%
  group_by(Sample) %>%
  mutate(medRO = median(as.numeric(RO))) %>%
  mutate(allele_count = n()) %>%
  mutate(Title = gsub("_contigs$", "", Sample)) %>%
  ungroup()

median_SFS <- SFS %>%
  group_by(Sample) %>%
  summarize(
    median_SFS = median(as.numeric(RO), na.rm = TRUE),
    Total = sum(as.numeric(RO)),
    No0 = sum(as.numeric(RO) != 0, na.rm = TRUE),
    percent = No0/Total) %>%
    dplyr::rename(name = Sample)

ggplot(median_SFS, aes(x = name, y = percent)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Non-Zero Values per Sample",
       x = "Sample Name",
       y = "Number of variant sites") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Open a graphics device (e.g., PNG)
png(filename = paste0("SFS_",reply,".png"))

# Plot the SFS for each isolate
SFS_plot <- ggplot(SFS, aes(x = RO)) +
                geom_histogram(binwidth = 0.1, boundary = 0, fill = "blue", color = "black") +
                facet_wrap(~ Title, scales = "free_y") +
                geom_vline(aes(xintercept = medRO), color = "red", linetype = "dashed") +
                labs(title = paste0(reply,", SFS (Folded)"), x = "Allele Frequency (0-0.5)", y = "Count") +
                theme_minimal() +
                theme(
                  #axis.title.x = element_blank(),
                  #axis.title.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  #axis.ticks.x = element_blank(),
                  #axis.ticks.y = element_blank()
                ) +
                facet_wrap(~ Title, scales = "free_y", labeller = labeller(Title = function(x) {
                  paste0(x," (",SFS %>% filter(Title == x) %>% pull(allele_count) %>% unique(),")")
                }))

# Save the plot using ggsave
ggsave(filename = paste0("SFS_",reply,".png"), plot = SFS_plot, path = visuals_path, width = 8, height = 6, dpi = 300)

# Close the graphics device
dev.off()

# Plot the SFS also in R
SFS_plot

#######################################


###{r Import metadata}

# Import the metadata of all the isolates and subset those present in the PCA
metadata_isolates <- fread(metadata_isolates_path, header=T)
metadata_isolates <- metadata_isolates %>%
                        dplyr::select(sample, Temperature, Generation) %>%
                        dplyr::rename(name = sample) %>%
                        mutate(name = paste0("i",name))
metadata_isolates$Source <- "Isolate"

# Remove four isolates that I know that are contaminated and have intermediate allele frequencies
if (reply == "s__Lactiplantibacillus_plantarum") {
  metadata_isolates <- metadata_isolates %>%
    filter(!name %in% c("iX79", "iX81", "iX595", "iX371"))
}

metadata_pools <- as.data.frame(rownames(freq_asin)[grepl("^cF|^hF", rownames(freq_asin))])
colnames(metadata_pools) <- "name"
metadata_pools <- metadata_pools %>%
  mutate(
    Temperature = case_when(
      str_detect(name, "F0") ~ "BASE",
      str_detect(name, "^h") ~ "HOT",
      str_detect(name, "^c") ~ "COLD"),
    Generation = str_extract(name, "F\\d+"),
    Source = "Pool")

metadata <- rbind(metadata_isolates,metadata_pools)
metadata <- merge(metadata, median_SFS, by = "name")

#######################################


###{r Calculate PCA}

# Do the PCA 
pca.res <- prcomp(freq_asin, retx=TRUE, center = FALSE, scale. =FALSE)

# The variable pca.res$x contains the coordinates of each sample in each PC
pca_data <- as.data.frame(pca.res$x)

# Add sample information (assuming row names of freq_asin are sample names)
pca_data$name <- rownames(freq_asin)

# How much variance is explained by each PC? 
# The variance explained by each PC is the square of the standard deviations divided by the total variance.
sdev <- pca.res$sdev
variance_explained <- sdev^2 / sum(sdev^2)
percentage_variance_explained <- variance_explained * 100

#######################################


###{r Plot PCAs}

# Merge the eigenvectors and the metadata
pca_merged <- merge(pca_data, metadata, by = "name")

# Remove "_contigs$" from the end of some genomes
pca_merged$name <- gsub("_contigs$", "", pca_merged$name)

pca_pool <- subset(pca_merged, pca_merged$Source=="Pool")
pca_isolate <- subset(pca_merged, pca_merged$Source=="Isolate")

# Open a graphics device (e.g., PNG)
png(filename = paste0("PCA_temp_",reply,".png"))

# Plot the PCA colouring the samples by isolation temperature regime
PCA_temp <- ggplot() +
                geom_point(data = pca_isolate, aes(x = PC1, y = PC2, color = Temperature, shape = Source), size = 2, alpha = 0.9) +
                #geom_text(data = pca_isolate, aes(x = PC1, y = PC2, label = name, color = Temperature), size = 1.5, hjust = -0.5) +
                geom_point(data = pca_pool, aes(x = PC1, y = PC2, color = Temperature, shape = Source), size = 2) +
                geom_text(data = pca_pool, aes(x = PC1, y = PC2, label = name, color = Temperature), size = 1.5, hjust = -0.5) +
                labs(
                  title = paste0(reply," (",snps," variable SNPs)"),
                  x = paste0("PC1 (", round(percentage_variance_explained[1], 2), "% variance)"),
                  y = paste0("PC2 (", round(percentage_variance_explained[2], 2), "% variance)"),
                  color = "Temperature", shape = "Source") +
                scale_colour_manual(values = c("HOT" = "red", "BASE" = "grey", "COLD" = "lightblue", "COLD-CONSTANT" = "blue")) +
                theme_minimal()

# Save the plot using ggsave
ggsave(filename = paste0("PCA_temp_",reply,".png"), plot = PCA_temp, path = visuals_path, width = 6, height = 4, dpi = 300)

# Close the graphics device
dev.off()

# Open a graphics device (e.g., PNG)
png(filename = paste0("PCA_SFS_",reply,".png"))

# Plot the PCA colouring the samples by median AF
PCA_SFS <- ggplot() +
              geom_point(data = pca_isolate, aes(x = PC1, y = PC2, color = median_SFS, shape = Source), size = 2, alpha = 0.4) +
              geom_text(data = pca_isolate, aes(x = PC1, y = PC2, label = name, color = median_SFS), size = 1.5, hjust = -0.5) +
              geom_point(data = pca_pool, aes(x = PC1, y = PC2, color = median_SFS, shape = Source), size = 2) +
              geom_text(data = pca_pool, aes(x = PC1, y = PC2, label = name, color = median_SFS), size = 1.5, hjust = -0.5) +
              labs(
                title = paste0(reply," (",snps," variable SNPs)"),
                x = paste0("PC1 (", round(percentage_variance_explained[1], 2), "% variance)"),
                y = paste0("PC2 (", round(percentage_variance_explained[2], 2), "% variance)"),
                color = "median_SFS", shape = "Source") +
              scale_color_gradient(low = "purple", high = "yellow") +
              theme_minimal()

# Save the plot using ggsave
ggsave(filename = paste0("PCA_SFS_",reply,".png"), plot = PCA_SFS, path = visuals_path, width = 6, height = 4, dpi = 300)

# Close the graphics device
dev.off()

# Plot here the PCAs as well
PCA_temp
PCA_SFS

#######################################

