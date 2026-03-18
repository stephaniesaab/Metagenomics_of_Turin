#Alpha and Beta diversity analysis and visualizations of Metagenomics

#TODO ====
# 1) Visualization of Taxonomic abundance (stacked bar plot or circumference plot)
# 2) Alpha diversity: Shannon index and something else?
# 3) Beta Diversity : PCoA plot (Bray-Curtis distance)
# 4) Differential Abundance: Compare groups
# 5) Fix plot of abundance to show differences between omnivores and vegans, label X axis with diet
#---
#Load libraries ====
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)

#1) Making Metadata table ====
metadata_table <- data.frame(
  Accession = c("SRR8146938", "SRR8146936", "SRR8146935", "SRR8146968", "SRR8146963", "SRR8146944"),
  Diet = c("Ominvore", "Ominvore", "Ominvore", "Vegan", "Vegan", "Vegan")
)
metadata_table <- column_to_rownames(metadata_table, "Accession")

#Sample column matches 
SAM = sample_data(metadata_table)

#2) Manually building phyloseq object
#Load Kraken and Bracken reports and merge into one table

#Loading .bracken files and merging them into one table

#Get all the bracken file names in the kraken_results folder
bracken_files <- list.files(path = "./kraken_results/", pattern = "*.bracken$", full.names = TRUE)

#Get all the tables, only take the columns with species names and read counts
data_list <- lapply(bracken_files, function(x) {
  accession_name <- basename(x)
  
  read_tsv(x) %>% 
    select(name, new_est_reads) %>% 
    rename(!!accession_name := new_est_reads)
})

#Join the tables by the species column
abundance_table <- reduce(data_list, full_join, by = "name") %>% 
  replace(is.na(.), 0)

#Add relative proportions of counts
relative_abundance <- apply(abundance_table[, -1], 2, function(x) x / sum(x))

#Add back species names
relative_abundance <- cbind(name = abundance_table$name, as.data.frame(relative_abundance))

#Remove .bracken from column names
colnames(abundance_table) <- gsub(".bracken", "", colnames(abundance_table))
colnames(relative_abundance) <- gsub(".bracken", "", colnames(abundance_table))

#Check first 5 rows
head(abundance_table)
head(relative_abundance)

#OTU Table (Counts) ====
#Make species as row-names to make a numeric matrix
counts_matrix <- abundance_table %>% 
  column_to_rownames("name") %>% 
  as.matrix()

#Make the OTU table, the taxa are row names so set taxa_are_rows = TRUE
OTU = otu_table(counts_matrix, taxa_are_rows = TRUE)

#Taxonomy table ====

#Previously set bracken results to list species names
#Create taxonomy table where species is the name
tax_matrix <- matrix(rownames(counts_matrix), nrow = nrow(counts_matrix), ncol = 1)
rownames(tax_matrix) <- rownames(counts_matrix) #Set species as row-names
colnames(tax_matrix) <- "Species"
TAX = tax_table(tax_matrix)

#Make phyloseq object ====
physeq <- phyloseq(OTU, TAX, SAM)

#See the object
physeq

otu_table <- as.data.frame(t(otu_table(phyloseq)))
rare_curve <- rarecurve(otu_table, step = 10000)

#Taxanomic abundance (Top 10 most abundant species) ====

#Get relative abundance percentages
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

#Get top 10 most abundant species
#Calculate total abundance of each species across all samples
top10 <- names(sort(taxa_sums(physeq_rel), decreasing = TRUE))[1:10]

#Subset data to get top 10
physeq_top10 <- prune_taxa(top10, physeq_rel)

#Dataframe for plotting
df_abundance <- psmelt(physeq_top10)

#Plotting top 10 most abundant species
ggplot(df_abundance, aes(x = Sample, y = Abundance, fill = Species)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Paired") + 
  labs(title = "Top 10 Species Relative Abundance", y = "Relative Abundance (0-1)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Alpha diversity ====

#Quick view
plot_richness(physeq)

#View by diet (omnivore of Vegan) -> Vegans have lower diversity measures in both indices
plot_richness(physeq, x="Diet", measures=c("Shannon", "Simpson")) + 
  geom_boxplot() + 
  labs(title = "Alpha Diversity Comparison")
