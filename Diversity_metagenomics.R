#Alpha and Beta diversity analysis and visualizations of Metagenomics

#TODO ====
# 1) Visualization of Taxonomic abundance (stacked bar plot or circumference plot)
# 2) Alpha diversity: Shannon index and something else?
# 3) Beta Diversity : PCoA plot (Bray-Curtis distance)
# 4) Differential Abundance: Compare groups
# 5) Fix plot of abundance to show differences between omnivores and vegans, label X axis with diet
# 6) Colour alpha diversity plot by diet
# 7) T tests don't take into account compositionality in the same way?
# 8) Run the permanova
# 9) try to get phylum level and plot genus and phylum levels
# 10) To plot the relative abundance you want to do some statistics to filter for statistically significant ones, then get top 20 and plot them, find the significant ones and then do a box and whisker plot for significant species
# 11) Converting to BIOM is easier -> Not available as a cluster in narval
# 12) Use the species report files to get the taxonomic levels and get the phylum plots
# 13) Need to adjust p values?

# Bray-curtis takes into account abundance, while Jaccard only cares about presence/absence. 
# What does it mean if they both do a good job of splitting the data?
##### ====
#There's something weird here because three samples have a lot of sagatella copri 
#There's individual variation is very high 



#---
#Load libraries ====
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(biomformat)
library(plotly)
library(forcats)

#1) Making Metadata table ====
metadata_table <- data.frame(
  Accession = c("SRR8146938", "SRR8146936", "SRR8146935", "SRR8146968", "SRR8146963", "SRR8146944"),
  Diet = c("Omnivore", "Omnivore", "Omnivore", "Vegan", "Vegan", "Vegan")
)
metadata_table <- column_to_rownames(metadata_table, "Accession")

#Sample column matches 
SAM = sample_data(metadata_table)

#2) Manually building phyloseq object ====
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

## OTU Table (Counts) ====
#Make species as row-names to make a numeric matrix
counts_matrix <- abundance_table %>% 
  column_to_rownames("name") %>% 
  as.matrix()

#Make the OTU table, the taxa are row names so set taxa_are_rows = TRUE
OTU = otu_table(counts_matrix, taxa_are_rows = TRUE)

## Taxonomy table ====

#Previously set bracken results to list species names
#Create taxonomy table where species is the name
tax_matrix <- matrix(rownames(counts_matrix), nrow = nrow(counts_matrix), ncol = 1)
rownames(tax_matrix) <- rownames(counts_matrix) #Set species as row-names
colnames(tax_matrix) <- "Species"
TAX = tax_table(tax_matrix)

## Make phyloseq object ====
physeq <- phyloseq(OTU, TAX, SAM)

#See the object
physeq

otu_table <- as.data.frame(t(otu_table(physeq)))
rare_curve <- rarecurve(otu_table, step = 10000)

#3) Phylum Plot ====

#Read the BIOM file
biom_data <- read_biom("./kraken_bioms/merged_samples.biom")
phybiom <- import_biom(biom_data)

#Clean taxonomy names (Kraken names them k__, p__...)
#Phylum column
colnames(tax_table(phybiom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_table(phybiom) <- gsub("[a-z]__", "", tax_table(phybiom))
# tax_table(phybiom)[, "Phylum"] <- gsub("p__", "", tax_table(phybiom)[, "Phylum"])
# tax_table(phybiom)[, "Kingdom"] <- gsub("k__", "", tax_table(phybiom)[, "Kingdom"])
# tax_table(phybiom)[, "Class"] <- gsub("c__", "", tax_table(phybiom)[, "Class"])
# tax_table(phybiom)[, "Order"] <- gsub("o__", "", tax_table(phybiom)[, "Order"])
# tax_table(phybiom)[, "Family"] <- gsub("f__", "", tax_table(phybiom)[, "Family"])
# tax_table(phybiom)[, "Genus"] <- gsub("g__", "", tax_table(phybiom)[, "Genus"])
# tax_table(phybiom)[, "Species"] <- gsub("s__", "", tax_table(phybiom)[, "Species"])

#Transform to relative abundance (percentages)
phybiom_rel <- transform_sample_counts(phybiom, function(x) x / sum(x)*100)

#Filter to Phylum level (all should be in kingdom bacteria) to keep it clean
#Merge by phylum
phybiom_class <- tax_glom(phybiom_rel, taxrank = "Order")

sun_df <- psmelt(phybiom_class) %>% 
  filter(Abundance > 5) %>% #Remove low-abundance (< 1%) species
  mutate(Kingdom = "Bacteria") #Make a common root for the center

#Level 1 is kingdom (root)
lvl1 <- sun_df %>% 
  group_by(labels = Kingdom) %>% 
  summarize(values = sum(Abundance), .groups = "drop") %>% 
  mutate(ids = labels, parents = "")

#Level 2 is phyla
lvl2 <- sun_df %>% 
  group_by(ids = Phylum, labels = Phylum, parents = Kingdom) %>% 
  summarize(values = sum(Abundance), .groups = "drop")

#Level 3 is classes  -> want unique ones (Phylum-Class)
lvl3 <- sun_df %>% 
  group_by(Phylum, Class) %>% 
  summarize(values = sum(Abundance), .groups = "drop") %>% 
  mutate(ids = paste(Phylum, Class, sep = " - "), 
         labels = Class,
         parents = Phylum)

#Level 4 is Order
lvl4 <- sun_df %>% 
  group_by(Phylum, Class, Order) %>% 
  summarize(values = sum(Abundance), .groups = "drop") %>% 
  mutate(ids = paste(Phylum, Class, Order, sep = " - "), 
         labels = Order, 
         parents = paste(Phylum, Class, sep = " - "))

#Combine levels into one table
sun_data <- bind_rows(lvl1, lvl2, lvl3, lvl4)

#Make sunburst plot with plotly
sun_fig <- plot_ly(
  data = sun_data,
  ids = ~ids, #Deepest level being Order
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  type = 'sunburst',
  branchvalues = 'remainder', #Because filtered for taxa > 5% abundance
  textinfo = 'label+percent entry', #Show name and % of parent
  insidetextorientation = 'radial', #Curve text for better fit
  hoverinfo = 'label+value+percent root' #Give details when you hover
  ) %>% 
  layout(title = "Taxanomic Hierarchy Merged Samples", 
         margin = list(l = 0, r = 0, b = 0, t = 40))

sun_fig

#Taxanomic abundance (Top 10 most abundant species) ====

#Get relative abundance percentages
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

#Get top 10 most abundant species
#Calculate total abundance of each species across all samples
top10 <- names(sort(taxa_sums(physeq_rel), decreasing = TRUE))[1:10]
top20 <- names(sort(taxa_sums(physeq_rel), decreasing = TRUE))[1:20]

#Subset data to get top 10
physeq_top10 <- prune_taxa(top10, physeq_rel)
physeq_top20 <- prune_taxa(top20, physeq_rel)

#Dataframe for plotting
df_abundance <- psmelt(physeq_top10)
df_abundance20 <- psmelt(physeq_top20)

#Plotting top 10 most abundant species
ggplot(df_abundance, aes(x = Diet, y = Abundance, fill = Species)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Paired") + 
  labs(title = "Top 10 Species Relative Abundance", y = "Relative Abundance (0-1)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Ugly plot ====
ggplot(df_abundance20, aes(x = Sample, y = Abundance, fill = Species)) + 
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

#Extract alpha diversity metrics
alpha_metrics <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson"))

#Add diet data
diet_data <- data.frame(sample_data(SAM))
alpha_data <- cbind(alpha_metrics, diet_data)

#Check normality for Shannon index
shapiro.test(alpha_metrics$Shannon) #p = 0.3717, normal

#Checking between the two diets (Vegan and Omnivore)
t.test(Shannon ~ Diet, data = alpha_data)

#Should do permanova



#gracken --input_dir bracken_reports --bac_taxonomy bac120_taxonomy_r95.tsv --ar_taxonomy ar122_taxonomy_r95.tsv --bac_tree bac120_r95.tree --ar_tree ar122_r95.tree --out_prefix my_tree --taxonomy gtdb

#Beta diversity  ====

#Calculate Bray-Curtis distance
# PCoA with bray-curtis
sample_data(phybiom) <- SAM
bray_dist <- phyloseq::distance(physeq, method = "bray")
ord_pcoa_bray <- ordinate(physeq, method="PCoA")

#Plot Beta diversity Bray-Curtis distance
plot_ordination(
  physeq, 
  ord_pcoa_bray, 
  color="Diet", 
  title="Bray PCoA") +
  geom_point(size = 4) +
  labs(title = "Beta Diversity: PCoA of Bray-Curtis Distances") +
  theme_minimal()

# NMDS with the same distance measure. 
ord_nmds_bray <- ordinate(physeq, method="NMDS", distance="bray")
plot_ordination(physeq, ord_nmds_bray, color="Time", title="Bray NMDS") + geom_point(size = 4)

# Jaccard distance
ord_pcoa_jaccard <- ordinate(physeq, method="PCoA", distance="jaccard")
plot_ordination(physeq, ord_pcoa_jaccard, color="Time", title="Bray PCoA") + geom_point(size = 4)


#Adonis2 call for permanova

metadata <- as(sample_data(physeq), "data.frame")
adonis2(phyloseq::distance(physeq, method = "bray") ~ Diet,
        data = metadata)

#Differential abundance -> ANCOM-BC

#Simple approach: Wilcoxon test on top speciesS with p-value adjustment
species_counts <- as.data.frame(t(otu_table(physeq_rel)))
species_counts$Diet <- metadata_table$Diet

#Tests and BH Adjustment
p_vals <- sapply(top20, function(sp) {
  wilcox.test(species_counts[[sp]] ~ species_counts$Diet)$p.value
  })
adj_p <- p.adjust(p_vals, method = "BH")
min(adj_p) #0.8


#
