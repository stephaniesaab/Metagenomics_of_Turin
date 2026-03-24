#Alpha and Beta diversity analysis and visualizations of Metagenomics
#Load libraries ====
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(biomformat)
library(plotly)
library(forcats)
library(ggpubr)
library(ALDEx2)

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
    dplyr::select(name, new_est_reads) %>% 
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

#Get rarefaction curve
otu_table <- as.data.frame(t(otu_table(physeq)))
rare_curve <- rarecurve(otu_table, step = 10000, label = TRUE,
                        main = "Rarefaction Curves by Sample",
                        col = as.numeric(as.factor(metadata_table$Diet)))

#3) Phylum Plot ====

#Read the BIOM file
biom_data <- read_biom("../kraken_bioms/merged_samples.biom")
phybiom <- import_biom(biom_data)

#Clean taxonomy names (Kraken names them k__, p__...)
#Phylum column
colnames(tax_table(phybiom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_table(phybiom) <- gsub("[a-z]__", "", tax_table(phybiom))

#Transform to relative abundance (percentages)
phybiom_rel <- transform_sample_counts(phybiom, function(x) x / sum(x)*100)

#Filter to Phylum level (all should be in kingdom bacteria) to keep it clean
#Merge by phylum
phybiom_class <- tax_glom(phybiom_rel, taxrank = "Order")

sun_df <- psmelt(phybiom_class) %>% 
  filter(Abundance > 5) %>% #Remove low-abundance (< 5%) species
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
  labs(title = "Top 10 Species by Relative Abundance", y = "Relative Abundance (0-1)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Abundance by sample ====
ggplot(df_abundance, aes(x = Sample, y = Abundance, fill = Species)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Paired") + 
  labs(title = "Top 10 Species by Relative Abundance", y = "Relative Abundance (0-1)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#Alpha diversity ====

#Quick view
plot_richness(physeq, x = "Diet")

#View by diet (omnivore of Vegan) -> Vegans have lower diversity measures in both indices
plot_richness(physeq, x="Diet", measures=c("Shannon", "Simpson")) + 
  geom_boxplot() + 
  labs(title = "Alpha Diversity Comparison")

#View by sample 
plot_richness(physeq, measures = c("Shannon", "Simpson")) +
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

# Enhanced Alpha Diversity Plot
plot_richness(physeq, x="samples", measures=c("Observed", "Shannon", "Simpson"), color="Diet") + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_jitter(width = 0.2, size = 3) + 
  theme_bw() +
  labs(title = "Alpha Diversity Comparison: Vegan vs. Omnivore Gut Microbiomes")

#Add plot with t-test values between Diets
alpha_plot <- plot_richness(physeq, x = 'Diet', measures = c("Observed", "Shannon", "Simpson"), color = "Diet") +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(width = 0.2) +
  
  #Adding t-test values in the plot
  stat_compare_means(
    comparisons = list(c("Omnivore", "Vegan")),
    method = "t.test",
    label = "p.signif" #Use symbols (*, **, ns)
  )
  theme_bw() + 
  labs(title = "Alpha Diversity by Diet Group")

#Add the t-test values to the plot
alpha_plot + stat_compare_means(method = "t.test", label.y = 3.5) #Get the plot with the t-test p-values
alpha_plot #Get the plot with the symbols only


#Beta diversity  ====

## Bray-Curtis PCoA (Weighted by abundance) ====
bray_dist <- phyloseq::distance(physeq, method = "bray")
ord_pcoa_bray <- ordinate(physeq, method="PCoA", distance = bray_dist)

#Get PCoA coordinates manually because plot_ordination wasn't capturing Diet to colour points
pcoa_df <- as.data.frame(ord_pcoa_bray$vectors[, 1:2]) #Get PC1 and PC2
colnames(pcoa_df) <- c("PC1", "PC2")

#Add diet metadata manually
pcoa_df$Diet <- as.character(sample_data(physeq)[rownames(pcoa_df), ][["Diet"]])
#Double check samples and diets match the physeq data
stopifnot(all(rownames(pcoa_df) == rownames(sample_data(physeq)[rownames(pcoa_df), ])))

#Get variance for axis labels
var_pc1 <- round(ord_pcoa_bray$values$Relative_eig[1] * 100, 1)
var_pc2 <- round(ord_pcoa_bray$values$Relative_eig[2] * 100, 1)

#Get PCoA plot
bc_plot <- ggplot(data = pcoa_df,
       aes(x = PC1, y = PC2, color = Diet, shape = Diet))+
  geom_point(size = 5, alpha = 0.8) + 
  labs(
    title = "Beta Diversity: PCoA of Gut Microbiomes of Vegans and Omnivores",
    x = paste0("PC1 (", var_pc1, "%)"), 
    y = paste0("PC2 (", var_pc2, "%)")
  ) +
  theme_bw()

#Plot Beta diversity Bray-Curtis distance
plot_ordination(
  physeq, 
  ord_pcoa_bray, 
  color="Diet") +
  geom_point(size = 5) +
  labs(title = "Beta Diversity: PCoA of Bray-Curtis Distances", 
       x = paste0("PC1 (", var_pc1, "%)"), 
       y = paste0("PC2 (", var_pc2, "%)"))+
  theme_minimal()

## NMDS Bray-Curtis (uses stress to represent how well the plot represents the distances) ====
# Trymax ensures the algorithm tries hard enough to find the best fit
ord_nmds_bray <- ordinate(physeq, method = "NMDS", distance = bray_dist, trymax = 100)

#Get Stress value
nmds_stress <- ord_nmds_bray$stress
nmds_stress #Get 7e-05, almost zero

#Get coordinates
nmds_df <- as.data.frame(ord_nmds_bray$points)
colnames(nmds_df) <- c("NMDS1", "NMDS2")

#Mapping sample ID to diet
nmds_df$Diet <- as.character(sample_data(physeq)[rownames(nmds_df), ][["Diet"]])
#Double check samples and diets match the physeq data
stopifnot(all(rownames(nmds_df) == rownames(sample_data(physeq)[rownames(nmds_df), ])))


#First plot
plot_ordination(
  physeq, 
  ord_nmds_bray,
  color = "Diet") + #Not recognizing Diet in metadata
  geom_point(size = 4) +
  labs(title = "NMDS of Bray-Curtis Distances")


#Get NMDS plot with colour
nmds_plot <- ggplot(data = nmds_df,
       aes(x = NMDS1, y = NMDS2, color = Diet, shape = Diet))+
  geom_point(size = 5, alpha = 0.8) +
  labs(
    title = "Beta Diversity: NMDS of Gut Microbiomes of Vegans and Omnivores",
    subtitle = paste("Bray-Curtis Distance (Stress =", nmds_stress, ")"),
    caption = "Figure 4: Non-metric Multidimensional Scaling. Stress < 0.2 indicates a good fit."
  ) +
  theme_bw()

## Jaccard PCoA (Binary: Presence / Absence matrix) ====
#Use binary = TRUE because Jaccard uses binary matrix
jac_dist <- phyloseq::distance(physeq, method = "jaccard", binary = TRUE)
ord_pcoa_jaccard <- ordinate(physeq, method = "PCoA", distance = jac_dist)

#Get PCoA coordinates
jac_df <- as.data.frame(ord_pcoa_jaccard$vectors[, 1:2])
colnames(jac_df) <- c("PC1", "PC2")

#Map metadata manually
jac_df$Diet <- as.character(sample_data(physeq)[rownames(jac_df), ][["Diet"]])

#Get % variance
jac_pc1 <- round(ord_pcoa_jaccard$values$Relative_eig[1] * 100, 1)
jac_pc2 <- round(ord_pcoa_jaccard$values$Relative_eig[2] * 100, 1)

#Make plot
plot_ordination(
  physeq, 
  ord_pcoa_jaccard,
  color = "Diet") +
  geom_point(size = 4) +
  labs(title = "PCoA of Jaccard Distances")

#Get Jaccard plot
jac_plot <- ggplot(data = jac_df,
       aes(x = PC1, y = PC2, color = Diet, shape = Diet))+
  geom_point(size = 5, alpha = 0.8) +
  labs(
    title = "Beta Diversity: Jaccard of Gut Microbiomes of Vegans and Omnivores",
    x = paste0("PC1 (", jac_pc1, "%)"), 
    y = paste0("PC2 (", jac_pc2, "%)")
  ) +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

bc_plot | nmds_plot | jac_plot
## Permanova : capture results to interpret clustering ====
metadata <- as(sample_data(physeq), "data.frame")
permanova_res <- adonis2(bray_dist ~ Diet, data = metadata)
print(permanova_res)

#Get R-squared and pvalue
r2_val <- round(permanova_res$R2[1], 3)
p_val  <- permanova_res$`Pr(>F)`[1]

#Differential Abundance ====

#Extract counts (raw integers needed for ALDEx2)
counts_matrix <- as.matrix(otu_table(physeq))
diets <- as.character(sample_data(physeq)$Diet)

#Run ALDEx2 core (Monte Carlo Sampling)
aldex_out <- aldex(counts_matrix, diets, mc.samples = 128, test = "t", effect = TRUE)

#Run CLR
clr <- aldex.clr(counts_matrix, diets)
effect_obj <- aldex.effect(clr)

#Get effect size results
effect_df <- data.frame(effect_obj)


#Get results
# 'we.eBH' is the Expected Benjamini-Hochberg p-value
# 'diff.btw' is the median difference between groups
aldex_res <- aldex_out %>% 
  rownames_to_column(var = "Species") %>% 
  arrange(we.eBH)
head(aldex_res[order(aldex_res$we.eBH), ])

aldex_res2 <- effect_obj %>% 
  rownames_to_column(var = "Species") %>% 
  arrange(effect)
head(aldex_res2[order(aldex_res2$effect), ])

# abs(effect) > 2, difference between groups is larger than the noise within groups
significant_taxa <- aldex_res2 %>% 
  filter(abs(effect) > 2) %>%
  dplyr::select(Species, effect, diff.btw, diff.win) %>%
  arrange(desc(abs(effect)))

#View the list of significant taxa based on effect size
print(significant_taxa)

#Checking order of levels -> Omnivore is reference (alphabetical)
#Effect size direction: Positive = enriched in vegan, negative = enriched in omnivore
levels(as.factor(diets))

#Plot differential abundance
# 'diff.btw' is difference between groups; 'diff.win' is variation within
ggplot(aldex_res, aes(x = diff.win, y = diff.btw, 
                      color = abs(effect) > 2)) + #All are not significant, colour by effect size instead
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = c("black", "red"), 
                     labels = c("|Effect| < 1 (Noisy)", "|Effect| > 1 (Biologically Significant)")) +
  theme_bw() +
  labs(
    title = "ALDEx2 Differential Abundance: Vegan vs. Omnivore",
    subtitle = "Red points indicate trends toward Diet-specific enrichment",
    x = "Median Dispersion (Within-group)",
    y = "Median Difference (Between-group)",
    color = "Biological Significance"
  )

#Simple approach: Wilcoxon test on top 20 species with p-value adjustment
species_counts <- as.data.frame(t(otu_table(physeq_rel)))
species_counts$Diet <- metadata_table$Diet

#Tests and BH Adjustment
p_vals <- sapply(top20, function(sp) {
  wilcox.test(species_counts[[sp]] ~ species_counts$Diet)$p.value
  })
adj_p <- p.adjust(p_vals, method = "BH")
min(adj_p) #0.8

#Species with the lowest raw p-value
best_sp <- names(which.min(p_vals)) #"Alistipes onderdonkii"

#Boxplot for the most different taxon
ggplot(species_counts, aes(x = Diet, y = .data[[best_sp]], fill = Diet)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2) +
  labs(title = paste("Abundance of", best_sp),
       subtitle = paste("Raw p-value:", round(min(p_vals), 4), "(Not significant after BH adjustment)"),
       y = "Relative Abundance") +
  theme_classic()
