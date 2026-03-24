# Metagenomics_of_Turin
Assignment 3 for Genomic Methods of Bioinformatics, assessing metagenomic sequence data of 6 omnivores and vegans in Turin, Italy.

## Table of Contents
* [1.Introduction](#1-introduction)
* [1.1 Diet and the Gut Microbiome](#11-diet-and-the-gut-microbiome)
* [1.2 Shotgun Metagenomics](#12-shotgun-metagenomics)
* [1.3 Turin Cohort](#13-turin-cohort) 
* [1.4 Comparison of Tools](#14-comparison-of-tools)
* [2. Methods](#2-methods)
* [2.1 Data Acquisition](#21-data-acquisition)
* [2.2 QC and Assessment](#22-qc-and-assessment)
* [2.3 Taxonomic Profiling and Abundance Re-estimation](#23-taxonomic-profiling-and-abundance-re-estimation)
* [2.4 Pre-processing](#24-pre-processing)
* [2.5 Alpha and Beta Diversity Analysis](#25-alpha-and-beta-diversity-analysis)
* [2.6 Differential Abundance](#26-differential-abundance)
* [3. Results](#3-results)
* [3.1 QC](#31-qc)
* [3.2 Taxonomic Classification](#32-taxonomic-classification)
* [3.3 Microbial Community Composition](#33-microbial-community-composition)
* [3.4 Alpha and Beta Diversity](#34-alpha-and-beta-diversity)
* [3.5 Differential Abundance](#35-differential-abundance)
* [4. Discussion](#4-discussion)
* [4.1 Variability of *Segatella copri*](#41-variability-of-segatella-copri)
* [4.2 Alpha and Beta Diversity](#42-alpha-and-beta-diversity)
* [4.3 Candidate Biomarkers](#43-candidate-biomarkers)
* [4.4 Conclusion](#44-conclusion)
* [5. References](#5-references)

## 1. Introduction(#1-introduction)
## 1.1 Diet and the Gut Microbiome

The human gut microbiome is a complex metabolic organ essential for host well-being. It performs critical functions including compound breakdown and vitamin biosynthesis(Filippis et al., 2019). When this ecosystem is imbalanced, it reaches “dysbiosis”. Dysbiosis of the gut microbiome is linked to several pathologies, including neurological disorders, obesity, inflammatory bowel disease, and cancer(Filippis et al., 2019). A number of factors can influence the composition and diversity of the gut microbiome including: diet, age, geography, drugs, environmental substances(David et al., 2014). The “Westernization” of global diets is characterized by high fat and protein intake and low fiber. It is proposed to have lead to a loss of microbial diversity in Western populations. Specifically, dietary habits can dictate the dominance of certain enterotypes: Prevotella is typically associated with fiber-rich diets, while Bacteroides abundance is linked to high-protein and high-fat diets. Even short-term macronutrient shifts can alter microbial composition and diversity. Hence, comparing long-term dietary cohorts, like vegans and omnivores, provides insight into human microbe co-evolution and the effects of dietary habits on gut microbiome composition and diversity(David et al., 2014). 

## 1.2 Shotgun Metagenomics

This study utilizes shotgun metagenomic sequencing to achieve species-level taxonomic resolution. Shotgun metagenomics involves extracting and randomly shearing the total DNA of all microbes in a sample. This allows for a comprehensive species profile and provides insights into the diversity of the microbiome. This gives insight into functional potential of bacteria in the gut and helps answer what do these species do instead of just which species are present(Nearing et al., 2022).
Metagenomics has significant methodological trade-offs, including time, costs, and computational resources. The process requires much higher sequencing coverage than 16S analysis, leading to increased costs for sequencing and higher computational time to process reads. Human DNA contamination occurs in about 50-90% of raw sequences, necessitating robust bioinformatic filtering before reads can be analyzed. Inter-study comparisons are often hindered by variations in DNA extraction methods and the high level of inter-individual variation within the gut microbiome, even within a single genus(Wang et al., 2015).

## 1.3 Turin Cohort

This study analyzes fecal metagenomes from six health adults in Turin, Italy. Turin provides a unique geographic context, ensuring that variations in diet are controlled for geographical differences. Turin is an urban center in Northern Italy, with diets representing a blend of traditional Mediterranean cuisines and modern “Westernized” habits. Researchers used Illumina NextSeq 500 to sequence the gut microbiome of 6 subjects from Turin, Italy, that have omnivorous or vegan diets. Total DNA was fragmented and sequenced and human reads were removed from the total reads(Filippis et al., 2019). Examining this cohort will help evaluate whether dietary choices (Vegan vs. Omnivore) can significantly impact microbiome composition and diversity(Filippis et al., 2019).

The primary goal of this study is to utilize a bioinformatic pipeline to compare the gut microbiomes of vegan and omnivore individuals. The goal of this study was to use a bioinformatic pipeline to analyze and compare the metagenomic sequencing data from 6 individuals in Turin, 3 omnivores and 3 vegans. This study aims to: 
1)	Characterize Diversity: Evaluate the Alpha (within-sample diversity, richness/evenness) diversity between the two dietary groups.
2)	Identify community differences: Determine if the diets cluster distinctly via Bray-Curtis (abundance) and Jaccard (presence/absence) metrics.
3)	Identify differentially abundant species: Use ALDEx2 to identify species that are differentially abundant between vegans and omnivores.

## 1.4 Comparison of Tools

The taxonomic profiler was selected based on the tool’s precision and classification sensitivity. Kraken2 and MetaPhlAn are both standard tools used in the literature. Kraken2 uses a k-mer based alignment that classifies reads based on matches of k-mers of length 31-35 bp in a pre-installed or user-built database(Sun et al., 2021). It uses a Lowest Common Ancestor (LCA) approach for reads that may map to multiple taxa, where a read is assigned to the most specific taxonomic node shared by all the genomes matching that k-mer. MetaPhlAn is marker-gene based, it uses a database of clade-specific marker genes and matches sequences that are unique to a clade(Wright et al., 2023). Comparative studies found that Kraken2 classifies a significantly higher proportion of reads than MetaPhlAn3. However, at a default confidence threshold of 0, Kraken2 is prone to false positives. Researchers show that the Shannon diversity index in Kraken2-classified samples decreases steeply as the confidence threshold (user defined, between 0-1) moves from 0 to 0.20, suggesting that this range effectively filters false positives(Wright et al., 2023). They found that Kraken2 classified a higher proportion of reads than MetPhlAn3 and the Shannon diversity is closer to the true diversity for all Kraken2 classifications. Kraken2 was faster than MetaPhlAn3 but used a larger maximum memory. This can be bypassed using a supercomputer such as the Canada Compute clusters, that allow for large RAM and memory usage(Liu et al., 2024; Wright et al., 2023). 

Taxonomic profiling is limited by distinction between sequence abundance and taxonomic abundance. Larger genomes produce more k-mers, leading to overestimation. It is recommended that a Kraken2 pipeline also uses Bracken to apply a Bayesian re-estimation algorithm(Liu et al., 2024). Bracken uses the known k-mer distribution of genomes to redistribute reads from higher taxon levels to species level, providing a more accurate estimate of relative abundance(Lu et al., 2017a).

Microbiome datasets are compositional, and hence the abundance of one taxon is dependent on others. Two standard tools used for differential abundance are ALDEx2 and ANCOM-BC2(Pelto et al., 2025). ALDEx2 handles compositionality by performing a Centered Log-Ratio (CLR) transformation. It uses Monte Carlo Dirichlet sampling to model technical uncertainty of the sequencing process, this is useful for small datasets such as this study (n = 3 per group), as it prevents sampling noise from influencing results and appearing as biological significance. ANCOM-BC2 is a robust, complex tool that relies on bias-correction factors that are less consistent in small cohorts(Pelto et al., 2025). Studies comparing tools for identifying amplicon sequence variants (ASVs) and operational taxonomic units (OTUs) between two groups of 16S rRNA gene datasets. They found that ALDEx2 and ANCOM-II produce the most consistent results across studies and intersect with results of different approaches. Both tools use counts as input and perform a Wilcoxon rank-sum test. Both ALDEx2 and ANCOM-II were more conservative in identifying a small number of ASVs as significant, having higher precision but with potential loss of sensitivity. ALDEx2 uses a Monte Carlo Dirichlet sampling approach which down weights low abundance ASVs, contributing to a conservative approach, while ANCOM-II uses large multiple testing and thus relies on lots of information. Additionally, ANCOM-II results were less consistent than ALDEx2 in identifying the same genera as significantly differentially abundant across case-related datasets. ALDEx2 is a recommended tool for its more conservative method, as it has a higher replication rate of 79%, while ANCOM-BC2 had a much lower rate at 35%(Pelto et al., 2025). Since this study uses only three samples from each diet group, small samples are more prone to false positives due to sampling noise. ALDEx2 uses more conservative methods, which help it maintain consistency across small sample sizes. ANCOM-BC2 is a more complex model that is more prone to noise-effects. Overall, this study uses ALDEx2 for its consistent and conservative approach to identifying differentially abundant taxa. 

The final diversity analysis was done using Phyloseq package in R. It corresponds the feature table, feature annotation, metadata, representative sequences, and evolutionary tree in one S4 class object. It includes microbial diversity analysis which can be used to find the alpha and beta-diversity of the samples. Phyloseq provides microbiome analysis and data filtering, normalization, and diversity calculations(Wen et al., 2023).

# 2. Methods

## 2.1 Data Acquisition

Raw metagenomic sequence data were retrieved from the NCBI Sequence Read Archive (SRA) using the SRA Toolkit (v3.4.0). FastQ files were generated using fasterq-dump with the --split-3 flag to ensure paired-end integrity and isolate orphaned reads to prevent downstream processing errors. High-performance computing tasks were executed on the Narval cluster (Calcul Québec / Compute Canada) using SLURM scheduling for reproducibility, utilizing 8 threads and 16GB of memory for initial file conversion.

## 2.2 QC and Assessment

Raw reads were assessed with FastQC (v0.12.1) as that is the version available on the Narval cluster(Babraham Bioinformatics - FastQC A Quality Control Tool for High Throughput Sequence Data, n.d.). Quality control reports created by FastQC were aggregated into one report using MultiQC(v1.14). The process was given 8 threads. Individual reports from the 12 files were aggregated using MultiQC to generate a single summary report using default parameters. Based on the report by MultiQC (multiqc_report.html), the sequence reads were high quality with lengths of 150bp. Each sample had a high mean quality score along its sequence (Q > 30). Read counts varied between 28M-46M, with lengths of 150 bp, which are standard for shotgun metagenomics(Wen et al., 2023). All samples had a low adapter content (<2.15%). Overall, the sequences all showed high quality status for adapter content, per base sequence quality, per sequence quality scores, GC content, and per base N content.  No trimming or filtering were done as the reads were already short lengths (150bp) and additional trimming may remove important sequences and reduce the number of k-mers available for taxonomic profiling.

## 2.3 Taxonomic Profiling and Abundance Re-estimation

Taxonomic profiling was conducted using Kraken2 (v2.1.6) against a Standard 16GB RefSeq database (dated 2025-10-15) (available at: https://genome- idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20251015.tar.gz). Due to the large memory and computing resources required for Kraken, it was assigned 32 threads and 64GB memory. The database was specified using –db, the –paired flag was used to indicate that the reads were paired, the confidence threshold was set to 0.15 to minimize the chance of false-positives(Wright et al., 2023). The confidence threshold ensures that k-mer matches are only accepted if they reach a significant fraction of the query read, to prioritize precision of sensitivity in the comparisons. No trimming and filtering were done as the data was already high-quality, but as a safety measure the –minimum-base-quality flag was set to 20. The –use-names flag was set to indicate Kraken to include scientific names instead of just taxon IDs as this makes the report more readable. Bracken (v3.0) (Bayesian Re-estimation of Abundance with kraken) was used to redistribute reads form higher taxonomies to species level(Lu et al., 2017b). The database and output repositories were indicated, as well as the .report files produced from Kraken2. The -r flag was set to 150 (reads are 150bp), the -l flag was set to S to indicate estimations at the species level.

## 2.4 Pre-processing

Taxonomic reports were converted to BIOM format using kraken-biom (v1.2.0) with JSON formatting and imported into R (v4.5.1) via the phyloseq (v1.52.0) package(Lu & Salzberg, 2020; McMurdie & Holmes, 2013). The phyloseq object was constructed manually by integrating the OTU (count table), a taxonomic table (species), and the sample metadata (Diet: Omnivore vs. Vegan). Rarefaction curves were generated using the vegan package (v2.7-2) with a step size of 10,000 to visualize sequencing depth saturation(Dixon, 2003). For visualization, data were transformed to relative abundance, with phyla filtered for (>5\%) abundance for hierarchy plots and top 10 species selected for community composition bar plots, for legibility. The species-level abundance estimates from Bracken for all samples were merged. Missing values were replaced by zeroes to make a complete feature table. The phyla were filtered to remove low abundance species (< 5%) for the sunburst plot. Community composition was assessed by transforming data into relative abundances and filtering for the top 10 most abundant species for legibility. All plots were made with ggplot2 (v4.0.2) and t-tests for the alpha diversity plots were run with ggpubr (v0.6.3)(Create Elegant Data Visualisations Using the Grammar of Graphics • Ggplot2, n.d.; Ggplot2 Based Publication Ready Plots • Ggpubr, n.d.).

## 2.5 Alpha and Beta Diversity Analysis

Alpha-diversity was characterized using Observed, Shannon, and Simpson indices estimate_richness() function of the phyloseq package. Distribution normality was verified via the Shapiro-Wilk test, and dietary group differences in Shannon indices were evaluated using a parametric Student’s t-test (ggpubr v0.6.3)(Ggplot2 Based Publication Ready Plots • Ggpubr, n.d.). Beta-diversity was evaluated using two distances using the distance() function of phyloseq and methods “bray” and “jaccard”. Bray-Curtis distance was calculated via PCoA and NMDS (trymax = 100, to ensure the algorithm finds the best fit) to assess abundance-weighted shifts. Jaccard distance was implemented with the binary = TRUE parameter to identify differences in the presence of rare taxa (presence/absence). Statistical significance of dietary clustering was determined via a permanova of the Bray-curtis distances using the adonis2 function of the vegan package and default parameters.

## 2.6 Differential Abundance

To assess the composition of the microbiome, differential abundance was assessed using ALDEx2 (v1.40.0)(Fernandes et al., 2014). The aldex() function was used with a Monte Carlo sampling (128 instances) followed by Wilcoxon Rank-Sum tests with Benjamini-Hochberg (BH) p-value adjustment to control the false discovery rate (FDR). The effect plot was generated by mapping median between-group differences (diff.btw) against median within-group dispersion (diff.win). Taxa were considered significant if they displayed an effect size > 1. A Wilcoxon Rank-sum test was applied to the top 20 most abundant species, with BH adjustment. 

# 3. Results
## 3.1 QC

Based on the report by MultiQC (multiqc_report.html), each sample had a high mean quality score along its sequence (Q > 30) with standard lengths of 150bp(Metagenomics Sequencing Guide, n.d.). Read counts varied between 28M-46M. All samples had a low adapter content (<2.15%) and stable GC content. Overall, the sequences all showed high quality status for adapter content, per base sequence quality, per sequence quality scores, GC content, and per base N content. Hence, no trimming or filtering were done. This preserved the maximum information and number of k-mers for taxonomic classification.

## 3.2 Taxonomic Classification

Taxonomic assignment via Kraken2 yielded a classification rate between 7-11% at the root node (Table 1). This is a conservative classification, which is a result of the 0.15 confidence threshold applied during processing to prioritize precision and minimize false-positive assignments. 

**Table 1: Kraken2 Classification Rate Percentage by Sample and Diet**
| Sample|Diet|Classification of reads (%)|
| :--- | :---: | ---: |
| SRR8146935 | Omnivore | 9.47 |
| SRR8146936 | Omnivore | 11.58 |
| SRR8146938 | Omnivore | 7.13 |
| SRR8146944 |Vegan |	11.02
| SRR8146963 |Vegan |	8.09
| SRR8146968 |Vegan	|9.86

## 3.3 Microbial Community Composition

Figure 1 revealed that the Turin cohort is dominated by the phylum Bacteroidota (59%) followed by Bacillota (12%), with minor contributions from Verrucomicrobiota, Actinomycetota, and Pseudomonadota. At the class level, Bacteroidia was the most prevalent (39%). Species-level analysis identified Segatella copri as the most abundant taxon in both dietary groups. Notably, the top 10 species profiles were nearly identical between Vegan and Omnivores (Figure 2). However, individual samples (Figure 3) showed that S. copri dominance was driven by three samples (two vegans, one omnivore), suggesting that the differences in species abundance are driven by high inter-individual variability independent of dietary patterns. Trends in lower-abundance taxa were seen with Faecalibacterium prausnitzii appearing more abundant in omnivores and Alistipes putredinis in vegans.

## 3.4 Alpha and Beta Diversity

Alpha diversity metrics (Shannon and Simpson indices) showed higher median diversity in omnivores (Figure 4), but did not find statistically significant differences between groups (p = 0.5187). When plotting by sample, it is clear that there is large variation between samples as two omnivore samples appear to have high Shannon and Simpson indices, while one vegan sample has a high Shannon index (> 3) and Simpson index (>0.75) (Figure 5). Observed alpha-diversity showed more similarities between diets. This lack of significance is likely due to the limited statistical power of the cohort size (n = 3) and high individual variation. 

Beta diversity was assessed using two distances. Bray-curtis distance is weighted by abundance. Principal Coordinate Analysis (PCoA) and Non-metric multidimensional scaling (NMDS) ordinations showed overlapping clusters, likely due to similarities driven by the shared dominance of S. copri. The NMDS achieved an excellent fit (Stress = 7 x 10-5), though this is likely due to the small sample size. Jaccard distance accounts for the presence or absence of rare species and revealed clearer group separation, particularly along PC1 (Figure 7). This suggests that the dietary signal in the cohorts is driven by the presence of diet-specific taxa rather than shifts in the abundance of species. PERMANOVA analysis confirmed that dietary group explains approximately 17% of total community variance (R2 = 0.17), but this is not significant (p = 0.4).

## 3.5 Differential Abundance

Following BH adjustment, no taxa reached significance (padj_min = 0.8). Since the sample has low statistical power, biological significance was assessed using effect size (Figure 8). Five species were identified as having a high biological effect (|effect| > 2), where positive values indicate enrichment in vegans and negative values indicate enrichment in omnivores (Table 2). 

|Species|	Effect|	diff.btw|	diff.win|
| :--- | :---: | :---:| ---: |
|*Holdemania massiliensis*|	-3.733167054|	-9.545938632	|2.707845205
|*Ruminococcus bicirculans* |	3.367441867	|6.171349013	|1.893124784
|*Phascolarctobacterium faecium*|	-2.798243892	|-8.410286769	|3.225362334
|*Ruminococcus callidus*|	2.275096889	|10.93511796|	4.917630134
|*Gordonibacter pamelaeae*|	2.202334875|	4.443706221	|1.959022937

# 4. Discussion
## 4.1 Variability of *Segatella copri*

Individual samples showed significant inter-individual variation in S. copri abundance (Figure 3). As a key member of the human gut microbiome, S. copri is defined as a complex of 13 distinct, yet closely related species with vast intra-species genetic variation(Blanco-Míguez et al., 2023). These species often co-exist within a single host. Strains have similar 16S rRNA gene sequences, making them difficult to distinguish using amplicon-based microbiome analyses. By utilizing k-mer based metagenomic classification (Kraken2/Bracken), this study captured the broad presence of this taxon. S. copri species are highly prevalent in non-Westernized populations, with multiple species commonly co-occurring within individuals, and strains can differ between Western and non-Western lifestyles(Blanco-Míguez et al., 2023). Its high variability in abundance in this Turin cohort suggests a transitional microbial state that may be independent of dietary patterns and instead reflect large inter-individual variation  in S. copri (Figure 5). 

## 4.2 Alpha and Beta Diversity

Literature suggests that vegans and vegetarians have significantly greater richness in alpha diversity compared to omnivores, especially for Bacteroidetes OTUs(Tomova et al., 2019). A study on the impact of vegan, vegetarian, and omnivore diets on the gut microbiome of adults in Italy found statistically significant differences in alpha diversity between vegetarians and omnivores, but only marginally significant differences between vegans and omnivores, with diversity being lower in omnivores compared to vegetarians. This trend was not significant in this study (p = 0.5187). The lack of significance in Shannon and Simpson indices likely reflects the limited statistical power of the small cohort (n = 3 per group) and the high inter-individual noise seen in dominant species like S. copri. However, the observed median trend aligns with the hypothesis that dietary fiber consumption in plant-based diets supports niche diversity, though the metrics are obscured by individual variation in this small dataset. 

The discrepancy between the Bray-Curtis and Jaccard distances aligned with the literature on microbiome profiles in Italians. The abundance-weighted Bray-Curtis metric failed to cluster diets, likely due to shared dominance of a few highly abundant taxa. This is similar to results in the Italian cohort as the beta diversity of the samples did not show differences, which the authors attribute to similarities in nutrient intake rather than differences in food items(Losasso et al., 2018). Common factors between the diets was high fat content and reduced protein and carbohydrate intakes, which may drive microbial communities towards a westernized profile. 

The Jaccard distance (presence/absence) showed some distinct clustering (Figure 7). These results align with the Italian cohort as they found differences in composition in the rare microbiota(Losasso et al., 2018). This suggests that the vegan microbiome is defined by the recruitment of rare, diet-specific taxa. Specific food choices may drive the presence of rare bacteria that have distinct functions adapted for the vegan gut. 

## 4.3 Candidate Biomarkers

Despite no significant p-values, the high ALDEx2 effect sizes identified potentially biologically relevant biomarkers (Table 2). A study investigating the causal association between new gut microbiota and Alzheimers disease found that a higher abundance of *Holdemania massiliensis* could increase the risk of AD(Lin et al., 2026). Since this was found to be more prevalent in omnivore diets (Table 2), this suggests that long-term omnivorous diets may harbor pathogens not found in vegan diets. Conversely, enrichment of *Ruminococcus bicirculans* in vegans highlights a functional adaptation to plant-based fibers. *R. bicirculans* was found to have a carbohydrate active enzyme that may be able to use Beta-glucans and xyloglucans for growth(Wegmann et al., 2014). Its presence suggests that the vegan diet actively supports a microbial niche specialized in fermenting complex plant polysaccharides, which are have been shown to support gut health and control cholesterol levels(Wegmann et al., 2014).

## 4.4 Conclusion

This study demonstrates that there may be some microbial differences due to dietary patterns of omnivores and vegans. This study is limited by its small sample size and low statistical power (n = 6). One strength of the pipeline is the use of high-performance computing to utilize resource-intensive programs (Kraken2). Future studies should investigate differences in larger sample sizes to obtain increase statistical power and incorporate metatranscriptomics to determine if diet-specific species actively express specialized enzymes.

## 5. References

*Babraham Bioinformatics—FastQC A Quality Control tool for High Throughput Sequence Data*. (n.d.). Retrieved March 15, 2026, from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Blanco-Míguez, A., Gálvez, E. J. C., Pasolli, E., De Filippis, F., Amend, L., Huang, K. D., Manghi, P., Lesker, T.-R., Riedel, T., Cova, L., Punčochář, M., Thomas, A. M., Valles-Colomer, M., Schober, I., Hitch, T. C. A.,
Clavel, T., Berry, S. E., Davies, R., Wolf, J., … Strowig, T. (2023). Extension of the Segatella copri complex to 13 species with distinct large extrachromosomal elements and associations with host conditions. *Cell Host & Microbe*, 31(11), 1804-1819.e9. https://doi.org/10.1016/j.chom.2023.09.013

*Create Elegant Data Visualisations Using the Grammar of Graphics • ggplot2*. (n.d.). Retrieved March 15, 2026, from https://ggplot2.tidyverse.org/

David, L. A., Maurice, C. F., Carmody, R. N., Gootenberg, D. B., Button, J. E., Wolfe, B. E., Ling, A. V., Devlin, A. S., Varma, Y., Fischbach, M. A., Biddinger, S. B., Dutton, R. J., & Turnbaugh, P. J. (2014). Diet rapidly and reproducibly alters the human gut microbiome. *Nature*, 505(7484), 559–563. https://doi.org/10.1038/nature12820

Dixon, P. (2003). VEGAN, a package of R functions for community ecology. *Journal of Vegetation Science*, 14(6), 927–930. https://doi.org/10.1111/j.1654-1103.2003.tb02228.x

Fernandes, A. D., Reid, J. N., Macklaim, J. M., McMurrough, T. A., Edgell, D. R., & Gloor, G. B. (2014). Unifying the analysis of high-throughput sequencing datasets: Characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome*, 2, 15. https://doi.org/10.1186/2049-2618-2-15

Filippis, F. D., Pasolli, E., Tett, A., Tarallo, S., Naccarati, A., Angelis, M. D., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets. *Cell Host & Microbe*, 25(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004

*Ggplot2 Based Publication Ready Plots • ggpubr*. (n.d.). Retrieved March 15, 2026, from https://rpkgs.datanovia.com/ggpubr/

Lin, F., Yang, X., Xie, X., Qiu, L., Zheng, C., Zheng, X., Li, Z., & Nie, B. (2026). The causal association between the new gut microbiota and Alzheimer’s disease: A Mendelian randomization study. *Journal of Alzheimer’s Disease Reports*, 10, 25424823261422629. https://doi.org/10.1177/25424823261422629

Liu, Y., Ghaffari, M. H., Ma, T., & Tu, Y. (2024). Impact of database choice and confidence score on the performance of taxonomic classification using Kraken2. *aBIOTECH*, 5(4), 465–475. https://doi.org/10.1007/s42994-024-00178-0

Losasso, C., Eckert, E. M., Mastrorilli, E., Villiger, J., Mancin, M., Patuzzi, I., Di Cesare, A., Cibin, V., Barrucci, F., Pernthaler, J., Corno, G., & Ricci, A. (2018). Assessing the Influence of Vegan, Vegetarian and Omnivore Oriented Westernized Dietary Styles on Human Gut Microbiota: A Cross Sectional Study. *Frontiers in Microbiology*, 9, 317. https://doi.org/10.3389/fmicb.2018.00317

Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: Estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104. https://doi.org/10.7717/peerj-cs.104

Lu, J., & Salzberg, S. L. (2020). Ultrafast and accurate 16S rRNA microbial community analysis using Kraken 2. *Microbiome*, 8, 124. https://doi.org/10.1186/s40168-020-00900-2

McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. *PLoS ONE*, 8(4), e61217. https://doi.org/10.1371/journal.pone.0061217

*Metagenomics Sequencing Guide*. (n.d.). Retrieved March 15, 2026, from https://genohub.com/shotgun-metagenomics-sequencing/

Nearing, J. T., Douglas, G. M., Hayes, M. G., MacDonald, J., Desai, D. K., Allward, N., Jones, C. M. A., Wright, R. J., Dhanani, A. S., Comeau, A. M., & Langille, M. G. I. (2022). Microbiome differential abundance methods produce different results across 38 datasets. *Nature Communications*, 13, 342. https://doi.org/10.1038/s41467-022-28034-z

Pelto, J., Auranen, K., Kujala, J. V., & Lahti, L. (2025). Elementary methods provide more replicable results in microbial differential abundance analysis. *Briefings in Bioinformatics*, 26(2), bbaf130. https://doi.org/10.1093/bib/bbaf130

Sczyrba, A., Hofmann, P., Belmann, P., Koslicki, D., Janssen, S., Dröge, J., Gregor, I., Majda, S., Fiedler, J., Dahms, E., Bremges, A., Fritz, A., Garrido-Oter, R., Jørgensen, T. S., Shapiro, N., Blood, P. D., Gurevich, A., Bai, Y., Turaev, D., … McHardy, A. C. (2017). Critical Assessment of Metagenome Interpretation – a benchmark of computational metagenomics software. *Nature Methods*, 14(11), 1063–1071. https://doi.org/10.1038/nmeth.4458

Sun, Z., Huang, S., Zhang, M., Zhu, Q., Haiminen, N., Carrieri, A. P., Vázquez-Baeza, Y., Parida, L., Kim, H.-C., Knight, R., & Liu, Y.-Y. (2021). Challenges in Benchmarking Metagenomic Profilers. *Nature Methods*, 18(6), 618–626. https://doi.org/10.1038/s41592-021-01141-3

Tomova, A., Bukovsky, I., Rembert, E., Yonas, W., Alwarith, J., Barnard, N. D., & Kahleova, H. (2019). The Effects of Vegetarian and Vegan Diets on Gut Microbiota. *Frontiers in Nutrition*, 6, 47. https://doi.org/10.3389/fnut.2019.00047

Wang, W.-L., Xu, S.-Y., Ren, Z.-G., Tao, L., Jiang, J.-W., & Zheng, S.-S. (2015). Application of metagenomics in the human gut microbiome. *World Journal of Gastroenterology : WJG*, 21(3), 803–814. https://doi.org/10.3748/wjg.v21.i3.803

Wegmann, U., Louis, P., Goesmann, A., Henrissat, B., Duncan, S. H., & Flint, H. J. (2014). Complete genome of a new Firmicutes species belonging to the dominant human colonic microbiota (‘Ruminococcus bicirculans’) reveals two chromosomes and a selective capacity to utilize plant glucans. *Environmental Microbiology*, 16(9), 2879–2890. https://doi.org/10.1111/1462-2920.12217

Wen, T., Niu, G., Chen, T., Shen, Q., Yuan, J., & Liu, Y.-X. (2023). The best practice for microbiome analysis using R. *Protein & Cell*, 14(10), 713–725. https://doi.org/10.1093/procel/pwad024

Wright, R. J., Comeau, A. M., & Langille, M. G. I. (2023). From defaults to databases: Parameter and database choice dramatically impact the performance of metagenomic taxonomic classification tools. *Microbial Genomics*, 9(3), 000949. https://doi.org/10.1099/mgen.0.000949
