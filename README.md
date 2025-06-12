# Molecular Additivity of QTNs

This repository contains codes for performing transcriptomics and proteomics analysis to study the molecular additivity of quantitative trait nucleotides (QTNs).

## Repository Structure

The repository is organized into the following main folders:

1. **Temporal_transcriptomic_analysis**
2. **Proteomic_analysis**
3. **Metabolic_modeling**
4. **Amino acid analysis**

### Description of Folders

#### 1. Temporal_transcriptomic_analysis
This folder contains an R script `LRT_feature_counts.R`, used to perform differential gene expression analysis using the DESeq2 Likelihood Ratio Test (LRT) method. The LRT method helps identify genes that show significant changes in expression across different conditions. For a detailed explanation of the LRT method, please take a look at [this resource](https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html).

#### 2. Proteomic_analysis
The absolute quantification of proteins from raw mass spectrometry (MS) files were done using autoprot [this resource](https://github.com/biosustain/autoprot)

This folder is further divided into two subfolders:
- **mRNA-protein_correlation_analysis**: Contains scripts for performing correlation analysis between mRNA and protein abundance data.
- **protein_abundance_allocation_analysis**: Includes scripts for performing differential protein abundance analysis and protein allocation analysis to assess how proteins are distributed across different cellular processes.

#### 3. Metabolic_modeling
This folder contains:
- `Build_iMAT.MAT`: A MATLAB file for generating context-specific metabolic models using the iMAT algorithm based on proteome data.
- `jaccard index and optgp sampling.ipynb`: A Jupyter notebook for calculating the Jaccard index and performing OptGP sampling to analyze model similarity and metabolic flux distributions.

### 4. Amino Acid Analysis
This folder contains scripts for analyzing intracellular amino acid profiles in yeast strains:
- **`Z_score_heatmap.R`**  
  Generates a Z-score heatmap to cluster amino acids, using the MMTT strain as a reference for comparison across conditions.

- **`trajectory_analysis.R`**  
  Analyzes selected amino acids in the MMTT strain. Focuses on:
  - Amino acids exhibiting a **biphasic trend**.
  - Amino acids showing an **early surge** followed by a **sustained level** during the later phase of sporulation.

---


### Acknowledgement
* [Centre for Integrative Biology and Systems medicinE](https://ibse.iitm.ac.in/)
* [Robert Bosch Centre for Data Science and Artificial Intelligence (RBCDSAI)](https://rbcdsai.iitm.ac.in/)
