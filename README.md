# Molecular Additivity of QTNs

This repository contains codes for performing transcriptomics and proteomics analysis to study the molecular additivity of quantitative trait nucleotides (QTNs).

## Repository Structure

The repository is organized into the following main folders:

1. **Temporal_transcriptomic_analysis**
2. **Proteomic_analysis**
3. **Metabolic_modeling**

### Description of Folders

#### 1. Temporal_transcriptomic_analysis
This folder contains an R script `LRT_feature_counts.R`, used to perform differential gene expression analysis using the DESeq2 Likelihood Ratio Test (LRT) method. The LRT method helps identify genes that show significant changes in expression across different conditions. For a detailed explanation of the LRT method, please take a look at [this resource](https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html).

#### 2. Proteomic_analysis
This folder is further divided into two subfolders:
- **mRNA-protein_correlation_analysis**: Contains scripts for performing correlation analysis between mRNA and protein abundance data.
- **protein_abundance_allocation_analysis**: Includes scripts for performing differential protein abundance analysis and protein allocation analysis to assess how proteins are distributed across different cellular processes.

#### 3. Metabolic_modeling
This folder contains:
- `Build_iMAT.MAT`: A MATLAB file for generating context-specific metabolic models using the iMAT algorithm based on proteome data.
- `jaccard index and optgp sampling.ipynb`: A Jupyter notebook for calculating the Jaccard index and performing OptGP sampling to analyze model similarity and metabolic flux distributions.

---


### Acknowledgement
* [Centre for Integrative Biology and Systems medicinE](https://ibse.iitm.ac.in/)
* [Robert Bosch Centre for Data Science and Artificial Intelligence (RBCDSAI)](https://rbcdsai.iitm.ac.in/)
