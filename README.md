# TRECOR: Bayesian Covariance Regression for Microbiome Networks

This repository contains the R code and reproducible scripts for the method **TRECOR** (Tree-based Covariance Regression), as presented in the manuscript:

> [cite_start]**BAYESIAN COVARIANCE REGRESSION FOR DIFFERENTIAL NETWORK
ANALYSIS OF ZERO-INFLATED MICROBIOME DATA** [cite: 1]  
> [cite_start]*Zichun Xu and Jing Ma* [cite: 2]  

## Overview

[cite_start]TRECOR is a Bayesian covariance regression framework designed specifically for zero-inflated compositional count data[cite: 45]. [cite_start]By representing microbiome counts through a latent multivariate normal distribution defined on the internal nodes of a phylogenetic tree, TRECOR allows both the mean and the covariance of the latent variables to depend on clinical or environmental covariates[cite: 46]. 

## Repository Structure

* **`R/`**: Contains all base functions and the core algorithm.
  * `Gibbs_R_LTNM.R`: The primary implementation of the TRECOR Gibbs sampling algorithm.
* **`scripts/`**: Executable R scripts for running the models. 
  * [cite_start]Includes scripts used to fit the model on simulated datasets and the real gut microbiome cohort (from Yatsunenko et al., 2012)[cite: 15].
* **`notebook/`**: R Markdown files documenting the downstream analysis.
  * Data quality control and preprocessing.
  * Evaluation of posterior estimation and prediction accuracy.
  * Generation of the figures presented in the manuscript (e.g., ROC curves, differential correlation networks, and tree-based visualizations).