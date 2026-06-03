# Modeling Allele Frequency Distributions in Human Populations

This repository contains the final project for the **Statistical Data Analysis** course. The project implements a statistical inference pipeline to fit, evaluate, and compare parametric distributions against empirical human allele frequency (AF) data sourced from the **gnomAD** database across different ancestral backgrounds.

## Project Overview

The core objective of this study is to analyze how varying evolutionary selection pressures shape the frequency distribution of human genomic variants across different loci. We focus on three genes with highly distinct biological profiles and evolutionary contexts:
* **APOE & CYP2C19**: Highly conserved loci governed primarily by purifying selection, resulting in an abundance of ultra-rare variants.
* **HLA-B**: A highly polymorphic locus belonging to the major histocompatibility complex (MHC), governed by balancing selection to maintain immune diversity.

---

## What Was Done

1.  **Maximum Likelihood Estimation (MLE)**:
    * Fitted two-parameter **Beta** and **Lognormal** distributions to overall population data.
    * Executed a stratified analysis across **10 non-overlapping global sub-populations** to observe shifts in the shape parameters ($\alpha$ and $\beta$).
2.  **Theoretical Population Genetics Alignment**:
    * Simulated the theoretical stationary distribution of the neutral **Wright-Fisher model** to evaluate the geometric mismatch between idealized models and real-world census sequencing data.
3.  **Goodness-of-Fit & Model Diagnostics**:
    * Evaluated model performance quantitatively using the **Kolmogorov-Smirnov (KS) Test** and compared observed empirical moments (mean and variance) against theoretical expected values.
    * Diagnosed local structural blind spots (such as zero-inflation effects) utilizing **Empirical Cumulative Distribution Function (ECDF)** overlays and **Quantile-Quantile (Q-Q) plots**.
4.  **Cross-Platform Visualization Review**:
    * Analyzed visual rendering trade-offs between environments. Identified log-scale coordinate warping issues within R's graphing engine and successfully resolved them by engineering clear, unwarped parameter scatter grids in Python.

---

## Code Architecture (Dual-Language Implementation)

This project contains two fully completed parallel architectures, allowing the entire statistical pipeline to be run, replicated, and verified in both language ecosystems:

* 📂 **`R_version/`**: Core pipeline implemented using `tidyverse`, `fitdistrplus`, `gridExtra`, and `ggplot2` for robust data processing, summary moment-matching, and formal KS-testing.
* 📂 **`Python_version/`**: Parallel implementation leveraging `pandas`, `scipy.stats`, and `matplotlib` to handle numerical data arrays and render unwarped, publication-quality parameter space annotations without geometric distortion.

The underlying maximum likelihood optimization routines achieve absolute mathematical consensus between both environments, providing an informative benchmark on cross-platform visualization behavior.

## Project Report
📄 [Read the full project report here](report/Statistical_Data_Analysis_Final_Project.pdf)