# Modeling Allele Frequency Distributions in Human Populations

This repository contains the final project for the **Statistical Data Analysis** course. The project implements a statistical inference pipeline to fit, evaluate, and compare parametric distributions against empirical human allele frequency (AF) data sourced from the **gnomAD** database across different ancestral backgrounds.

# Background: Modeling Allele Frequency Distributions in Human Populations

Genetic variation is the foundation of biological diversity. In human populations, much of this variation occurs at specific locations in the genome known as genetic variants, often in the form of Single Nucleotide Polymorphisms (SNPs). Each variant is characterized by an **allele frequency (AF)**, which measures how common a particular genetic variant is within a population. For example, if a variant is observed in 1 out of 100 individuals, its allele frequency is 0.01.

## Why Study Allele Frequency Distributions?

Allele frequencies are not uniformly distributed across the genome. Most genetic variants are extremely rare, while only a small fraction reach high frequencies. These patterns emerge from fundamental evolutionary forces such as mutation, natural selection, and genetic drift.

Studying allele frequency distributions is important because it:

- Provides insight into genetic diversity within and between populations.
- Helps identify variants that may influence disease susceptibility or drug response.
- Supports genetic and biomedical research, including personalized medicine.
- Forms the basis for realistic synthetic genome simulation and population genetics modeling.

## Why Is This Useful?

Accurate statistical models of allele frequency distributions are essential for:

- Simulating realistic synthetic genomic datasets.
- Designing and evaluating algorithms for rare variant detection.
- Understanding the evolutionary processes that shape human genetic variation.

---

## What Was Done


## What Was Done

1. **Data Preparation**
   - Collected allele frequency data for APOE, HLA-B, and CYP2C19 from gnomAD.
   - Removed missing values and filtered invalid allele frequencies.
   - Created both global and population-specific datasets.

2. **Exploratory Data Analysis**
   - Calculated summary statistics and quantiles.
   - Visualized allele frequency distributions using histograms, ECDFs, and boxplots.
   - Compared distributions across genes and populations.

3. **Beta Distribution Modeling**
   - Fitted Beta distributions to allele frequency data using Maximum Likelihood Estimation (MLE).
   - Estimated shape parameters (\(\alpha\) and \(\beta\)) for each gene and population.

4. **Model Evaluation**
   - Assessed goodness-of-fit using the Kolmogorov–Smirnov (KS) test.
   - Compared empirical and theoretical distributions using Q–Q plots.
   - Examined deviations between the Beta model and observed data.

5. **Population Analysis**
   - Repeated the fitting procedure for multiple ancestry groups.
   - Compared fitted parameters across populations.
   - Investigated population-specific allele frequency patterns.

6. **Interpretation**
   - Related the results to population genetics concepts, including mutation, genetic drift, and selection.
   - Compared observed patterns with expectations from the Wright–Fisher framework.
   - Discussed limitations of the Beta distribution for modeling real genomic data.
---

## Code Architecture (Dual-Language Implementation)

This project contains two fully completed parallel architectures, allowing the entire statistical pipeline to be run, replicated, and verified in both language ecosystems:

* 📂 **`R_version/`**: Core pipeline implemented using `tidyverse`, `fitdistrplus`, `gridExtra`, and `ggplot2` for robust data processing, summary moment-matching, and formal KS-testing.
* 📂 **`Python_version/`**: Parallel implementation leveraging `pandas`, `scipy.stats`, and `matplotlib` to handle numerical data arrays and render unwarped, publication-quality parameter space annotations without geometric distortion.

## Project Report
📄 [Read the full project report here](Report/gene_af_report_R.pdf)