# Awesome Repeated Measures "omics" Analysis

A [community-maintained](https://github.com/smdabdoub/awesome-rm-omics/graphs/contributors) list of software packages for analysis of longitudinal/repeated measures studies utilizing high-dimensional "omics" data.

While many of the packages here are marketed for a single type of "omics" data (metagenomics, metabolome, etc.). Many can be utilized across the -omics universe.

The common thread among the methods listed here is that the same samples are measured across multiple time-points and may be measured with two or more modalities (e.g.: microbiome, metabolome, and miRNA). The data can be described as multiple matrices/tables with the same number of samples (not a strict requirement for all methods) and varying number of features.

The repo is in the style of Mike Love's [awesome-multi-omics](https://github.com/mikelove/awesome-multi-omics) repository for multi-omics analysis methods (which, in turn, is in the style of Sean Davis' [awesome-single-cell](https://github.com/seandavi/awesome-single-cell) repo for single-cell analysis methods). Note that there is a Multi-omics section for software handling repeated measures studies with multiple measurement modalities, so there will be some overlap with [awesome-multi-omics](https://github.com/mikelove/awesome-multi-omics).

[Contributions welcome](https://github.com/smdabdoub/awesome-rm-omics/blob/master/CONTRIBUTING.md)...

For brevity, entries for software and scientific articles only list the last name of the first author.

## Software packages and methods

### Spline Modeling

- 2015 - Straube - [A Linear Mixed Model Spline Framework for Analysing Time Course ‘Omics’ Data](https://doi.org/10.1371/journal.pone.0134540)
- 2017 - [**metaDprof**](https://cals.arizona.edu/~anling/software/metaDprof.htm) - Luo - [An informative approach on differential abundance analysis for time-course metagenomic sequencing data](https://doi.org/10.1093/bioinformatics/btw828)
- 2018 - [**SplinectomeR**](https://rrshieldscutler.github.io/splinectomeR/) - Shields-Cutler - [SplinectomeR Enables Group Comparisons in Longitudinal Microbiome Studies](https://doi.org/10.3389/fmicb.2018.00785)
- 2018 - [**MetaLonDA**](https://github.com/aametwally/MetaLonDA) - Metwally - [MetaLonDA: a flexible R package for identifying time intervals of differentially abundant features in metagenomic longitudinal studies](https://doi.org/10.1186/s40168-018-0402-y) - [vignette](https://cran.r-project.org/web/packages/MetaLonDA/vignettes/MetaLonDA.html)

### Mixed Modeling

- 2016 - Chen - zero-inflated Beta regression (ZIBR) - [A two-part mixed-effects model for analyzing longitudinal microbiome compositional data](https://doi.org/10.1093/bioinformatics/btw308)
- 2018 - Tipton - [Measuring associations between the microbiota and repeated measures of continuous clinical variables using a lasso-penalized generalized linear mixed model](https://doi.org/10.1186/s13040-018-0173-9)
- 2020 - [**NBZIMM**](https://github.com/nyiuab/NBZIMM) - Zhang - [NBZIMM: negative binomial and zero-inflated mixed models, with application to microbiome/metagenomics data analysis](https://doi.org/10.1186/s12859-020-03803-z) - [vignette](https://abbyyan3.github.io/NBZIMM-tutorial/)
- 2020 - Zhang - [Fast zero-inflated negative binomial mixed modeling approach for analyzing longitudinal metagenomics data](https://doi.org/10.1093/bioinformatics/btz973)
- 2021 - [**MaAsLin 2**](https://huttenhower.sph.harvard.edu/maaslin/) - Mallick - [Multivariable association discovery in population-scale meta-omics studies](https://doi.org/10.1371/journal.pcbi.1009442) - [vignette](https://github.com/biobakery/biobakery/wiki/maaslin2)
- 2022 - [**MiRKAT-MC**](https://github.com/Zhiwen-Owen-Jiang/MiRKATMC) - Jiang - [MiRKAT-MC: A Distance-Based Microbiome Kernel Association Test With Multi-Categorical Outcomes](https://doi.org/10.3389/fgene.2022.841764) - [vignette](https://cran.r-project.org/web/packages/MiRKAT/vignettes/MiRKAT_Vignette.html)

### Multi-level Modeling

- 2012 - Liquet - [A novel approach for biomarker selection and the integration of repeated measures experiments from two assays](https://doi.org/10.1186/1471-2105-13-325)

### Network Analysis

- 2011 - Xia - [Extended local similarity analysis (eLSA) of microbial community and other time series data with replicates](https://doi.org/10.1186%2F1752-0509-5-S2-S15)
- 2018 - Coenen - [Limitations of Correlation-Based Inference in Complex Virus-Microbe Communities](https://doi.org/10.1128/msystems.00084-18)
- 2019 - Ai - [Constructing the Microbial Association Network from large-scale time series data using Granger causality](https://doi.org/10.3390%2Fgenes10030216)

### Machine Learning/Artificial Neural Networks (Deep Learning) Methods

- 2021 - Sharma - [phyLoSTM: a novel deep learning model on disease prediction from longitudinal microbiome data](https://doi.org/10.1093/bioinformatics/btab482)

### Phylogeny-based

- 2020 - Darcy - [A phylogenetic model for the recruitment of species into microbial communities and application to studies of the human microbiome](https://doi.org/10.1038/s41396-020-0613-7)

### Causal Effect Estimation

- 2020 - **SparseMCMM** - Wang - linear log-contrast regression/Dirichlet regression - [Estimating and testing the microbial causal mediation effect with high-dimensional and compositional microbiome data](https://doi.org/10.1093/bioinformatics/btz565)

### Longitudinal Multi-omics Integration

- 2019 - Bodein - [A Generic Multivariate Framework for the Integration of Microbiome Longitudinal Studies With Other Data Types](https://www.frontiersin.org/articles/10.3389/fgene.2019.00963)
- 2019 - [**DIABLO**](https://github.com/singha53-zz/diablo) - Singh - Multi-block Sparse PLS-Discriminant Analysis - [DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays](https://doi.org/10.1093/bioinformatics/bty1054)
- 2022 - **MEFISTO** - Velten - [Identifying temporal and spatial patterns of variation from multi-modal data using MEFISTO](https://www.nature.com/articles/s41592-021-01343-9)
- 2023 - [**PALMO**](https://github.com/aifimmunology/PALMO) - Vasaikar - [A comprehensive platform for analyzing longitudinal multi-omics data](https://doi.org/10.1038/s41467-023-37432-w)

## Reviews / Evaluations / Opinion

- 2018 - Liang - [Dynamic modeling and network approaches for omics time course data](https://doi.org/10.1093/bib/bbx036)
- 2021 - Lv - [Causal Inference in Microbiome Medicine: Principles and Applications](https://doi.org/10.1016/j.tim.2021.03.015)
- 2022 - Corander - [Causal discovery for the microbiome](https://doi.org/10.1016/S2666-5247(22)00186-0)
- 2022 - Kodikara - [Statistical challenges in longitudinal microbiome data analysis](https://doi.org/10.1093/bib/bbac273)

## Longitudinal/Repeated Measures Application Papers
- 2021 - Duran-Pinedo - [Long-term dynamics of the human oral microbiome during clinical disease progression](https://doi.org/10.1186/s12915-021-01169-z)

## Meetings and workshops

