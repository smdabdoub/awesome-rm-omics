# Awesome Repeated Measures "omics" Analysis

A [community-maintained](https://github.com/smdabdoub/awesome-rm-omics/graphs/contributors) list of software packages for analysis of longitudinal/repeated measures studies utilizing high-dimensional "omics" data.

While many of the packages here are marketed for a single type of "omics" data (metagenomics, metabolome, etc.). Many can be utilized across the -omics universe.

The common thread among the methods listed here is that the same samples are measured across multiple time-points and may be measured with two or more modalities (e.g.: microbiome, metabolome, and miRNA). The data can be described as multiple matrices/tables with the same number of samples (not a strict requirement for all methods) and varying number of features.

The repo is in the style of Mike Love's [awesome-multi-omics](https://github.com/mikelove/awesome-multi-omics) repository for multi-omics analysis methods (which, in turn, is in the style of Sean Davis' [awesome-single-cell](https://github.com/seandavi/awesome-single-cell) repo for single-cell analysis methods). Note that there is a Multi-omics section for software handling repeated measures studies with multiple measurement modalities, so there will be some overlap with [awesome-multi-omics](https://github.com/mikelove/awesome-multi-omics).

[Contributions are welcome and encouraged](https://github.com/smdabdoub/awesome-rm-omics/blob/master/CONTRIBUTING.md). I can't read all the papers by myself.

For brevity, entries for software and scientific articles only list the last name of the first author.

**NOTE**: Entries with ** before the year are non-peer reviewed preprints.

## Software packages and methods

### Spline Modeling

- 2015 - Straube - [A Linear Mixed Model Spline Framework for Analysing Time Course ‘Omics’ Data](https://doi.org/10.1371/journal.pone.0134540)
- 2017 - [**metaDprof**](https://cals.arizona.edu/~anling/software/metaDprof.htm) - Luo - [An informative approach on differential abundance analysis for time-course metagenomic sequencing data](https://doi.org/10.1093/bioinformatics/btw828)
- 2018 - [**SplinectomeR**](https://rrshieldscutler.github.io/splinectomeR/) - Shields-Cutler - [SplinectomeR Enables Group Comparisons in Longitudinal Microbiome Studies](https://doi.org/10.3389/fmicb.2018.00785)
- 2018 - [**MetaLonDA**](https://github.com/aametwally/MetaLonDA) - Metwally - [MetaLonDA: a flexible R package for identifying time intervals of differentially abundant features in metagenomic longitudinal studies](https://doi.org/10.1186/s40168-018-0402-y) - [vignette](https://cran.r-project.org/web/packages/MetaLonDA/vignettes/MetaLonDA.html)

### Linear/Mixed Modeling

- 2016 - [**ZIBR**](https://github.com/chvlyl/ZIBR) - Chen - zero-inflated Beta regression - [A two-part mixed-effects model for analyzing longitudinal microbiome compositional data](https://doi.org/10.1093/bioinformatics/btw308)
- 2018 - [**LassoGLMMforMicrobiomes**](https://github.com/ghedin-lab/LassoGLMMforMicrobiomes) - Tipton - lasso-penalized generalized linear mixed model (LassoGLMM) with variable selection - [Measuring associations between the microbiota and repeated measures of continuous clinical variables using a lasso-penalized generalized linear mixed model](https://doi.org/10.1186/s13040-018-0173-9)
- 2018 - [**MALLARD**](https://github.com/LAD-LAB/MALLARD-Paper-Code) - Silverman - [Dynamic linear models guide design and analysis of microbiota studies within artificial human guts](https://doi.org/10.1186/s40168-018-0584-3)
- 2020 - [**NBZIMM**](https://github.com/nyiuab/NBZIMM) - Zhang - [NBZIMM: negative binomial and zero-inflated mixed models, with application to microbiome/metagenomics data analysis](https://doi.org/10.1186/s12859-020-03803-z) - [vignette](https://abbyyan3.github.io/NBZIMM-tutorial/)
- 2020 - Zhang - [Fast zero-inflated negative binomial mixed modeling approach for analyzing longitudinal metagenomics data](https://doi.org/10.1093/bioinformatics/btz973)
- 2021 - [**MaAsLin 2**](https://huttenhower.sph.harvard.edu/maaslin/) - Mallick - [Multivariable association discovery in population-scale meta-omics studies](https://doi.org/10.1371/journal.pcbi.1009442) - [vignette](https://github.com/biobakery/biobakery/wiki/maaslin2)
- 2022 - [**MiRKAT-MC**](https://github.com/Zhiwen-Owen-Jiang/MiRKATMC) - Jiang - [MiRKAT-MC: A Distance-Based Microbiome Kernel Association Test With Multi-Categorical Outcomes](https://doi.org/10.3389/fgene.2022.841764) - [vignette](https://cran.r-project.org/web/packages/MiRKAT/vignettes/MiRKAT_Vignette.html)
- 2023 - [**coda4microbiome**](https://malucalle.github.io/coda4microbiome/) - Calle - [coda4microbiome: compositional data analysis for microbiome cross-sectional and longitudinal studies](https://doi.org/10.1186/s12859-023-05205-3)

### Multi-level Modeling

- 2012 - [**mixOmics**](http://mixomics.org/) - Liquet - [A novel approach for biomarker selection and the integration of repeated measures experiments from two assays](https://doi.org/10.1186/1471-2105-13-325)

### Bayesian Modeling

- 2015 - [**BioMiCo**](https://sourceforge.net/projects/biomico/) - Shafiei - [BioMiCo: a supervised Bayesian model for inference of microbial community structure](https://doi.org/10.1186/s40168-015-0073-x)

### Stochastic/Probabilistic Modeling

- 2018 - [**TGP-CODA**](https://github.com/tare/GPMicrobiome) - Äijö - [Temporal probabilistic modeling of bacterial compositions derived from 16S rRNA sequencing](https://doi.org/10.1093/bioinformatics/btx549)
- 2020 - [**LUMINATE**](https://github.com/tyjo/luminate) - Joseph - [Efficient and Accurate Inference of Mixed Microbial Population Trajectories from Longitudinal Count Data](https://doi.org/10.1016/j.cels.2020.05.006)

### Network Analysis

- 2011 - [**eLSA**](https://bitbucket.org/charade/elsa/wiki/Home) - Xia - [Extended local similarity analysis (eLSA) of microbial community and other time series data with replicates](https://doi.org/10.1186%2F1752-0509-5-S2-S15)
- 2018 - Coenen - [Limitations of Correlation-Based Inference in Complex Virus-Microbe Communities](https://doi.org/10.1128/msystems.00084-18)
- 2019 - Ai - bivariate autoregressive model with Granger causality - [Constructing the Microbial Association Network from large-scale time series data using Granger causality](https://doi.org/10.3390%2Fgenes10030216)

### Machine Learning/Artificial Neural Networks (Deep Learning) Methods

- 2021 - [**phyLoSTM**](https://github.com/divya031090/phyLoSTM) - Sharma - combined modeling using CNN for feature extraction and LSTM for temporal dependency analysis - [phyLoSTM: a novel deep learning model on disease prediction from longitudinal microbiome data](https://doi.org/10.1093/bioinformatics/btab482)
- 2023 - Hu - [A review on longitudinal data analysis with random forest](https://doi.org/10.1093/bib/bbad002)
- 2023 - [**EMBED**](https://github.com/mayar-shahin/EMBED) - Shahin - [EMBED: Essential MicroBiomE Dynamics, a dimensionality reduction approach for longitudinal microbiome studies](https://doi.org/10.1038/s41540-023-00285-6)

### Ecology Literature

- 2011 - Lindenmayer - [Cross-sectional vs. longitudinal research: a case study of trees with hollows and marsupials in Australian forests](https://doi.org/10.1890/11-0279.1)
- 2019 - [**pldist**](https://github.com/aplantin/pldist) - Plantinga - [pldist: ecological dissimilarities for paired and longitudinal microbiome association analysis](https://doi.org/10.1093/bioinformatics/btz120)
- 2019 - McGeoch - [Measuring continuous compositional change using decline and decay in zeta diversity](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2832)
- 2019 - Buckley - [Multivariate methods for testing hypotheses of temporal community dynamics](https://doi.org/10.1101/362822) - [paper 2](https://doi.org/10.7717/peerj.11250)
- 2020 - Yang - [Toward a more temporally explicit framework for community ecology](https://doi.org/10.1111/1440-1703.12099)

### Phylogeny-based

- 2020 - Darcy - [A phylogenetic model for the recruitment of species into microbial communities and application to studies of the human microbiome](https://doi.org/10.1038/s41396-020-0613-7)

### Causal Effect Estimation
Not technically longitudinal, but heading in that direction and good to know about.

- 2020 - [**SparseMCMM**](https://github.com/chanw0/SparseMCMM) - Wang - linear log-contrast regression/Dirichlet regression - [Estimating and testing the microbial causal mediation effect with high-dimensional and compositional microbiome data](https://doi.org/10.1093/bioinformatics/btz565)
- 2022 - [**CMMB**](https://github.com/mbsohn/cmmb) - Sohn - [A compositional mediation model for a binary outcome: Application to microbiome studies](https://doi.org/10.1093/bioinformatics/btab605)
- 2023 - [**SparseMCMM_HD**](https://github.com/chanw0/SparseMCMM) - - [A microbial causal mediation analytic tool for health disparity and applications in body mass index](https://doi.org/10.21203/rs.3.rs-2463503/v1)

### Pathway Analysis

- 2012 - **DIA** - Bionaz - [A Novel Dynamic Impact Approach (DIA) for Functional Analysis of Time-Course Omics Studies](https://doi.org/10.1371/journal.pone.0032455)
- 2016 - Du - [A decision analysis model for KEGG pathway analysis](https://doi.org/10.1186/s12859-016-1285-1)

### Longitudinal Multi-omics Integration

- 2019 - [**timeOmics**](https://github.com/abodein/timeOmics) - Bodein - [A Generic Multivariate Framework for the Integration of Microbiome Longitudinal Studies With Other Data Types](https://www.frontiersin.org/articles/10.3389/fgene.2019.00963) [vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/timeOmics/inst/doc/vignette.html)
- 2019 - [**DIABLO**](https://github.com/singha53-zz/diablo) - Singh - Multi-block Sparse PLS-Discriminant Analysis - [DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays](https://doi.org/10.1093/bioinformatics/bty1054)
- 2020 - [**PyIOmica**](https://pypi.python.org/pypi/pyiomica) - Domanskyi - [PyIOmica: longitudinal omics analysis and trend identification](https://doi.org/10.1093/bioinformatics/btz896)
- 2022 - [**MEFISTO**](https://biofam.github.io/MOFA2/MEFISTO.html) - Velten - [Identifying temporal and spatial patterns of variation from multi-modal data using MEFISTO](https://www.nature.com/articles/s41592-021-01343-9) - [vignette](https://biofam.github.io/MOFA2/tutorials.html)
- 2023 - [**PALMO**](https://github.com/aifimmunology/PALMO) - Vasaikar - [A comprehensive platform for analyzing longitudinal multi-omics data](https://doi.org/10.1038/s41467-023-37432-w)

### Data Simulation

- ** 2022 - [**MTIST**](https://github.com/jsevo/mtist) - Hussey - [The MTIST platform: a microbiome time series inference standardized test simulation, dataset, and scoring systems](https://doi.org/10.1101/2022.10.18.512783)
- ** 2024 - [**SimMiL**](https://github.com/nweaver111/SimMiL) - Weaver & Hendricks - [SimMiL: Simulating Microbiome Longitudinal Data](https://doi.org/10.1101/2024.03.18.585571)

## Reviews / Evaluations / Opinion

- 2018 - Liang - [Dynamic modeling and network approaches for omics time course data](https://doi.org/10.1093/bib/bbx036)
- 2021 - Lv - [Causal Inference in Microbiome Medicine: Principles and Applications](https://doi.org/10.1016/j.tim.2021.03.015)
- 2022 - Corander - [Causal discovery for the microbiome](https://doi.org/10.1016/S2666-5247(22)00186-0)
- 2022 - Kodikara - [Statistical challenges in longitudinal microbiome data analysis](https://doi.org/10.1093/bib/bbac273)
- 2023 - Grieneisen - [How longitudinal data can contribute to our understanding of host genetic effects on the gut microbiome](https://doi.org/10.1080%2F19490976.2023.2178797)
- 2024 - Lyu - [Methodological Considerations in Longitudinal Analyses of Microbiome Data: A Comprehensive Review](https://doi.org/10.3390/genes15010051)

## Longitudinal/Repeated Measures Application Papers
- 2021 - Duran-Pinedo - [Long-term dynamics of the human oral microbiome during clinical disease progression](https://doi.org/10.1186/s12915-021-01169-z)
- 2022 - Armoni - [Temporal Alignment of Longitudinal Microbiome Data](https://doi.org/10.3389/fmicb.2022.909313)
- 2022 - Murillo - [Assessing the drivers of gut microbiome composition in wild redfronted lemurs via longitudinal metacommunity analysis](https://doi.org/10.1038/s41598-022-25733-x)

## Meetings and workshops

