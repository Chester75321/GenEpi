# GenEpi
GenEpi is a package to uncover epistasis associated with phenotypes by a machine learning approach, developed by Yu-Chuan Chang at [c4Lab](http://bioinfo.bime.ntu.edu.tw/c4lab/) of National Taiwan University and Taiwan AI Labs

<img src="https://github.com/Chester75321/GenEpi/raw/master/GenEpi.png" width="90%" height="90%">

The architecture and modules of GenEpi.

## Introduction

GenEpi is designed to group SNPs by a set of loci in the gnome. For examples, a locus could be a gene. In other words, we use gene boundaries to group SNPs. A locus can be generalized to any particular regions in the genome, e.g. promoters, enhancers, etc. GenEpi first considers the genetic variants within a particular region as features in the first stage, because it is believed that SNPs within a functional region might have a higher chance to interact with each other and to influence molecular functions. 

GenEpi adopts two-element combinatorial encoding when producing features and models them by L1-regularized regression with stability selection In the first stage (STAGE 1) of GenEpi, the genotype features from each single gene will be combinatorically encoded and modeled independently by L1-regularized regression with stability selection. In this way, we can estimate the prediction performance of each gene and detect within-gene epistasis with a low false positive rate. In the second stage (STAGE 2), both of the individual SNP and the within-gene epistasis features selected by STAGE 1 are pooled together to generate cross-gene epistasis features, and modeled again by L1-regularized regression with stability selection as STAGE 1. Finally, the user can combine the selected genetic features with environmental factors such as clinical features to build the final prediction models.

## Standalone App
(Latest Update!) The standalone and installation free app - AppGenEpi (v.beta) is now released. Just download it then 1) chmod +x ./AppGenEpi_MacOS_beta and 2) sudo ./AppGenEpi_MacOS_beta for executing. Have fun~

| OS    |  Version |                                                Link                                                        |
|-------|:--------:|-----------------------------------------------------------------------------------------------------------:|
| MacOS | Catalina | [AppGenEpi_MacOS_beta](https://drive.google.com/file/d/1jOn3gOCpXW_OlpxNkpDzuRxt2_fIM8UM/view?usp=sharing) |
| Linux | CentOS 7 | [AppGenEpi_Linux_beta](https://drive.google.com/file/d/1hDpaoOgVXXYtHlBwx3QZ111OiyCqXSEE/view?usp=sharing) |


<img src="https://github.com/Chester75321/GenEpi/raw/master/AppGenEpi.png" width="100%" height="100%">

The snapshot of AppGenEpi.

## Citing

Please considering cite the following paper (currently avaliables as a pre-print in BioRxiv) if you use GenEpi in a scientific publication:

[1] Yu-Chuan Chang, June-Tai Wu, Ming-Yi Hong, Yi-An Tung, Ping-Han Hsieh, Sook Wah Yee, Kathleen M. Giacomini, Yen-Jen Oyang, and Chien-Yu Chen. "Genepi: Gene-Based Epistasis Discovery Using Machine Learning." bioRxiv  (2019): 421719.

## Quickstart
This section gets you started quickly. The completed GenEpi's documentation please find on [Welcome to GenEpi’s docs!](https://genepi.readthedocs.io/en/latest/)

### Installation
```
$ pip install GenEpi
```

>**NOTE:** GenEpi is a memory-consuming package, which might cause memory errors when calculating the epistasis of a gene containing a large number of SNPs. We recommend that the memory for running GenEpi should be over 256 GB.

### Running a quick test
Please use following command to run a quick test, you will obtain all the outputs of GenEpi in your current folder.
```
$ GenEpi -g example -p example -o ./
```

### Interpreting the main results table
GenEpi will automatically generate three folders (snpSubsets, singleGeneResult, crossGeneResult) beside your .GEN file. You could go to the folder **crossGeneResult** directly to obtain your main table for episatasis in **Result.csv**.

| RSID                        | -Log<sub>10</sub>(&chi;<sup>2</sup> p-value) | Odds Ratio | Genotype Frequency | Gene Symbol |
|-----------------------------|---------------------------------------------:|-----------:|-------------------:|-------------|
| rs157580_BB rs2238681_AA    |                                       8.4002 |     9.3952 |             0.1044 | TOMM40      |
| rs449647_AA rs769449_AB     |                                       8.0278 |     5.0877 |             0.2692 | APOE        |
| rs59007384_BB rs11668327_AA |                                       8.0158 |    12.0408 |             0.0824 | TOMM40      |
| rs283811_BB rs7254892_AA    |                                       8.0158 |    12.0408 |             0.0824 | PVRL2       |
| rs429358_AA                 |                                       5.7628 |     0.1743 |             0.5962 | APOE        |
| rs73052335_AA rs429358_AA   |                                       5.6548 |     0.1867 |             0.5714 | APOC1\*APOE |

>The first column lists each feature by its RSID and the genotype (denoted as RSID_genotype), the pairwise epistatis features are represented using two SNPs. The last column describes the genes where the SNPs are located according to the genomic coordinates. We used a star sign to denote the epistasis between genes. The p-values of the &chi;<sup>2</sup> test (the quantitative task will use student t-test) are also included. The odds ratio significantly away from 1 also indicates whether the features are potential causal or protective genotypes. Since low genotype frequency may cause unreliable odds ratios, we also listed this information in the table.

### Options
For checking all the optional arguments, please use --help:
```
$ GenEpi --help
```

You will obtain the following argument list:
```
usage: GenEpi [-h] -g G -p P [-s S] [-o O] [-m {c,r}] [-k K] [-t T]
              [--updatedb] [-b {hg19,hg38}] [--compressld] [-d D] [-r R]

optional arguments:
  -h, --help      show this help message and exit
  -g G            filename of the input .gen file
  -p P            filename of the input phenotype
  -s S            self-defined genome regions
  -o O            output file path
  -m {c,r}        choose model type: c for classification; r for regression
  -k K            k of k-fold cross validation
  -t T            number of threads

update UCSC database:
  --updatedb      enable this function
  -b {hg19,hg38}  human genome build

compress data by LD block:
  --compressld    enable this function
  -d D            threshold for compression: D prime
  -r R            threshold for compression: R square
```

## Meta
Chester (Yu-Chuan Chang) - chester75321@gmail.com  
Distributed under the MIT license. See ``LICENSE`` for more information.  
[https://github.com/Chester75321/GenEpi/](https://github.com/Chester75321/GenEpi/)
