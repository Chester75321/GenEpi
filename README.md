# GenEpi
GenEpi is a package to uncover epistasis associated with phenotypes by a machine learning approach, developed by Yu-Chuan Chang at [c4Lab](http://bioinfo.bime.ntu.edu.tw/c4lab/) of National Taiwan University.

<img src="https://github.com/Chester75321/GenEpi/raw/master/GenEpi.png" width="75%" height="75%">

The architecture and modules of GenEpi.

## Getting Started
### Installation
```
$ pip install GenEpi
```
>**NOTE:** GenEpi is a memory-consuming package, which might cause memory errors when calculating the epistasis of a gene containing a large number of SNPs. We recommend that the memory for running GenEpi should be over 256 GB.

### Inputs
**1\. Genotype Data:**

GenEpi takes the [Genotype File Format](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format_new.html) (.GEN) used by Oxford statistical genetics tools, such as IMPUTE2 and SNPTEST as the input format for genotype data. If your files are in [PLINK format](https://www.cog-genomics.org/plink/1.9/formats) (.BED/.BIM/.FAM) or [1000 Genomes Project text Variant Call Format](https://www.cog-genomics.org/plink/1.9/formats#vcf) (.VCF), you could use [PLINK](https://www.cog-genomics.org/plink/1.9/) with the following command to convert the files to the .GEN file.

If your files are in the **.BED/.BIM/.FAM** format.
```
$ plink --bfile prefixOfTheFilename --recode oxford --out prefixOfTheFilename
```
If your file is in the **.VCF** format.
```
$ plink --vcf filename.vcf --recode oxford --out prefixOfTheFilename
```

**2\. Phenotype & Environmental Factor Data**

GenEpi takes the .csv file without header line as the input format for phenotype and environmental factor data. The last column of the file will be considered as the phenotype data and the other columns will be considered as the environmental factor (covariates) data.
>**NOTE:** The sequential order of the phenotype data should be the same as that in the .GEN file.

## Usage Example
### Running a Test
We provided an [example script](https://github.com/Chester75321/GenEpi/tree/master/genepi/example/example.py) in [example folder](https://github.com/Chester75321/GenEpi/tree/master/genepi/example). Please use the following command for running a quick test.
```
$ python example.py
```

### Applying on Your Data
You may use this example script as a recipe and modify the input file names in Line 14 and 15 for running your data.
```python
str_inputFileName_genotype = "../sample.gen" # full path of the .GEN file.
str_inputFileName_phenotype = "../sample.csv" # full path of the .csv file.
```

### Options
For changing the build of USCS genome browser, please modify the parameter of the step one:
```python
genepi.DownloadUCSCDB(str_hgbuild="hg38") # for example: change to build hg38.
```
You could modify the threshold for Linkage Disequilibrium dimension reduction in the step two:
```python
#default: float_threshold_DPrime=0.9 and float_threshold_RSquare=0.9
genepi.EstimateLDBlock(str_inputFileName_genotype, float_threshold_DPrime=0.8, float_threshold_RSquare=0.8)
```

## Interpreting the Results
### The Main Table
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

### Other Details
**1\. Linkage Disequilibrium**

After performing linkage disequilibrium (LD) dimension reduction, GenEpi will generate two files, a dimension-reduced .GEN file and a file containing LD blocks (.LDBlock file). Each row in the .LDBlock file indicates a LD block (see below for example.). The SNPs in front of colon signs are the representative SNPs of each LD block, and only these SNPs will be retained in the dimension-reduced .GEN file.
```
rs429358:rs429358
rs7412:rs7412
rs117656888:rs117656888
rs1081105:rs1081105
rs1081106:rs1081106,rs191315680
```

**2\. Single-gene .GEN Files**

The subsets of the .GEN file for each gene will be stored in the folder **snpSubsets**.

**3\. Single-gene Results**

All of the within-gene epistasis selected by sinlge-gene models will be stored in the folder **singleGeneResult**, of which the format is the same as that in the **Result.csv** of cross-gene result. The performance of each single-gene model will be shown in **All_Logistic/Lasso_k-Fold.csv** in the same folder (see below for examples.).

| Gene Symbol | F1 Score |
|-------------|---------:|
| APOE        |   0.6109 |
| TOMM40      |   0.5900 |
| PVRL2       |   0.5745 |
| APOC1       |   0.5736 |

**4\. Model Persistance**

The final models of the step five and step six will be persisted in the folder **crossGeneResult** as **RFClassifier/Regressor.pkl** and **RFClassifier/Regressor_Covariates.pkl**, respectively. You could keep these models for future use without reconstructing them.

## Meta
Chester (Yu-Chuan Chang) - chester75321@gmail.com  
Distributed under the MIT license. See ``LICENSE`` for more information.  
[https://github.com/Chester75321/GenEpi/](https://github.com/Chester75321/GenEpi/)
