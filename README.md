# GenEpi
GenEpi is a package to uncover epistasis associated with phenotypes by a machine learning approach, developed by Yu-Chuan Chang at [c4Lab](http://bioinfo.bime.ntu.edu.tw/c4lab/) of National Taiwan University and Taiwan AI Labs

<img src="https://github.com/Chester75321/GenEpi/raw/master/GenEpi.png" width="90%" height="90%">

The architecture and modules of GenEpi.

## Getting Started
### Installation
```
$ pip install GenEpi
```
>**NOTE:** GenEpi is a memory-consuming package, which might cause memory errors when calculating the epistasis of a gene containing a large number of SNPs. We recommend that the memory for running GenEpi should be over 256 GB.

### Inputs
We provided test data [sample.gen](https://github.com/Chester75321/GenEpi/raw/master/genepi/example/sample.gen) and [sample.csv](https://github.com/Chester75321/GenEpi/raw/master/genepi/example/sample.csv) in [example folder](https://github.com/Chester75321/GenEpi/raw/master/genepi/example). Please see the following detail about input data.

**1\. Genotype Data:**

GenEpi takes the [Genotype File Format](http://www.cog-genomics.org/plink/1.9/formats#gen) (.GEN) used by Oxford statistical genetics tools, such as IMPUTE2 and SNPTEST as the input format for genotype data. If your files are in [PLINK format](http://www.cog-genomics.org/plink/1.9/formats) (.BED/.BIM/.FAM) or [1000 Genomes Project text Variant Call Format](http://www.cog-genomics.org/plink/1.9/formats#vcf) (.VCF), you could use [PLINK](http://www.cog-genomics.org/plink/1.9/) with the following command to convert the files to the .GEN file.

If your files are in the **.BED/.BIM/.FAM** format.
```
$ plink --bfile prefixOfTheFilename --recode oxford --out prefixOfTheFilename
```
If your file is in the **.VCF** format.
```
$ plink --vcf filename.vcf --recode oxford --out prefixOfTheFilename
```

**2\. Phenotype & Environmental Factor Data**

GenEpi takes the .CSV file without header line as the input format for phenotype and environmental factor data. The last column of the file will be considered as the phenotype data and the other columns will be considered as the environmental factor (covariates) data.
>**NOTE:** The sequential order of the phenotype data should be the same as that in the .GEN file.

## Usage Example
### Running a Quick Test
You will obtain all the outputs of GenEpi in current folder.
```
$ GenEpi -g example -p example -o ./
```

### Applying on Your Data
```
$ GenEpi -g full_path_of_your_.GEN_file -p full_path_of_your_.CSV_file -o ./
```

### Applying Seld-defined Genome Regions on Your Data
Prepare your genome regions in .TXT with the columns [chromosome, start, end, strand, geneSymbol], for eample:
```
1,10873,14409,+,DDX11L1
1,14361,30370,-,WASH7P
1,34610,37081,-,FAM138F
1,68090,70008,+,OR4F5
...
```

Then, use the parameter -s for applying it on your data
```
$ GenEpi -s full_path_of_your_genome_region_file -g full_path_of_your_.GEN_file -p full_path_of_your_.CSV_file -o ./
```

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

For changing the build of USCS genome browser, please modify the parameter -b:
```
$ GenEpi -g example -p example -o ./ --updatedb -b hg38
```

You could modify the threshold for Linkage Disequilibrium dimension reduction by following command:
```
$ GenEpi -g example -p example -o ./ --compressld -d 0.9 -r 0.9
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

After performing linkage disequilibrium (LD) dimension reduction, GenEpi will generate two files, a dimension-reduced .GEN file and a file containing LD blocks (.LDBlock file). Each row in the .LDBlock file indicates a LD block (see below for examples). The SNPs in front of colon signs are the representative SNPs of each LD block, and only these SNPs will be retained in the dimension-reduced .GEN file.
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

All of the within-gene epistasis selected by sinlge-gene models will be stored in the folder **singleGeneResult**, of which the format is the same as that in the **Result.csv** of cross-gene result. The performance of each single-gene model will be shown in **All_Logistic/Lasso_k-Fold.csv** in the same folder (see below for examples).

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
