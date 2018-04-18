# GenEpi
GenEpi is a package to uncover epistasis associated with phenotypes by machine learning approach, developed by Yu-Chuan Chang at [c4Lab](http://bioinfo.bime.ntu.edu.tw/c4lab/) of National Taiwan University.

![](https://github.com/Chester75321/GenEpi/blob/master/GenEpi.png)

## Getting Started
### Installation
```
$ pip install GenEpi
```
>**NOTE:** GenEpi is a memory-consuming package, which might cause memory errors when calculating the epistasis of a gene containing a large number of SNPs. We recommend that the memory for running GenEpi should be over 256 GB.

### Inputs
1. Genotype Data
GenEpi takes [Genotype File Format](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format_new.html) (.GEN) used by Oxford statistical genetics tools, such as IMPUTE2 and SNPTEST as input format for genotype data. If your files is in PLINK text format (.PED and .MAP), you could use [GTOOL](http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html) with following command to convert .PED files to .GEN file.
>```
>$ gtool -P --ped example.ped --map example.map --og out.gen --os out.sample
>```
2. Phenotype & Environmental Factor Data
GenEpi takes .csv file without header line as input format for phenotype and environmental factor data. The last column of the file would be considered as phenotype data and the others would be considered as covariate (environmental factor data).

>**NOTE:** The order of the phenotype data should be same as .GEN file.


