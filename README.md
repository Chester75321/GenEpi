# GenEpi
GenEpi is a package to uncover epistasis associated with phenotypes by machine learning approach, developed by Yu-Chuan Chang at [c4Lab](http://bioinfo.bime.ntu.edu.tw/c4lab/) of National Taiwan University.

<img src="https://github.com/Chester75321/GenEpi/raw/master/GenEpi.png" width="75%" height="75%">

The architecture and modules of GenEpi.

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

## Usage example
### Running a test
We provided an [example script](https://github.com/Chester75321/GenEpi/tree/master/genepi/example/example.py) in [example folder](https://github.com/Chester75321/GenEpi/tree/master/genepi/example). Please use following command for running a quick test.
```
$ python example.py
```

### Applying on your data
You may use this example script as a recipe and modify the input file names in Line 14 and 15 for running your data.
```python
str_inputFileName_genotype = "../sample.gen" # full path of the .GEN file.
str_inputFileName_phenotype = "../sample.csv" # full path of the .csv file.
```

### Options
For changing the build of USCS genome browser, please modify parameter of step one:
```python
genepi.DownloadUCSCDB(str_hgbuild="hg38") # for example: change to build hg38.
```
You could modify the threshold for Linkage Disequilibrium dimension reduction in step two:
```python
#default float_threshold_DPrime=0.9 and float_threshold_RSquare=0.9
genepi.EstimateLDBlock(str_inputFileName_genotype, float_threshold_DPrime=0.8, float_threshold_RSquare=0.8)
```

## Meta
Chester (Yu-Chuan Chang) - chester75321@gmail.com  
Distributed under the MIT license. See ``LICENSE`` for more information.  
[https://github.com/Chester75321/GenEpi/](https://github.com/Chester75321/GenEpi/)
