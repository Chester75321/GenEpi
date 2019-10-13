.. _format:

I/O File Fomats
===============

We provided test data `sample.gen`_ and `sample.csv`_ in `example folder`_. After running a `quick test <quickstart\.html#Running a Quick Test>`_, GenEpi will automatically copy these test data to output path and generate all of the output folders and files as following tree structure. Please see the following detail about the format of these input and output data.

.. code-block:: none

   ./
   ├── GenEpi_Log_DATE-TIME.txt
   ├── crossGeneResult
   │   ├── Classifier.pkl
   │   ├── Classifier_Covariates.pkl
   │   ├── Feature.csv
   │   └── Result.csv
   ├── sample.LDBlock
   ├── sample.csv
   ├── sample.gen
   ├── sample_LDReduced.gen
   ├── singleGeneResult
   │   ├── All_Logistic_k2.csv
   │   ├── APOC1_Feature.csv
   │   ├── APOC1_Result.csv
   │   ├── APOE_Feature.csv
   │   ├── APOE_Result.csv
   │   ├── PVRL2_Feature.csv
   │   ├── PVRL2_Result.csv
   │   ├── TOMM40_Feature.csv
   │   └── TOMM40_Result.csv
   └── snpSubsets
       ├── APOC1_23.gen
       ├── APOE_11.gen
       ├── PVRL2_48.gen
       └── TOMM40_67.gen


.. _sample.gen: https://github.com/Chester75321/GenEpi/raw/master/genepi/example/sample.gen
.. _sample.csv: https://github.com/Chester75321/GenEpi/raw/master/genepi/example/sample.csv
.. _example folder: https://github.com/Chester75321/GenEpi/raw/master/genepi/example

Input: Genotype Data
--------------------

`Sample.gen`_ is an example of genotype data. GenEpi takes the `Genotype File Format`_ (.GEN) used by Oxford statistical genetics tools, such as IMPUTE2 and SNPTEST as the input format for genotype data. If your files are in `PLINK format`_ (.BED/.BIM/.FAM) or `1000 Genomes Project text Variant Call Format`_ (.VCF), you could use `PLINK`_ with the following command to convert them to the .GEN file.

If your files are in the .BED/.BIM/.FAM format.

.. code-block:: none

   $ plink --bfile prefixOfTheFilename --recode oxford --out prefixOfTheFilename

If your file is in the .VCF format.

.. code-block:: none

   $ plink --vcf filename.vcf --recode oxford --out prefixOfTheFilename

.. _Sample.gen: https://github.com/Chester75321/GenEpi/raw/master/genepi/example/sample.gen
.. _Genotype File Format: http://www.cog-genomics.org/plink/1.9/formats#gen
.. _PLINK format: http://www.cog-genomics.org/plink/1.9/formats
.. _1000 Genomes Project text Variant Call Format: http://www.cog-genomics.org/plink/1.9/formats#vcf
.. _PLINK: http://www.cog-genomics.org/plink/1.9/

Input: Phenotype Data
---------------------

`Sample.gen`_ is an example of phenotype data with environmental factor (e.g. the age of each sample). GenEpi takes the common .CSV file without header line as the input format for phenotype and environmental factor data. The last column of the file will be considered as the phenotype data (e.g. 1 or 0, which indicate case/control, respectively) and the other columns will be considered as the environmental factor data.

.. warning::

   The sequential order of the phenotype data should be the same as that in the .GEN file.

Output: LDBlock File
--------------------

GenEpi has a moudle, which can reduce the dimension of the input feature by estimating the linkage disequilibrium (LD). After performing dimension reduction, GenEpi will generate two files, a dimension-reduced .GEN file and a file containing LD blocks (.LDBlock file). Each row in the .LDBlock file indicates a LD block (see below for examples). The SNPs in front of colon signs are the representative SNPs of each LD block, and only these SNPs will be retained in the dimension-reduced .GEN file.

.. code-block:: none

   rs429358:rs429358
   rs7412:rs7412
   rs117656888:rs117656888
   rs1081105:rs1081105
   rs1081106:rs1081106,rs191315680

Output: SnpSubsets Folder
-------------------------

Since GenEpi is a gene-based epistasis discovering method, the input genotype will first be splited into group of each gene. The subsets of the .GEN file for each gene will be stored in the folder snpSubsets. The naming rule of the filename is GeneSymbol_NumberOfVariantsOnGene.gen

Output: SingleGeneResult Folder
-------------------------------

In first stage of GenEpi, all the .GEN file for each gene will be modeled gene by gene. Every models will output two kins of data the GeneSymbol_Result.csv and GeneSymbol_Feature.csv. The format of GeneSymbol_Result.csv please refer to the section `Interpreting the Main Result Table <quickstart\.html#Interpreting the Main Result Table>`_. The only difference is the GeneSymbol_Result.csv is a result table for single gene. Moreover, GeneSymbol_Feature.csv are the raw features corresponds to the episatsis in GeneSymbol_Result.csv. Beside these two types of files, there is a file named All_Logistic/Lasso.csv, which contains all the prediction scores of each gene, please see as following example.

.. code-block:: none

   $ head All_Logistic_k2.csv

   GeneSymbol,F1Score
   PVRL2,0.5745454545454547
   APOC1,0.5681818181818181
   TOMM40,0.602510460251046
   APOE,0.592

Output: CrossGeneResult Folder
------------------------------

The results of the second stage of GenEpi - cross gene modeling will be generated in this folder. The formats are as same as the description in previous section `SingleGeneResult Folder <#SingleGeneResult Folder>`_. Moreover, the final models will be persisted in this folder as Classifier/Regressor.pkl and Classifier/Regressor_Covariates.pkl, respectively. You could keep these models for future use without reconstructing them.

Output: Porcess Log
-------------------

The performance of genetic feature only and genetic + environmental factor models will be logged into GenEpi_Log_DATE-TIME.txt. Other process information such as the setting of arguments, the time cost will also be recorded.

.. code-block:: none

   start analysis at: 20191009-17:54:45

   Arguments in effect:
           -g (input genotype filename): ./sample.gen
           -p (input phenotype filename): ./sample.csv
           -s (self-defined genome regions): None
           -o (output filepath): ./

           -m (model type): Classification -k (k-fold cross validation): 2
           -t (number of threads): 4

           --updatedb (enable function of update UCSC database): True
           -b (human genome build): hg19

           --compressld (enable function of LD data compression): True
           -d (D prime threshold): 0.9
           -r (R square threshold): 0.9

   Number of variants: 223
   Number of samples: 364
   Overall genetic feature performance (F1 score)
   Training: 0.6307053941908715
   Testing (2-fold CV): 0.6134453781512604

   Ensemble with co-variate performance (F1 score)
   Training: 0.632
   Testing (2-fold CV): 0.6016260162601627

   end analysis at: 20191009-17:54:58

Seld-defined Genome Regions
---------------------------

GenEpi supports seld-defined genome regions for first stage to subset the data. Please prepare your genome regions in .TXT with the columns [chromosome, start, end, strand, geneSymbol], for eample:

.. code-block:: none

   1,10873,14409,+,DDX11L1
   1,14361,30370,-,WASH7P
   1,34610,37081,-,FAM138F
   1,68090,70008,+,OR4F5
   ...