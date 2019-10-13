.. _quickstart:

Quickstart
==========

This section gets you started quickly, the I/O described in `I/O File Fomats <format\.html>`_, more usage examples please find in `More Usage Examples <example\.html>`_, discussing each of the sub-modules introduced in `How it Work <workflow\.html>`_.

Running a Quick Test
--------------------

Please use following command to run a quick test, you will obtain all the outputs of GenEpi in your current folder.

.. code-block:: none

   $ GenEpi -g example -p example -o ./

The progress will print on console:

.. code-block:: none

   step1: Down load UCSC Database. DONE!
   step2: Estimate LD. DONE!
   Warning of step3: .gen file should be sorted by chromosome and position
   step3: Split by gene. DONE!
   step4: Detect single gene epistasis. DONE!
   step5: Detect cross gene epistasis. DONE! (Training score:0.63; 2-fold Test Score:0.61)
   step6: Ensemble with covariates. DONE! (Training score:0.63; 2-fold Test Score:0.60)

GenEpi will automatically generate three folders (snpSubsets, singleGeneResult, crossGeneResult) in output path (arg: -o). The following tree structure is the contents of the output folder. You could go to the folder **crossGeneResult** directly to obtain your main result table for episatasis in **Result.csv**.

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


Interpreting the Main Result Table
----------------------------------

Here is the contents of Result.csv, which mean the episatasis seleted by GenEpi.

=========================== ================== ========== ================== ===========
RSID                        -Log10(χ2 p-value) Odds Ratio Genotype Frequency Gene Symbol
=========================== ================== ========== ================== ===========
rs157580_BB rs2238681_AA    8.4002             9.3952     0.1044             TOMM40
rs449647_AA rs769449_AB     8.0278             5.0877     0.2692             APOE
rs59007384_BB rs11668327_AA 8.0158             12.0408    0.0824             TOMM40
rs283811_BB rs7254892_AA    8.0158             12.0408    0.0824             PVRL2
rs429358_AA                 5.7628             0.1743     0.5962             APOE
rs73052335_AA rs429358_AA   5.6548             0.1867     0.5714             APOC1*APOE
=========================== ================== ========== ================== ===========

We listed the statistical significance of the selected genetic features in Result.csv. The first column lists each feature by its RSID and the genotype (denoted as RSID_genotype), the pairwise epistatis features are represented using two SNPs. The last column describes the genes where the SNPs are located according to the genomic coordinates. We used a star sign to denote the epistasis between genes. The p-values of the χ2 test (the quantitative task will use student t-test) are also included. The odds ratio significantly away from 1 also indicates whether the features are potential causal or protective genotypes. Since low genotype frequency may cause unreliable odds ratios, we also listed this information in the table.