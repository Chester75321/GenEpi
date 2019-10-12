.. _example:

More Usage Examples
===================

For checking all the optional arguments, please use --help:

.. code-block:: none

   $ GenEpi --help

You will obtain the following argument list:

.. code-block:: none

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

Applying on Your Data
---------------------

.. code-block:: none

   $ GenEpi -g full_path_of_your_.GEN_file -p full_path_of_your_.CSV_file -o ./

For Quantatative Study
----------------------

GenEpi can support both case/control and quantitative studies, for quantitative studies please modify the parameter -m. You could download the `test data`_ and excute the following command.

.. code-block:: none

   $ GenEpi -g sample.gen -p sample_q.csv -o ./ -m r

.. _test data: https://github.com/Chester75321/GenEpi/raw/master/genepi/example

Changing the Genome Build
-------------------------

For changing the build of USCS genome browser, please modify the parameter -b:

.. code-block:: none

   $ GenEpi -g example -p example -o ./ --updatedb -b hg38

Threshold for Dimension Reduction  
---------------------------------

You could modify the threshold for Linkage Disequilibrium dimension reduction by following command:

.. code-block:: none

   $ GenEpi -g example -p example -o ./ --compressld -d 0.9 -r 0.9

Seld-defined Genome Region
--------------------------

Please prepare your seld-defined genome region in this `format <format\.html#Seld-defined Genome Regions>`_. Then, use the parameter -s for applying it on your data.

.. code-block:: none

   $ GenEpi -s full_path_of_your_genome_region_file -g full_path_of_your_.GEN_file -p full_path_of_your_.CSV_file -o ./

