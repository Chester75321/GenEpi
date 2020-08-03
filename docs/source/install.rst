.. _install:

Installation 
============

How to Install
--------------

GenEpi supports Python 3.7 and up. Use pip to install GenEpi and its dependencies.

.. code-block:: none

   $ pip install GenEpi

Check that you installed the GenEpi sucessfully.

.. code-block:: none

   $ GenEpi --help

After executed previous command on console, you will see:

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

Dependencies
------------

Here is the dependency list for running GenEpi. pip takes care of these dependencies automatically when you install GenEpi.

   - numpy >= 1.13.0 
   - psutil >= 4.3.0
   - pymysql >= 0.8.0
   - scipy >= 0.19.0
   - scikit-learn == 0.21.2

System Requirements
-------------------

For running a quick test, you could install GenEpi on any laptop e.g. a MacBook. When applying GenEpi on a real whole genome-wide dataset, here are recommended system requirements:

:Processor: 2.3 GHz Intel XEON® E5-2673 v4 * 32
:RAM: 256 GiB
:Storage: 500 GiB

These requirements are refer to the specification of Microsoft Azure E32 v3.

.. note::

   GenEpi is a memory-consuming package, which might cause memory errors when calculating the epistasis of a gene containing a large number of SNPs. We recommend that the memory for running GenEpi should be over 256 GB.