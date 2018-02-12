# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester Chang
"""

import os
import genepi

### step1_downloadUCSCDB
genepi.DownloadUCSCDB(str_hgbuild="hg19")

### step2_estimateLD
genepi.EstimateLDBlock(os.path.dirname(os.path.abspath(__file__)) + "/sample.gen", float_threshold_DPrime=0.9, float_threshold_RSquare=0.9)

### step3_splitByGene
genepi.SplitByGene(os.path.dirname(os.path.abspath(__file__)) + "/sample_LDReduced.gen")

str_inputFilePath_genotype = "D:\\Phd\\Grade_05\\Alzheimer\\GenEpi\\genepi\\example\\snpSubsets\\APOE_11.gen"
str_inputFilePath_phenotype = "D:\\Phd\\Grade_05\\Alzheimer\\GenEpi\\genepi\\example\\Sample.csv"
genepi.SingleGeneEpistasisLogistic(str_inputFilePath_genotype,str_inputFilePath_phenotype)
    
