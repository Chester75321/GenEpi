# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester (Yu-Chuan Chang)
"""

import os
import genepi


def main():
    # give file names of "Genotype" and "Phenotype" here !
    str_inputFileName_genotype = os.path.dirname(os.path.abspath(__file__)) + "/sample.gen"
    str_inputFileName_phenotype = os.path.dirname(os.path.abspath(__file__)) + "/sample.csv"
    
    ### step1_downloadUCSCDB
    genepi.DownloadUCSCDB(str_hgbuild="hg19")
    
    ### step2_estimateLD
    genepi.EstimateLDBlock(str_inputFileName_genotype, float_threshold_DPrime=0.9, float_threshold_RSquare=0.9)
    
    ### step3_splitByGene
    genepi.SplitByGene(str_inputFileName_genotype.replace(".gen", "_LDReduced.gen"))
    
    ### step4_singleGeneEpistasis_Logistic (for case/control trial)
    genepi.BatchSingleGeneEpistasisLogistic(os.path.dirname(str_inputFileName_genotype) + "/snpSubsets/", str_inputFileName_phenotype, int_nJobs=1)
    ### step4_singleGeneEpistasis_Lasso (for quantitative trial)
    #genepi.BatchSingleGeneEpistasisLasso(os.path.dirname(str_inputFileName_genotype) + "/snpSubsets/", str_inputFileName_phenotype, int_nJobs=1)
    
    ### step5_crossGeneEpistasis_Logistic (for case/control trial)
    genepi.CrossGeneEpistasisLogistic(os.path.dirname(str_inputFileName_genotype) + "/singleGeneResult/", str_inputFileName_phenotype, int_nJobs=1)
    ### step5_crossGeneEpistasis_Lasso (for quantitative trial)
    #genepi.CrossGeneEpistasisLasso(os.path.dirname(str_inputFileName_genotype) + "/singleGeneResult/", str_inputFileName_phenotype, int_nJobs=1)
    
    ### step6_ensembleWithCovariates (for case/control trial)
    genepi.EnsembleWithCovariatesClassifier(os.path.dirname(str_inputFileName_genotype) + "/crossGeneResult/Feature.csv", str_inputFileName_phenotype, int_nJobs=1)
    ### step6_ensembleWithCovariates (for quantitative trial)
    #genepi.EnsembleWithCovariatesRegressor(os.path.dirname(str_inputFileName_genotype) + "/crossGeneResult/Feature.csv", str_inputFileName_phenotype, int_nJobs=1)
    
if __name__ == '__main__':
    main()
