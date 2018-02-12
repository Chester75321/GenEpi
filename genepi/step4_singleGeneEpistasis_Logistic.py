# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester Chang
"""

""""""""""""""""""""""""""""""
# importing libraries
""""""""""""""""""""""""""""""
import os
import sys
import numpy as np
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import chi2
from sklearn.preprocessing import PolynomialFeatures
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn.cross_validation import KFold
from sklearn import linear_model
from sklearn import grid_search
import sklearn.metrics as skMetric
import scipy.stats as stats
from psutil import virtual_memory
#import warnings
#warnings.filterwarnings("ignore")

""""""""""""""""""""""""""""""
# functions definition
""""""""""""""""""""""""""""""

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def SingleGeneEpistasisLogistic(str_inputFilePath_genotype, str_inputFilePath_phenotype, str_outputFilePath = ""):
    str_inputFilePath_genotype = "D:\\Phd\\Grade_05\\Alzheimer\\GenEpi\\genepi\\example\\snpSubsets\\APOE_11.gen"
    str_inputFilePath_phenotype = "D:\\Phd\\Grade_05\\Alzheimer\\GenEpi\\genepi\\example\\Sample.csv"
    str_outputFilePath = ""
    
    ### set path of output file
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFilePath_genotype)
    
    #-------------------------
    # load data
    #-------------------------
    ### count lines of input files
    int_num_genotype = sum(1 for line in open(str_inputFilePath_genotype))
    int_num_phenotype = sum(1 for line in open(str_inputFilePath_phenotype))
    
    ### get phenotype file
    list_phenotype = []
    with open(str_inputFilePath_phenotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_phenotype.append(line.strip().split(","))
    np_phenotype = np.array(list_phenotype)
    del list_phenotype
    
    ### get genotype file
    list_genotype_rsid = []
    ### declare a numpy array for one-hot-encoded genotype
    np_genotype = np.empty([int_num_phenotype, int_num_genotype * 3], dtype='int')
    with open(str_inputFilePath_genotype, 'r') as file_inputFile:
        idx_snp = 0
        for line in file_inputFile:
            list_thisSnp = line.strip().split(" ")
            list_genotype_rsid.append(list_thisSnp[1] + "_AA")
            list_genotype_rsid.append(list_thisSnp[1] + "_AB")
            list_genotype_rsid.append(list_thisSnp[1] + "_BB")
            for idx_subject in range(0, int_num_phenotype):
                list_allelType = [0, 0, 0]
                list_allelType[np.argmax(list_thisSnp[idx_subject * 3 + 5 : idx_subject * 3 + 8])] = 1
                np_genotype[idx_subject, idx_snp * 3 : idx_snp * 3 + 3] = list_allelType
            idx_snp = idx_snp + 1
    np_genotype_rsid = np.array(list_genotype_rsid)
    
    #-------------------------
    # data preprocessing
    #-------------------------    
    ### variance check (remove variance < 0.05)
    try:
        sk_variance = VarianceThreshold(threshold=(.95 * (1 - .95)))
        np_genotype = sk_variance.fit_transform(np_genotype)
        np_genotype_rsid = np_genotype_rsid[sk_variance.get_support()]
    except:
        return 0.0
    
    ### check memory leak
    mem = virtual_memory()
    print(sys.getsizeof(np_genotype) **2 / 2)
    print(mem.available)
    if sys.getsizeof(np_genotype) **2 / 2 > mem.available:
        return "MemErr"
    
    ### generate interaction term
    np_genotype_original = np_genotype
    np_genotype_original_rsid = np_genotype_rsid
    sklearn_poly = PolynomialFeatures(degree=2, interaction_only=True, include_bias=False)
    np_genotype = sklearn_poly.fit_transform(np_genotype)
    np_genotype_rsid = np.array(sklearn_poly.get_feature_names(np_genotype_rsid))
    ### variance check on ploynomial features
    try:
        sk_variance = VarianceThreshold(threshold=(.95 * (1 - .95)))
        np_genotype = sk_variance.fit_transform(np_genotype)
        np_genotype_rsid = np.array(np_genotype_rsid[sk_variance.get_support()])
    except:
        return 0.0
    #print("Generate Interaction Term: ")
    #print(str(np_genotype_rsid.shape[0]) + " with variance." + "\n")
    
    #-------------------------
    # feature selection
    #-------------------------
    """
    ### fisher exact test selection
    list_fisher = []
    for idx_feature in range(0, np_genotype.shape[1]):
        np_contingency = GenerateContingencyTable(np_genotype[:, idx_feature], np_phenotype[:, -1])
        oddsratio, pvalue = stats.fisher_exact(np_contingency)
        list_fisher.append(-np.log10(pvalue))
    np_fisher = np.array(list_fisher)
    np_selectedIdx = np.array([x > 5 for x in np_fisher])
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    print("Fisher Feature Selection: ")
    print(str(np_genotype_rsid.shape[0]) + " features been seleted.")
    """
    
    ### chi-square test selection
    np_chi2 = -np.log10(chi2(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
    np_selectedIdx = np.array([x > 2 for x in np_chi2])
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    #print("Chi2 Feature Selection: ")
    #print(str(np_genotype_rsid.shape[0]) + " features been seleted.")
    
    ### pre-modeling
    float_f1Score, np_weight = LogisticRegressionL1(np_genotype, np_phenotype[:, -1].astype(int))
    ### remove redundant polynomial features
    if float_f1Score > 0.0:
        try:
            for idx_original in range(0, np_genotype_original.shape[1]):
                list_selectedIdx = []
                for idx_polynomial in range(0, np_genotype.shape[1]):
                    list_selectedIdx.append(not(np.array_equal(np_genotype_original[:, idx_original], np_genotype[:, idx_polynomial])))
                np_selectedIdx = np.array(list_selectedIdx)
                np_genotype = np_genotype[:, np_selectedIdx]
                np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
            np_genotype = np.concatenate((np_genotype_original, np_genotype), axis=1)
            np_genotype_rsid = np.concatenate((np_genotype_original_rsid, np_genotype_rsid))
        except:
            return 0.0
    ### chi-square test selection
    np_chi2 = -np.log10(chi2(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
    np_selectedIdx = np.array([x > 2 for x in np_chi2])
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    #print("Remove Fake Interactive Features: ")
    #print(str(np_genotype_rsid.shape[0]) + " remained." + "\n")
    
    ### random logistic feature selection
    np_randWeight = np.array(RandomizedLogisticRegressionL1(np_genotype, np_phenotype[:, -1].astype(int)))
    np_selectedIdx = np.array([x >= 0.25 for x in np_randWeight])
    np_randWeight = np_randWeight[np_selectedIdx]
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    #print("Random Logistic Selection: ")
    #print(str(np_genotype_rsid.shape[0]) + " features been seleted." + "\n")
    
    #-------------------------
    # modeling
    #-------------------------
    ### final modeling
    float_f1Score, np_weight = LogisticRegressionL1(np_genotype, np_phenotype[:, -1].astype(int))
    #print("F1 Score: " + str(float_f1Score) + "\n")
    
    ### result analysis
    np_chi2 = -np.log10(chi2(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
    list_oddsRatio = []
    for idx_feature in range(0, np_genotype.shape[1]):
        np_contingency = GenerateContingencyTable(np_genotype[:, idx_feature], np_phenotype[:, -1])
        oddsratio, pvalue = stats.fisher_exact(np_contingency)
        list_oddsRatio.append(oddsratio)
    
    #-------------------------
    # outputting results
    #-------------------------
    ### statistics results
    with open(str_outputFilePath + str_inputFilePath_genotype.split("/")[-1].split("_")[0] + "_Result.csv", "w") as file_outputFile:
        file_outputFile.writelines("rsID,rand_weight,weight,chi2_p-value,odds_ratio" + "\n")
        for idx_feature in range(0, np_genotype_rsid.shape[0]):
            file_outputFile.writelines(str(np_genotype_rsid[idx_feature,]) + "," + str(np_randWeight[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_chi2[idx_feature,]) + "," + str(list_oddsRatio[idx_feature]) + "\n")
            
    ### output feature
    with open(str_outputFilePath + str_inputFilePath_genotype.split("/")[-1].split("_")[0] + "_Feature.csv", "w") as file_outputFile:
        file_outputFile.writelines(",".join(np_genotype_rsid) + "\n")
        for idx_subject in range(0, np_genotype.shape[0]):
            file_outputFile.writelines(",".join(np_genotype[idx_subject, :].astype(str)) + "\n")
    
    return float_f1Score