# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester (Yu-Chuan Chang)
"""

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import os
import numpy as np
from .step4_singleGeneEpistasis_Logistic import LogisticRegressionL1
from .step4_singleGeneEpistasis_Lasso import LassoRegression
from .step5_crossGeneEpistasis_Logistic import RFClassifier
from .step5_crossGeneEpistasis_Logistic import RFClassifierModelPersistence
from .step5_crossGeneEpistasis_Lasso import RFRegressor
from .step5_crossGeneEpistasis_Lasso import RFRegressorModelPersistence

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def LoadDataForEnsemble(str_inputFileName_feature, str_inputFileName_phenotype):
    ### get all selected snp ids
    list_genotype_rsid = []
    with open(str_inputFileName_feature, "r") as file_inputFile:
        ### grep the header
        list_rsids = file_inputFile.readline().strip().split(",")
        for rsid in list_rsids:
            list_genotype_rsid.append(rsid)
    np_genotype_rsid = np.array(list_genotype_rsid)
    
    ### count lines of input files
    int_num_genotype = len(np_genotype_rsid)
    int_num_phenotype = sum(1 for line in open(str_inputFileName_phenotype))
    
    ### get phenotype file
    list_phenotype = []
    with open(str_inputFileName_phenotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_phenotype.append(line.strip().split(","))
    np_phenotype = np.array(list_phenotype)
    del list_phenotype
    
    if np_phenotype.shape[1] < 2:
        print("step6: Error no other factos exist.")
        return None, None
    
    ### get genotype file
    np_genotype = np.empty([int_num_phenotype, int_num_genotype], dtype='int')
    with open(str_inputFileName_feature, "r") as file_inputFile:
        ### skip header
        file_inputFile.readline()
        idx_phenotype = 0
        ### read feaure and write into np_genotype
        for line in file_inputFile:
            np_genotype[idx_phenotype, :len(list_rsids)] = np.array([float(x) for x in line.strip().split(",")], dtype='int')
            idx_phenotype = idx_phenotype + 1
    
    ### concatenate genotype and other factors
    np_genotype = np.concatenate((np_genotype, np_phenotype[:, :-1]), axis=1).astype(float)
    
    return np_genotype, np_phenotype

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def EnsembleWithCovariatesClassifier(str_inputFileName_feature, str_inputFileName_phenotype, str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 4):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_feature)
    
    #-------------------------
    # load data
    #-------------------------
    np_genotype, np_phenotype = LoadDataForEnsemble(str_inputFileName_feature, str_inputFileName_phenotype)
    if np_genotype is None and np_phenotype is None:
        return 0.0
    
    #-------------------------
    # build model
    #-------------------------
    float_f1Score, np_weight = LogisticRegressionL1(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs)
    
    ### use random forest classifier for ensemble
    float_f1Score_rf_te = RFClassifier(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs)
    ### random forest classifier model persistence
    float_f1Score_rf_tr = RFClassifierModelPersistence(np_genotype, np_phenotype[:, -1].astype(int), str_outputFilePath, str_outputFileName="RFClassifier_Covariates.pkl")
    
    print("step6: Ensemble with covariates. DONE! (Training score:" + "{0:.2f}".format(float_f1Score_rf_tr) + "; " + str(int_kOfKFold) + "-fold Test Score:" + "{0:.2f}".format(float_f1Score_rf_te) + ")")
    
    return float_f1Score

def EnsembleWithCovariatesRegressor(str_inputFileName_feature, str_inputFileName_phenotype, str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 4):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_feature)
    
    #-------------------------
    # load data
    #-------------------------
    np_genotype, np_phenotype = LoadDataForEnsemble(str_inputFileName_feature, str_inputFileName_phenotype)
    if np_genotype is None and np_phenotype is None:
        return 0.0
    
    #-------------------------
    # build model
    #-------------------------
    float_AVG_S_P, np_weight = LassoRegression(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs)
    
    ### use random forest classifier for ensemble
    float_AVG_S_P_rf_te = RFRegressor(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs)
    ### random forest classifier model persistence
    float_AVG_S_P_rf_tr = RFRegressorModelPersistence(np_genotype, np_phenotype[:, -1].astype(int), str_outputFilePath, str_outputFileName="RFRegressor_Covariates.pkl")
    
    print("step6: Ensemble with covariates. DONE! (Training score:" + "{0:.2f}".format(float_AVG_S_P_rf_tr) + "; " + str(int_kOfKFold) + "-fold Test Score:" + "{0:.2f}".format(float_AVG_S_P_rf_te) + ")")
    
    return float_AVG_S_P