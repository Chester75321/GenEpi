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
import sys
import numpy as np
from psutil import virtual_memory
from sklearn.preprocessing import PolynomialFeatures
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import f_regression
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn import linear_model
from sklearn.cross_validation import KFold
from sklearn import grid_search
import sklearn.metrics as skMetric
import scipy.stats as stats
from skimage import filters

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def RandomizedLassoRegression(np_X, np_y):
    X = np_X
    y = np_y
    X_sparse = coo_matrix(X)
    X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
    estimator = linear_model.RandomizedLasso(n_jobs=1, n_resampling=500)
    estimator.fit(X, y)
    
    return estimator.scores_

def LassoRegression(np_X, np_y, int_kOfKFold = 2, int_nJobs = 4):
    X = np_X
    y = np_y
    X_sparse = coo_matrix(X)
    X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
    kf = KFold(X.shape[0], n_folds=int_kOfKFold)
    
    list_target = []
    list_predict = []
    list_weight = []
    for idxTr, idxTe in kf:
        alpha = np.logspace(-10, 10, 200)
        parameters = [{'alpha':alpha}]
        kf_estimator = KFold(len(idxTr), n_folds=2)
        estimator_lasso = linear_model.Lasso()
        estimator_grid = grid_search.GridSearchCV(estimator_lasso, parameters, scoring='neg_mean_squared_error')
        estimator_grid.set_params(n_jobs=int_nJobs, cv=kf_estimator)
        estimator_grid.fit(X[idxTr], y[idxTr])
        list_label = estimator_grid.best_estimator_.predict(X[idxTe])
        list_weight.append([float(item) for item in estimator_grid.best_estimator_.coef_])
        for idx_y, idx_label in zip(list(y[idxTe]), list_label):
            list_target.append(float(idx_y))
            list_predict.append(idx_label)
    np_weight = np.array(list_weight)
    np_weight = np.average(list_weight, axis=0)
    float_MSE = skMetric.mean_squared_error(list_target, list_predict)
    float_pearson = stats.stats.pearsonr(list_target, list_predict)[0]
    float_spearman = stats.stats.spearmanr(list_target, list_predict)[0]
    
    return (float_pearson + float_spearman) / 2, np_weight

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def SingleGeneEpistasisLasso(str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 4):    
    ### set path of output file
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype)
    
    #-------------------------
    # load data
    #-------------------------
    ### count lines of input files
    int_num_genotype = sum(1 for line in open(str_inputFileName_genotype))
    int_num_phenotype = sum(1 for line in open(str_inputFileName_phenotype))
    
    ### get phenotype file
    list_phenotype = []
    with open(str_inputFileName_phenotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_phenotype.append(line.strip().split(","))
    np_phenotype = np.array(list_phenotype)
    del list_phenotype
    
    ### get genotype file
    list_genotype_rsid = []
    ### declare a numpy array for one-hot-encoded genotype
    np_genotype = np.empty([int_num_phenotype, int_num_genotype * 3], dtype='int')
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
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
    # preprocess data
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
    if np_genotype.shape[1] **2 * (np_genotype.shape[0] * np_genotype.itemsize * 12) > mem.available:
        return "MemErr"
    
    ### generate interaction terms
    sklearn_poly = PolynomialFeatures(degree=2, interaction_only=True, include_bias=False)
    np_genotype = sklearn_poly.fit_transform(np_genotype)
    np_genotype_rsid = np.array(sklearn_poly.get_feature_names(np_genotype_rsid))
    
    ### variance check on polynomial features
    try:
        sk_variance = VarianceThreshold(threshold=(.95 * (1 - .95)))
        np_genotype = sk_variance.fit_transform(np_genotype)
        np_genotype_rsid = np.array(np_genotype_rsid[sk_variance.get_support()])
    except:
        return 0.0
    
    ### remove redundant polynomial features
    np_genotype, np_selectedIdx = np.unique(np_genotype, axis=1, return_index=True)
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    
    #-------------------------
    # select feature
    #-------------------------
    ### f regression feature selection
    np_fRegression = -np.log10(f_regression(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
    np_selectedIdx = np.array([x > 2 for x in np_fRegression])
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    
    ### random lasso feature selection
    np_randWeight = np.array(RandomizedLassoRegression(np_genotype, np_phenotype[:, -1].astype(int)))
    ### apply otsu method to decide threshold
    float_threshold = filters.threshold_otsu(np_randWeight)
    np_selectedIdx = np.array([x >= float_threshold for x in np_randWeight])
    np_randWeight = np_randWeight[np_selectedIdx]
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    
    #-------------------------
    # build model
    #-------------------------
    float_AVG_S_P, np_weight = LassoRegression(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs)
    
    ### filter out zero-weight features
    np_selectedIdx = np.array([x != 0.0 for x in np_weight])
    np_weight = np_weight[np_selectedIdx]
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    
    #-------------------------
    # analyze result
    #-------------------------
    ### calculate student t-test p-value
    np_fRegression = -np.log10(f_regression(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
        
    ### calculate genotype frequency
    np_genotypeFreq = np.sum(np_genotype, axis=0).astype(float) / np_genotype.shape[0]
    
    #-------------------------
    # output results
    #-------------------------
    ### output statistics of features
    with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).split("_")[0] + "_Result.csv"), "w") as file_outputFile:
        file_outputFile.writelines("rsID,weight,student-t-test_log_p-value,genotype_frequency" + "\n")
        for idx_feature in range(0, np_genotype_rsid.shape[0]):
            file_outputFile.writelines(str(np_genotype_rsid[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_fRegression[idx_feature,]) + "," + str(np_genotypeFreq[idx_feature]) + "\n")
            
    ### output feature
    with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).split("_")[0] + "_Feature.csv"), "w") as file_outputFile:
        file_outputFile.writelines(",".join(np_genotype_rsid) + "\n")
        for idx_subject in range(0, np_genotype.shape[0]):
            file_outputFile.writelines(",".join(np_genotype[idx_subject, :].astype(str)) + "\n")
    
    return float_AVG_S_P

def BatchSingleGeneEpistasisLasso(str_inputFilePath_genotype, str_inputFileName_phenotype, str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 4):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.abspath(os.path.join(str_inputFilePath_genotype, os.pardir)) + "/singleGeneResult/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    ### scan all of the gen file in path
    list_genotypeFileName = []
    for str_fileName in os.listdir(str_inputFilePath_genotype):
        if ".gen" in str_fileName:
            list_genotypeFileName.append(str_fileName)
    
    ### batch PolyLogisticRegression
    int_count_gene = 0
    with open(str_outputFilePath + "All_Lasso_k" + str(int_kOfKFold) + ".csv", "w") as file_outputFile:
        file_outputFile.writelines("GeneSymbol,AVG_S_P" + "\n")
        for item in list_genotypeFileName:
            int_count_gene = int_count_gene + 1
            str_genotypeFileName = str_inputFilePath_genotype + item
            float_AVG_S_P = SingleGeneEpistasisLasso(str_genotypeFileName, str_inputFileName_phenotype, str_outputFilePath, int_kOfKFold, int_nJobs)
            file_outputFile.writelines(item.split("_")[0] + "," + str(float_AVG_S_P) + "\n")
            str_print = "step4: Processing: " + "{0:.2f}".format(float(int_count_gene) / len(list_genotypeFileName) * 100) + "% - " + item + ": " + str(float_AVG_S_P) + "\t\t"
            sys.stdout.write('%s\r' % str_print)
            sys.stdout.flush()
    
    print("step4: Detect single gene epistasis. DONE! \t\t\t\t")