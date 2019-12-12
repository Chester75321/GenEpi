# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester (Yu-Chuan Chang)
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import os
import warnings
warnings.filterwarnings('ignore')
# ignore all warnings
warnings.simplefilter("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"

import os
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from sklearn.feature_selection import f_regression
from sklearn import linear_model
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
from sklearn.externals import joblib
import scipy.stats as stats

from genepi.step4_singleGeneEpistasis_Lasso import RandomizedLassoRegression
from genepi.step4_singleGeneEpistasis_Lasso import LassoRegressionCV
from genepi.step4_singleGeneEpistasis_Lasso import FeatureEncoderLasso

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def LassoRegression(np_X, np_y, int_nJobs = 4):
    """

    Implementation of the L1-regularized Lasso regression with k-fold cross validation.

    Args:
        np_X (ndarray): 2D array containing genotype data with `int8` type
        np_y (ndarray): 2D array containing phenotype data with `float` type
        int_nJobs (int): The number of thread (default: 4)

    Returns:
        (float): float_AVG_S_P
        
            The average of the Peason's and Spearman's correlation of the model
    
    """

    X = np_X
    y = np_y
    
    list_target = []
    list_predict = []
    
    alpha = np.logspace(-10, 10, 200)
    parameters = [{'alpha':alpha}]
    kf_estimator = KFold(n_splits=2)
    estimator_lasso = linear_model.Lasso(max_iter=1000)
    estimator_grid = GridSearchCV(estimator_lasso, parameters, scoring='neg_mean_squared_error', n_jobs=int_nJobs, cv=kf_estimator)
    estimator_grid.fit(X, y)
    list_label = estimator_grid.best_estimator_.predict(X)
    for idx_y, idx_label in zip(list(y), list_label):
        list_target.append(float(idx_y))
        list_predict.append(idx_label)
    float_pearson = stats.stats.pearsonr(list_target, list_predict)[0]
    float_spearman = stats.stats.spearmanr(list_target, list_predict)[0]
    
    return (float_pearson + float_spearman) / 2

def RegressorModelPersistence(np_X, np_y, str_outputFilePath = "", int_nJobs = 4):
    """

    Dumping regressor for model persistence

    Args:
        np_X (ndarray): 2D array containing genotype data with `int8` type
        np_y (ndarray): 2D array containing phenotype data with `float` type
        str_outputFilePath (str): File path of output file
        int_nJobs (int): The number of thread (default: 4)

    Returns:
        None
    
    """

    X = np_X
    y = np_y
    X_sparse = coo_matrix(X)
    X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
    
    alpha = np.logspace(-10, 10, 200)
    parameters = [{'alpha':alpha}]
    kf_estimator = KFold(n_splits=2)
    estimator_lasso = linear_model.Lasso(max_iter=1000)
    estimator_grid = GridSearchCV(estimator_lasso, parameters, scoring='neg_mean_squared_error', n_jobs=int_nJobs, cv=kf_estimator)
    estimator_grid.fit(X, y)
    
    joblib.dump(estimator_grid.best_estimator_, os.path.join(str_outputFilePath, "Regressor.pkl"))

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def CrossGeneEpistasisLasso(str_inputFilePath_feature, str_inputFileName_phenotype, str_inputFileName_score = "", str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 4):
    """

    A workflow to model a cross gene epistasis containing two-element combinatorial encoding, stability selection, filtering low quality varaint and  L1-regularized Lasso regression with k-fold cross validation.

    Args:
        str_inputFilePath_feature (str): File path of input feature files from stage 1 - singleGeneEpistasis
        str_inputFileName_phenotype (str): File name of input phenotype data
        str_inputFileName_score (str): File name of input score file from stage 1 - singleGeneEpistasis
        str_outputFilePath (str): File path of output file
        int_kOfKFold (int): The k for k-fold cross validation (default: 2)
        int_nJobs (int): The number of thread (default: 4)

    Returns:
        (tuple): tuple containing:

            - float_AVG_S_P_train (float): The average of the Peason's and Spearman's correlation of the model for training set
            - float_AVG_S_P_test (float): The average of the Peason's and Spearman's correlation of the model for testing set

        - Expected Success Response::

            "step5: Detect cross gene epistasis. DONE!"
    
    """
    
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.abspath(os.path.join(str_inputFilePath_feature, os.pardir)) + "/crossGeneResult/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    ### set default score file name
    if str_inputFileName_score == "":
        for str_fileName in os.listdir(str_inputFilePath_feature):
            if str_fileName.startswith("All_Lasso"):
                str_inputFileName_score = os.path.join(str_inputFilePath_feature, str_fileName)
    
    #-------------------------
    # load data
    #-------------------------
    ### scan score file and exclude useless genes
    dict_score = {}
    with open(str_inputFileName_score, "r") as file_inputFile:
        file_inputFile.readline()
        for line in file_inputFile:
            list_thisScore = line.strip().split(",")
            if list_thisScore[1] == "MemErr" or float(list_thisScore[1]) == 0.0:
                pass
            else:
                dict_score[list_thisScore[0]] = float(list_thisScore[1])
    
    ### get all the file names of feature file
    list_featureFileName = []
    for str_fileName in os.listdir(str_inputFilePath_feature):
        if "Feature.csv" in str_fileName:
            list_featureFileName.append(str_fileName)
    
    ### get all selected snp ids
    list_genotype_rsid = []
    for item in list_featureFileName:
        with open(os.path.join(str_inputFilePath_feature, item), "r") as file_inputFile:
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
    np_phenotype = np.array(list_phenotype, dtype=np.float)
    del list_phenotype
    
    ### get genotype file
    ### declare a dictionary for mapping snp and gene
    dict_geneMap ={}
    idx_genotype_rsid = 0
    np_genotype = np.empty([int_num_phenotype, int_num_genotype], dtype='int8')
    for item in list_featureFileName:
        with open(os.path.join(str_inputFilePath_feature, item), "r") as file_inputFile:
            ### grep feature from header of feature file
            list_rsids = file_inputFile.readline().strip().split(",")
            for rsid in list_rsids:
                ### key: rsIDs of a feature; value: gene symbol
                dict_geneMap[rsid] = item.split("_")[0]
            idx_phenotype = 0
            ### read feaure and write into np_genotype
            for line in file_inputFile:
                np_genotype[idx_phenotype, idx_genotype_rsid:idx_genotype_rsid + len(list_rsids)] = np.array([float(x) for x in line.strip().split(",")], dtype='int')
                idx_phenotype = idx_phenotype + 1
            idx_genotype_rsid = idx_genotype_rsid + len(list_rsids)
    
    #-------------------------
    # preprocess data
    #-------------------------
    ### f regression feature selection
    np_fRegression = -np.log10(f_regression(np_genotype.astype(int), np_phenotype[:, -1].astype(float))[1])
    np_selectedIdx = np.array([x > 5 for x in np_fRegression])
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0

    ### select degree 1 feature
    np_genotype_rsid_degree = np.array([str(x).count('*') + 1 for x in np_genotype_rsid])
    np_selectedIdx = np.array([x == 1 for x in np_genotype_rsid_degree])
    np_genotype_degree1 = np_genotype[:, np_selectedIdx]
    np_genotype_degree1_rsid = np_genotype_rsid[np_selectedIdx]
    
    ### remove redundant polynomial features
    np_genotype_degree1, np_selectedIdx = np.unique(np_genotype_degree1, axis=1, return_index=True)
    np_genotype_degree1_rsid = np_genotype_degree1_rsid[np_selectedIdx]
    
    ### generate cross gene interations
    np_genotype_crossGene_rsid, np_genotype_crossGene = FeatureEncoderLasso(np_genotype_degree1_rsid, np_genotype_degree1, np_phenotype, 1)
    
    ### remove degree 1 feature from dataset
    np_selectedIdx = np.array([x != 1 for x in np_genotype_rsid_degree])
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    
    ### concatenate cross gene interations
    if np_genotype_degree1.shape[1] > 0:
        np_genotype = np.concatenate((np_genotype, np_genotype_crossGene), axis=1)
        np_genotype_rsid = np.concatenate((np_genotype_rsid, np_genotype_crossGene_rsid))
    
    #-------------------------
    # select feature
    #-------------------------
    ### random lasso feature selection
    np_randWeight = np.array(RandomizedLassoRegression(np_genotype, np_phenotype[:, -1].astype(float)))
    np_selectedIdx = np.array([x >= 0.1 for x in np_randWeight])
    np_randWeight = np_randWeight[np_selectedIdx]
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    
    #-------------------------
    # build model
    #-------------------------
    float_AVG_S_P_test, np_weight = LassoRegressionCV(np_genotype, np_phenotype[:, -1].astype(float), int_kOfKFold, int_nJobs)
    float_AVG_S_P_train = LassoRegression(np_genotype, np_phenotype[:, -1].astype(float), int_nJobs)
    
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
    np_fRegression = -np.log10(f_regression(np_genotype.astype(int), np_phenotype[:, -1].astype(float))[1])
        
    ### calculate genotype frequency
    np_genotypeFreq = np.sum(np_genotype, axis=0).astype(float) / np_genotype.shape[0]
    
    #-------------------------
    # output results
    #-------------------------
    ### output statistics of features
    with open(os.path.join(str_outputFilePath, "Result.csv"), "w") as file_outputFile:
        file_outputFile.writelines("rsid,weight,student-t-test_log_p-value,genotype_frequency,geneSymbol,singleGeneScore" + "\n")
        for idx_feature in range(0, np_genotype_rsid.shape[0]):
            ### if this feature is single gene epistasis
            if np_genotype_rsid[idx_feature,] in dict_geneMap.keys():
                str_thisOutput = str(np_genotype_rsid[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_fRegression[idx_feature,]) + "," + str(np_genotypeFreq[idx_feature]) + "," + str(dict_geneMap[np_genotype_rsid[idx_feature,]]).split("@")[0] + "," + str(dict_score[dict_geneMap[np_genotype_rsid[idx_feature,]]]) + "\n"
                file_outputFile.writelines(str_thisOutput)
            ### else this feature is cross gene epistasis
            else:
                str_thisOutput = str(np_genotype_rsid[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_fRegression[idx_feature,]) + "," + str(np_genotypeFreq[idx_feature]) + "," + str(dict_geneMap[np_genotype_rsid[idx_feature,].split("*")[0]]).split("@")[0] + "*" + str(dict_geneMap[np_genotype_rsid[idx_feature,].split("*")[1]]).split("@")[0] + ", " + "\n"
                file_outputFile.writelines(str_thisOutput)

    ### output feature
    with open(os.path.join(str_outputFilePath, "Feature.csv"), "w") as file_outputFile:
        file_outputFile.writelines(",".join(np_genotype_rsid) + "\n")
        for idx_subject in range(0, np_genotype.shape[0]):
            file_outputFile.writelines(",".join(np_genotype[idx_subject, :].astype(str)) + "\n")

    #-------------------------
    # dump persistent model
    #-------------------------
    RegressorModelPersistence(np_genotype, np_phenotype[:, -1].astype(int), str_outputFilePath, int_nJobs)

    print("step5: Detect cross gene epistasis. DONE! (Training score:" + "{0:.2f}".format(float_AVG_S_P_train) + "; " + str(int_kOfKFold) + "-fold Test Score:" + "{0:.2f}".format(float_AVG_S_P_test) + ")")
    
    return float_AVG_S_P_train, float_AVG_S_P_test
    
