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

from sklearn import linear_model
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
from sklearn.externals import joblib

from genepi.step4_singleGeneEpistasis_Logistic import LogisticRegressionL1CV
from genepi.step4_singleGeneEpistasis_Lasso import LassoRegressionCV
from genepi.step5_crossGeneEpistasis_Logistic import LogisticRegressionL1
from genepi.step5_crossGeneEpistasis_Lasso import LassoRegression

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def LoadDataForEnsemble(str_inputFileName_feature, str_inputFileName_phenotype):
    """

    Loading genetic features for ensembling with covariates

    Args:
        str_inputFilePath_feature (str): File path of input feature files from stage 2 - crossGeneEpistasis
        str_inputFileName_phenotype (str): File name of input phenotype data
        
    Returns:
        (tuple): tuple containing:

            - np_genotype (ndarray): 2D array containing genotype data with `int8` type
            - np_phenotype (ndarray): 2D array containing phenotype data with `float` type
    
    """

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
        print("step6: There is no other factors exist.")
        return None, None
    
    ### get genotype file
    np_genotype = np.empty([int_num_phenotype, int_num_genotype], dtype=np.float16)
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

def ClassifierModelPersistence(np_X, np_y, str_outputFilePath = "", int_nJobs = 4):
    """

    Dumping ensemble classifier for model persistence

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
    
    cost = [2**x for x in range(-8, 8)]
    parameters = [{'C':cost, 'penalty':['l1'], 'dual':[False], 'class_weight':['balanced']}]
    kf_estimator = KFold(n_splits=2)
    estimator_logistic = linear_model.LogisticRegression(max_iter=100, solver='liblinear')
    estimator_grid = GridSearchCV(estimator_logistic, parameters, scoring='f1', n_jobs=int_nJobs, cv=kf_estimator)
    estimator_grid.fit(X, y)
    
    joblib.dump(estimator_grid.best_estimator_, os.path.join(str_outputFilePath, "Classifier_Covariates.pkl"))

def RegressorModelPersistence(np_X, np_y, str_outputFilePath = "", int_nJobs = 4):
    """

    Dumping ensemble regressor for model persistence

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
    
    joblib.dump(estimator_grid.best_estimator_, os.path.join(str_outputFilePath, "Regressor_Covariates.pkl"))

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def EnsembleWithCovariatesClassifier(str_inputFileName_feature, str_inputFileName_phenotype, str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 4):
    """

    A workflow to ensemble genetic features with covariates for L1-regularized Logistic regression.

    Args:
        str_inputFilePath_feature (str): File path of input feature files from stage 2 - crossGeneEpistasis
        str_inputFileName_phenotype (str): File name of input phenotype data
        str_outputFilePath (str): File path of output file
        int_kOfKFold (int): The k for k-fold cross validation (default: 2)
        int_nJobs (int): The number of thread (default: 4)

    Returns:
        (tuple): tuple containing:

            - float_f1Score_train (float): The F1 score of the model for training set
            - float_f1Score_test (float): The F1 score of the model for testing set

        - Expected Success Response::

            "step6: Ensemble with covariates. DONE!"
    
    """
    
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_feature)
    
    #-------------------------
    # load data
    #-------------------------
    np_genotype, np_phenotype = LoadDataForEnsemble(str_inputFileName_feature, str_inputFileName_phenotype)
    if np_genotype is None and np_phenotype is None:
        return 0.0, 0.0
    
    #-------------------------
    # build model
    #-------------------------
    float_f1Score_test, np_weight = LogisticRegressionL1CV(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs)
    float_f1Score_train = LogisticRegressionL1(np_genotype, np_phenotype[:, -1].astype(int), int_nJobs)
    
    #-------------------------
    # dump persistent model
    #-------------------------
    ClassifierModelPersistence(np_genotype, np_phenotype[:, -1].astype(int), str_outputFilePath, int_nJobs)
    
    print("step6: Ensemble with covariates. DONE! (Training score:" + "{0:.2f}".format(float_f1Score_train) + "; " + str(int_kOfKFold) + "-fold Test Score:" + "{0:.2f}".format(float_f1Score_test) + ")")

    return float_f1Score_train, float_f1Score_test

def EnsembleWithCovariatesRegressor(str_inputFileName_feature, str_inputFileName_phenotype, str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 4):
    """

    A workflow to ensemble genetic features with covariates for L1-regularized Lasso regression.

    Args:
        str_inputFilePath_feature (str): File path of input feature files from stage 2 - crossGeneEpistasis
        str_inputFileName_phenotype (str): File name of input phenotype data
        str_outputFilePath (str): File path of output file
        int_kOfKFold (int): The k for k-fold cross validation (default: 2)
        int_nJobs (int): The number of thread (default: 4)

    Returns:
        (tuple): tuple containing:

            - float_AVG_S_P_train (float): The average of the Peason's and Spearman's correlation of the model for training set
            - float_AVG_S_P_test (float): The average of the Peason's and Spearman's correlation of the model for testing set

        - Expected Success Response::

            "step6: Ensemble with covariates. DONE!"
    
    """
    
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_feature)
    
    #-------------------------
    # load data
    #-------------------------
    np_genotype, np_phenotype = LoadDataForEnsemble(str_inputFileName_feature, str_inputFileName_phenotype)
    if np_genotype is None and np_phenotype is None:
        return 0.0, 0.0
    
    #-------------------------
    # build model
    #-------------------------
    float_AVG_S_P_test, np_weight = LassoRegressionCV(np_genotype, np_phenotype[:, -1].astype(float), int_kOfKFold, int_nJobs)
    float_AVG_S_P_train = LassoRegression(np_genotype, np_phenotype[:, -1].astype(float), int_nJobs)
    
    #-------------------------
    # dump persistent model
    #-------------------------
    RegressorModelPersistence(np_genotype, np_phenotype[:, -1].astype(int), str_outputFilePath, int_nJobs)
    
    print("step6: Ensemble with covariates. DONE! (Training score:" + "{0:.2f}".format(float_AVG_S_P_train) + "; " + str(int_kOfKFold) + "-fold Test Score:" + "{0:.2f}".format(float_AVG_S_P_test) + ")")
    
    return float_AVG_S_P_train, float_AVG_S_P_test