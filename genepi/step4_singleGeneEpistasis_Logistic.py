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
import sys
import itertools
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import chi2
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn import linear_model
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
import sklearn.metrics as skMetric
import scipy.stats as stats
import multiprocessing as mp

from genepi.tools import randomized_l1

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def RandomizedLogisticRegression(np_X, np_y):
    """

    Implementation of the stability selection.

    Args:
        np_X (ndarray): 2D array containing genotype data with `int8` type
        np_y (ndarray): 2D array containing phenotype data with `float` type

    Returns:
        (ndarray): estimator.scores 
        
            1D array containing the scores of each genetic features with `float` type
    
    """

    X = np_X
    y = np_y
    X_sparse = coo_matrix(X)
    X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
    estimator = randomized_l1.RandomizedLogisticRegression(n_jobs=1, n_resampling=100)
    estimator.fit(X, y)
    
    return estimator.scores_

def LogisticRegressionL1CV(np_X, np_y, int_kOfKFold = 2, int_nJobs = 1):
    """

    Implementation of the L1-regularized Logistic regression with k-fold cross validation.

    Args:
        np_X (ndarray): 2D array containing genotype data with `int8` type
        np_y (ndarray): 2D array containing phenotype data with `float` type
        int_kOfKFold (int): The k for k-fold cross validation (default: 2)
        int_nJobs (int): The number of thread (default: 1)

    Returns:
        (ndarray): estimator.scores 
        
            1D array containing the scores of each genetic features with `float` type
    
    """

    X = np_X
    y = np_y
    X_sparse = coo_matrix(X)
    X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
    kf = KFold(n_splits=int_kOfKFold)
    
    list_target = []
    list_predict = []
    list_weight = []
    for idxTr, idxTe in kf.split(X):
        cost = [2**x for x in range(-8, 8)]
        parameters = [{'C':cost, 'penalty':['l1'], 'dual':[False], 'class_weight':['balanced']}]
        kf_estimator = KFold(n_splits=2)
        estimator_logistic = linear_model.LogisticRegression(max_iter=100, solver='liblinear')
        estimator_grid = GridSearchCV(estimator_logistic, parameters, scoring='f1', n_jobs=1, cv=kf_estimator)
        estimator_grid.fit(X[idxTr], y[idxTr])
        list_label = estimator_grid.best_estimator_.predict(X[idxTe])
        list_weight.append([float(item) for item in estimator_grid.best_estimator_.coef_[0]])
        for idx_y, idx_label in zip(list(y[idxTe]), list_label):
            list_target.append(float(idx_y))
            list_predict.append(idx_label)
    np_weight = np.array(list_weight)
    np_weight = np.average(list_weight, axis=0)
    float_f1Score = skMetric.f1_score(list_target, list_predict)
    
    return float_f1Score, np_weight

def FeatureEncoderLogistic(np_genotype_rsid, np_genotype, np_phenotype, int_dim):
    """

    Implementation of the two-element combinatorial encoding.

    Args:
        np_genotype_rsid (ndarray): 1D array containing rsid of genotype data with `str` type
        np_genotype (ndarray): 2D array containing genotype data with `int8` type
        np_phenotype (ndarray): 2D array containing phenotype data with `float` type
        int_dim (int): The dimension of a variant (default: 3. AA, AB and BB)

    Returns:
        (tuple): tuple containing:

            - list_interaction_rsid (ndarray): 1D array containing rsid of genotype data with `str` type
            - np_interaction (ndarray): 2D array containing genotype data with `int8` type
    
    """

    ### combinatorial encoding
    np_interaction = np_genotype
    list_interaction_rsid = list(np_genotype_rsid)    

    list_combs = list(itertools.combinations(range(int(np_interaction.shape[1]/int_dim)), 2))
    for idx_combs in range(len(list_combs)):
        try:
            ### generate interaction terms
            tuple_comb = list_combs[idx_combs]
            np_this_interaction = np.zeros([np_phenotype.shape[0], int_dim**2], dtype='int8')
            list_this_interaction_id = []
            for idx_x in range(int_dim):
                for idx_y in range(int_dim):
                    np_this_interaction_term = (np_genotype[:, tuple_comb[0] * int_dim + idx_x] * np_genotype[:, tuple_comb[1] * int_dim + idx_y]).astype(np.int8)
                    if not(np.array_equal(np_this_interaction_term, np_genotype[:, tuple_comb[0] * int_dim + idx_x])) and not(np.array_equal(np_this_interaction_term, np_genotype[:, tuple_comb[1] * int_dim + idx_y])):
                        np_this_interaction[:, idx_x * int_dim + idx_y] = np_this_interaction_term
                    list_this_interaction_id.append(np_genotype_rsid[tuple_comb[0] * int_dim + idx_x] + "*" + np_genotype_rsid[tuple_comb[1] * int_dim + idx_y])
            
            ### variance check (detect variance < 0.05)
            sk_variance = VarianceThreshold(threshold=(.95 * (1 - .95)))
            np_this_interaction = sk_variance.fit_transform(np_this_interaction)
            np_this_interaction_id = np.array(list_this_interaction_id)
            np_this_interaction_id = np.array(np_this_interaction_id[sk_variance.get_support()])
            
            ### chi-square test selection
            np_chi2 = -np.log10(chi2(np_this_interaction.astype(int), np_phenotype[:, -1].astype(int))[1])
            np_selectedIdx = np.array([x > 2 for x in np_chi2])
            np_this_interaction = np_this_interaction[:, np_selectedIdx]
            np_this_interaction_id = np_this_interaction_id[np_selectedIdx]
        
            ### append insteraction terms
            int_num_interaction = np_this_interaction.shape[1]
            if int_num_interaction == 0:
                continue
            np_interaction_append = np.empty((np_interaction.shape[0], np_interaction.shape[1] + int_num_interaction), dtype='int')
            np_interaction_append[:,:-(int_num_interaction)] = np_interaction
            np_interaction_append[:,-(int_num_interaction):] = np_this_interaction
            np_interaction = np_interaction_append
            list_interaction_rsid.extend(list(np_this_interaction_id))
        except:
            pass

    return np.array(list_interaction_rsid), np_interaction

def GenerateContingencyTable(np_genotype, np_phenotype):
    """

    Generating the contingency table for chi-square test.

    Args:
        np_X (ndarray): 2D array containing genotype data with `int8` type
        np_y (ndarray): 2D array containing phenotype data with `float` type

    Returns:
        (ndarray): np_contingency 
        
            2D array containing the contingency table with `int` type
    
    """

    np_contingency = np.array([[0, 0], [0, 0]])
    for idx_subject in range(0, np_genotype.shape[0]):
        np_contingency[int(np_genotype[idx_subject]), int(np_phenotype[idx_subject])] = np_contingency[int(np_genotype[idx_subject]), int(np_phenotype[idx_subject])] + 1
    np_contingency = np.rot90(np_contingency)
    np_contingency = np.rot90(np_contingency)
    
    return np_contingency

def FilterInLoading(np_genotype, np_phenotype):
    """

    This function is for filtering low quality varaint. Before modeling each subset of genotype features, two criteria were adopted to exclude low quality data. The first criterion is that the genotype frequency of a feature should exceed 5%, where the genotype frequency means the proportion of genotype among the total samples in the dataset. The second criterion is regarding the association between the feature and the phenotype. We used χ2 test to estimate the association between the feature and the phenotype, and the p-value should be smaller than 0.01.

    Args:
        np_genotype (ndarray): 2D array containing genotype data with `int8` type
        np_phenotype (ndarray): 2D array containing phenotype data with `float` type

    Returns:
        (ndarray): np_genotype
        
            2D array containing genotype data with `int8` type
    
    """

    try:
        ### variance check (detect variance < 0.05)
        sk_variance = VarianceThreshold(threshold=(.95 * (1 - .95)))
        np_genotype = sk_variance.fit_transform(np_genotype)
        
        ### chi-square test selection
        np_chi2 = -np.log10(chi2(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
        np_selectedIdx = np.array([x > 2 for x in np_chi2])
        np_genotype = np_genotype[:, np_selectedIdx]
        
        return np_genotype.shape[1]
    
    except:
        return 0

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def SingleGeneEpistasisLogistic(str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 1):      
    """

    A workflow to model a single gene containing two-element combinatorial encoding, stability selection, filtering low quality varaint and  L1-regularized Logistic regression with k-fold cross validation.

    Args:
        str_inputFileName_genotype (str): File name of input genotype data
        str_inputFileName_phenotype (str): File name of input phenotype data
        str_outputFilePath (str): File path of output file
        int_kOfKFold (int): The k for k-fold cross validation (default: 2)
        int_nJobs (int): The number of thread (default: 1)

    Returns:
        (float): float_f1Score
        
            The F1 score of the model
    
    """
    
    ### set path of output file
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype)
    
    #-------------------------
    # load data
    #-------------------------
    ### count lines of input files
    int_num_phenotype = sum(1 for line in open(str_inputFileName_phenotype))
    
    ### get phenotype file
    list_phenotype = []
    with open(str_inputFileName_phenotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_phenotype.append(line.strip().split(","))
    np_phenotype = np.array(list_phenotype, dtype=np.float)
    del list_phenotype
    
    ### get genotype file
    list_genotype = [[] for x in range(int_num_phenotype)]
    list_genotype_rsid = []
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_thisSnp = line.strip().split(" ")
            np_this_genotype = np.empty([int_num_phenotype, 3], dtype='int8')
            for idx_subject in range(0, int_num_phenotype):
                list_allelType = [0, 0, 0]
                list_allelType[np.argmax(list_thisSnp[idx_subject * 3 + 5 : idx_subject * 3 + 8])] = 1
                np_this_genotype[idx_subject, :] = list_allelType
            if FilterInLoading(np_this_genotype, np_phenotype) == 0:
                continue
            for idx_subject in range(0, int_num_phenotype):
                list_allelType = [0, 0, 0]
                list_allelType[np.argmax(list_thisSnp[idx_subject * 3 + 5 : idx_subject * 3 + 8])] = 1
                list_genotype[idx_subject].extend(list_allelType)
            list_genotype_rsid.append(list_thisSnp[1] + "_AA")
            list_genotype_rsid.append(list_thisSnp[1] + "_AB")
            list_genotype_rsid.append(list_thisSnp[1] + "_BB")
    np_genotype = np.array(list_genotype, dtype=np.int8)
    np_genotype_rsid = np.array(list_genotype_rsid)
    
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    
    #-------------------------
    # preprocess data
    #-------------------------        
    ### generate interaction terms
    np_genotype_rsid, np_genotype = FeatureEncoderLogistic(np_genotype_rsid, np_genotype, np_phenotype, 3)
    
    #-------------------------
    # select feature
    #-------------------------
    ### random logistic feature selection
    np_randWeight = np.array(RandomizedLogisticRegression(np_genotype, np_phenotype[:, -1].astype(int)))
    np_selectedIdx = np.array([x >= 0.25 for x in np_randWeight])
    np_randWeight = np_randWeight[np_selectedIdx]
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    
    #-------------------------
    # build model
    #-------------------------
    float_f1Score, np_weight = LogisticRegressionL1CV(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs)
    if float_f1Score == 0.0:
        return 0.0
    
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
    ### calculate chi-square p-value
    np_chi2 = -np.log10(chi2(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
    list_oddsRatio = []
    for idx_feature in range(0, np_genotype.shape[1]):
        np_contingency = GenerateContingencyTable(np_genotype[:, idx_feature], np_phenotype[:, -1])
        oddsratio, pvalue = stats.fisher_exact(np_contingency)
        list_oddsRatio.append(oddsratio)
        
    ### calculate genotype frequency
    np_genotypeFreq = np.sum(np_genotype, axis=0).astype(float) / np_genotype.shape[0]
    
    #-------------------------
    # output results
    #-------------------------
    ### output statistics of features
    with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).split("_")[0] + "_Result.csv"), "w") as file_outputFile:
        file_outputFile.writelines("rsid,weight,chi-square_log_p-value,odds_ratio,genotype_frequency" + "\n")
        for idx_feature in range(0, np_genotype_rsid.shape[0]):
            file_outputFile.writelines(str(np_genotype_rsid[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_chi2[idx_feature,]) + "," + str(list_oddsRatio[idx_feature]) + "," + str(np_genotypeFreq[idx_feature]) + "\n")
            
    ### output feature
    with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).split("_")[0] + "_Feature.csv"), "w") as file_outputFile:
        file_outputFile.writelines(",".join(np_genotype_rsid) + "\n")
        for idx_subject in range(0, np_genotype.shape[0]):
            file_outputFile.writelines(",".join(np_genotype[idx_subject, :].astype(str)) + "\n")
    
    return float_f1Score

def BatchSingleGeneEpistasisLogistic(str_inputFilePath_genotype, str_inputFileName_phenotype, str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = mp.cpu_count()):
    """

    Batch running for the single gene workflow.

    Args:
        str_inputFilePath_genotype (str): File path of input genotype data
        str_inputFileName_phenotype (str): File name of input phenotype data
        str_outputFilePath (str): File path of output file
        int_kOfKFold (int): The k for k-fold cross validation (default: 2)
        int_nJobs (int): The number of thread (default: 1)

    Returns:
        - Expected Success Response::

            "step4: Detect single gene epistasis. DONE!"
    
    """

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
    ### inital multiprocessing pool
    mp_pool = mp.Pool(int_nJobs)
    
    ### apply pool on the function that need be parallelizing
    dict_result = {}
    for int_count_gene, float_f1Score in enumerate(mp_pool.starmap(SingleGeneEpistasisLogistic, [(os.path.join(str_inputFilePath_genotype, gene), str_inputFileName_phenotype, str_outputFilePath, int_kOfKFold, int_nJobs) for gene in list_genotypeFileName]), 0):
        if list_genotypeFileName[int_count_gene] not in dict_result:
            dict_result[list_genotypeFileName[int_count_gene]] = float_f1Score
        str_print = "step4: Processing: " + "{0:.2f}".format(float(int_count_gene) / len(list_genotypeFileName) * 100) + "% - " + list_genotypeFileName[int_count_gene] + ": " + "\t\t"
        sys.stdout.write('%s\r' % str_print)
        sys.stdout.flush()

    mp_pool.close()

    ### output result
    with open(str_outputFilePath + "All_Logistic_k" + str(int_kOfKFold) + ".csv", "w") as file_outputFile:
        file_outputFile.writelines("GeneSymbol,F1Score" + "\n")
        for key, value in dict_result.items():
            file_outputFile.writelines(key.split("_")[0] + "," + str(value) + "\n")

    '''
    ### batch PolyLogisticRegression
    int_count_gene = 0
    with open(str_outputFilePath + "All_Logistic_k" + str(int_kOfKFold) + ".csv", "w") as file_outputFile:
        file_outputFile.writelines("GeneSymbol,F1Score" + "\n")
        for item in list_genotypeFileName:
            int_count_gene = int_count_gene + 1
            str_genotypeFileName = os.path.join(str_inputFilePath_genotype, item)
            float_f1Score = SingleGeneEpistasisLogistic(str_genotypeFileName, str_inputFileName_phenotype, str_outputFilePath, int_kOfKFold, int_nJobs)
            file_outputFile.writelines(item.split("_")[0] + "," + str(float_f1Score) + "\n")
            str_print = "step4: Processing: " + "{0:.2f}".format(float(int_count_gene) / len(list_genotypeFileName) * 100) + "% - " + item + ": " + str(float_f1Score) + "\t\t"
            sys.stdout.write('%s\r' % str_print)
            sys.stdout.flush()
    '''
    
    print("step4: Detect single gene epistasis. DONE! \t\t\t\t")