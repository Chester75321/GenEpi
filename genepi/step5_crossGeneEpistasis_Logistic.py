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
import scipy as sp
np.seterr(divide='ignore', invalid='ignore')
from sklearn.feature_selection import chi2
from sklearn import linear_model
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
from sklearn.externals import joblib
import sklearn.metrics as skMetric
import scipy.stats as stats

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from genepi.step4_singleGeneEpistasis_Logistic import RandomizedLogisticRegression
from genepi.step4_singleGeneEpistasis_Logistic import LogisticRegressionL1CV
from genepi.step4_singleGeneEpistasis_Logistic import FeatureEncoderLogistic

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def LogisticRegressionL1(np_X, np_y, int_nJobs = 4):
    """

    Implementation of the L1-regularized Logistic regression with k-fold cross validation.

    Args:
        np_X (ndarray): 2D array containing genotype data with `int8` type
        np_y (ndarray): 2D array containing phenotype data with `float` type
        int_nJobs (int): The number of thread (default: 4)

    Returns:
        (float): float_f1Score
        
            The F1 score of the model
    
    """

    X = np_X
    y = np_y
    
    list_target = []
    list_predict = []
    
    cost = [2**x for x in range(-8, 8)]
    parameters = [{'C':cost, 'penalty':['l1'], 'dual':[False], 'class_weight':['balanced']}]
    kf_estimator = KFold(n_splits=2)
    estimator_logistic = linear_model.LogisticRegression(max_iter=100, solver='liblinear')
    estimator_grid = GridSearchCV(estimator_logistic, parameters, scoring='f1', n_jobs=int_nJobs, cv=kf_estimator)
    estimator_grid.fit(X, y)
    list_label = estimator_grid.best_estimator_.predict(X)
    for idx_y, idx_label in zip(list(y), list_label):
        list_target.append(float(idx_y))
        list_predict.append(idx_label)
    float_f1Score = skMetric.f1_score(list_target, list_predict)
    
    return float_f1Score

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

def ClassifierModelPersistence(np_X, np_y, str_outputFilePath = "", int_nJobs = 4):
    """

    Dumping classifier for model persistence

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
    
    joblib.dump(estimator_grid.best_estimator_, os.path.join(str_outputFilePath, "Classifier.pkl"))

def fsigmoid(x, a, b):
    float_return = 1.0/(1.0+np.exp(-a*(x-b)))
    return float_return

def PlotPolygenicScore(np_X, np_y, int_kOfKFold = 2, int_nJobs = 1, str_outputFilePath=""):
    """

    Plot figure for polygenic score, including group distribution and prevalence to PGS

    Args:
        np_X (ndarray): 2D array containing genotype data with `int8` type
        np_y (ndarray): 2D array containing phenotype data with `float` type
        int_kOfKFold (int): The k for k-fold cross validation (default: 2)
        int_nJobs (int): The number of thread (default: 1)

    Returns:
        None
    
    """

    X = np_X
    y = np_y
    X_sparse = coo_matrix(X)
    X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
    kf = KFold(n_splits=int_kOfKFold)

    list_target = []
    list_predict = []
    list_proba = []
    for idxTr, idxTe in kf.split(X):
        cost = [2**x for x in range(-8, 8)]
        parameters = [{'C':cost, 'penalty':['l1'], 'dual':[False], 'class_weight':['balanced']}]
        kf_estimator = KFold(n_splits=2)
        estimator_logistic = linear_model.LogisticRegression(max_iter=100, solver='liblinear')
        estimator_grid = GridSearchCV(estimator_logistic, parameters, scoring='f1', n_jobs=1, cv=kf_estimator)
        estimator_grid.fit(X[idxTr], y[idxTr])
        list_predict_this = estimator_grid.best_estimator_.predict(X[idxTe])
        list_proba_this = estimator_grid.best_estimator_.predict_proba(X[idxTe])
        for idx_y, idx_predict, idx_proba in zip(list(y[idxTe]), list_predict_this, list_proba_this):
            list_target.append(float(idx_y))
            list_predict.append(idx_predict)
            list_proba.append(idx_proba)
        
    float_f1Score = skMetric.f1_score(list_target, list_predict)

    #-------------------------
    # group distribution
    #-------------------------
    pd_pgs = pd.concat([pd.DataFrame(list_target), pd.DataFrame(list_predict), pd.DataFrame(np.array(list_proba)[:,1])], axis=1)
    pd_pgs.columns = ['target', 'predict', 'proba']

    plt.figure(figsize=(5,5))
    pd_case = pd_pgs[pd_pgs.target == 1.0]
    sns.distplot(pd_case['proba'],  kde=True, bins=25, label='Case', color="#c9474b")
    pd_control = pd_pgs[pd_pgs.target == 0.0]
    sns.distplot(pd_control['proba'],  kde=True, bins=25, label='Control', color="#4e7e91")

    # plot formatting
    str_method = "GenEpi"
    plt.legend(prop={'size': 12})
    plt.title(str_method + ' Predicting F1 Score: ' + "%.4f" % float_f1Score + ' ')
    plt.xlim(0, 1)
    plt.ylim(0, 15)
    plt.xlabel('Polygenic Score')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(str_outputFilePath, "GenEpi_PGS.png"), dpi=300)
    plt.close('all')

    #-------------------------
    # prevalence to PGS
    #-------------------------
    int_step = 5
    np_percentile = np.percentile(pd_pgs['proba'], q=list(range(0, 100 + int_step, int_step)))
    pd_pgs['percentile'] = np.searchsorted(np_percentile, pd_pgs['proba'], side='left') * int_step

    pd_prevalence_obs = pd_pgs.groupby('percentile').sum()[['target']]
    pd_prevalence_obs_count = pd_pgs.groupby('percentile').count()[['target']]
    pd_prevalence_obs = pd_prevalence_obs/pd_prevalence_obs_count
    pd_prevalence_obs.columns = ['obs']
    pd_prevalence_pre = pd_pgs.groupby('percentile').mean()[['proba']]
    pd_prevalence_pre.columns = ['pre']

    pd_rr_ingroup_case = pd_pgs.groupby('percentile').sum()[['target']]
    pd_rr_ingroup_sum = pd_pgs.groupby('percentile').count()[['target']]
    pd_rr = (pd_rr_ingroup_case / pd_rr_ingroup_sum) / (pd_pgs.sum()[['target']] / pd_pgs.count()[['target']])
    pd_rr.columns = ['Relative Risk']

    plt.figure(figsize=(5,5))
    sns.scatterplot(x=pd_prevalence_obs.index, y=pd_prevalence_obs['obs'], hue=pd_rr['Relative Risk'], palette=sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True))
    popt, pcov = sp.optimize.curve_fit(fsigmoid, pd_prevalence_pre.index, pd_prevalence_pre['pre'], method='dogbox', bounds=([0., 0.],[1., 100.]))
    sns.lineplot(x=pd_prevalence_pre.index, y=fsigmoid(pd_prevalence_pre.index, *popt), color="black")

    plt.legend(prop={'size': 12}, loc='upper left')
    plt.title(str_method + ' Predicting F1 Score: ' + "%.4f" % float_f1Score + ' ')
    plt.xlim(0, 100)
    plt.ylim(0, 1)
    plt.xlabel('Polygenic Score Percentile')
    plt.ylabel('Prevalence of Percentile Group')
    plt.savefig(os.path.join(str_outputFilePath, "GenEpi_Prevalence.png"), dpi=300)
    plt.close('all')

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def CrossGeneEpistasisLogistic(str_inputFilePath_feature, str_inputFileName_phenotype, str_inputFileName_score = "", str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 4):   
    """

    A workflow to model a cross gene epistasis containing two-element combinatorial encoding, stability selection, filtering low quality varaint and  L1-regularized Logistic regression with k-fold cross validation.

    Args:
        str_inputFilePath_feature (str): File path of input feature files from stage 1 - singleGeneEpistasis
        str_inputFileName_phenotype (str): File name of input phenotype data
        str_inputFileName_score (str): File name of input score file from stage 1 - singleGeneEpistasis
        str_outputFilePath (str): File path of output file
        int_kOfKFold (int): The k for k-fold cross validation (default: 2)
        int_nJobs (int): The number of thread (default: 4)

    Returns:
        (tuple): tuple containing:

            - float_f1Score_train (float): The F1 score of the model for training set
            - float_f1Score_test (float): The F1 score of the model for testing set

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
            if str_fileName.startswith("All_Logistic"):
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
    ### chi-square test selection
    np_chi2 = -np.log10(chi2(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
    np_selectedIdx = np.array([x > 5 for x in np_chi2])
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
    np_genotype_crossGene_rsid, np_genotype_crossGene = FeatureEncoderLogistic(np_genotype_degree1_rsid, np_genotype_degree1, np_phenotype, 1)
    
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
    float_f1Score_test, np_weight = LogisticRegressionL1CV(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs)
    float_f1Score_train = LogisticRegressionL1(np_genotype, np_phenotype[:, -1].astype(int), int_nJobs)
    
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
    with open(os.path.join(str_outputFilePath, "Result.csv"), "w") as file_outputFile:
        file_outputFile.writelines("rsid,weight,chi-square_log_p-value,odds_ratio,genotype_frequency,geneSymbol,singleGeneScore" + "\n")
        for idx_feature in range(0, np_genotype_rsid.shape[0]):
            ### if this feature is single gene epistasis
            if np_genotype_rsid[idx_feature,] in dict_geneMap.keys():
                str_thisOutput = str(np_genotype_rsid[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_chi2[idx_feature,]) + "," + str(list_oddsRatio[idx_feature]) + "," + str(np_genotypeFreq[idx_feature]) + "," + str(dict_geneMap[np_genotype_rsid[idx_feature,]]).split("@")[0] + "," + str(dict_score[dict_geneMap[np_genotype_rsid[idx_feature,]]]) + "\n"
                file_outputFile.writelines(str_thisOutput)
            ### else this feature is cross gene epistasis
            else:
                str_thisOutput = str(np_genotype_rsid[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_chi2[idx_feature,]) + "," + str(list_oddsRatio[idx_feature]) + "," + str(np_genotypeFreq[idx_feature]) + "," + str(dict_geneMap[np_genotype_rsid[idx_feature,].split("*")[0]]).split("@")[0] + "*" + str(dict_geneMap[np_genotype_rsid[idx_feature,].split("*")[1]]).split("@")[0] + ", " + "\n"
                file_outputFile.writelines(str_thisOutput)            
 
    ### output feature
    with open(os.path.join(str_outputFilePath, "Feature.csv"), "w") as file_outputFile:
        file_outputFile.writelines(",".join(np_genotype_rsid) + "\n")
        for idx_subject in range(0, np_genotype.shape[0]):
            file_outputFile.writelines(",".join(np_genotype[idx_subject, :].astype(str)) + "\n")

    ### output figures
    PlotPolygenicScore(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs, str_outputFilePath)

    #-------------------------
    # dump persistent model
    #-------------------------
    ClassifierModelPersistence(np_genotype, np_phenotype[:, -1].astype(int), str_outputFilePath, int_nJobs)
    
    print("step5: Detect cross gene epistasis. DONE! (Training score:" + "{0:.2f}".format(float_f1Score_train) + "; " + str(int_kOfKFold) + "-fold Test Score:" + "{0:.2f}".format(float_f1Score_test) + ")")
    
    return float_f1Score_train, float_f1Score_test
