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
from sklearn.preprocessing import PolynomialFeatures
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import f_regression
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn.cross_validation import KFold
import scipy.stats as stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.externals import joblib
from .step4_singleGeneEpistasis_Lasso import RandomizedLassoRegression
from .step4_singleGeneEpistasis_Lasso import LassoRegression

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def CrossGeneInteractionGenerator(np_genotype_rsid, np_genotype):
    ### generate interaction terms
    sklearn_poly = PolynomialFeatures(degree=2, interaction_only=True, include_bias=False)
    np_genotype = sklearn_poly.fit_transform(np_genotype)
    np_genotype_rsid = np.array(sklearn_poly.get_feature_names(np_genotype_rsid))
    
    ### variance check on ploynomial features
    sk_variance = VarianceThreshold(threshold=(.95 * (1 - .95)))
    np_genotype = sk_variance.fit_transform(np_genotype)
    np_genotype_rsid = np.array(np_genotype_rsid[sk_variance.get_support()])
    
    ### remove redundant polynomial features
    np_genotype, np_selectedIdx = np.unique(np_genotype, axis=1, return_index=True)
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    
    return np_genotype_rsid, np_genotype

def RFRegressor(np_X, np_y, int_kOfKFold = 2, int_nJobs = 4):
    X = np_X
    y = np_y
    X_sparse = coo_matrix(X)
    X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
    kf = KFold(X.shape[0], n_folds=int_kOfKFold)
    
    list_target = []
    list_predict = []
    for idxTr, idxTe in kf:
        estimator_rf = RandomForestRegressor(n_estimators=1000, n_jobs=int_nJobs)
        estimator_rf.fit(X[idxTr], y[idxTr])
        list_label = estimator_rf.predict(X[idxTe])
        for idx_y, idx_label in zip(list(y[idxTe]), list_label):
            list_target.append(float(idx_y))
            list_predict.append(idx_label)
    float_pearson = stats.stats.pearsonr(list_target, list_predict)[0]
    float_spearman = stats.stats.spearmanr(list_target, list_predict)[0]

    return (float_pearson + float_spearman) / 2

def RFRegressorModelPersistence(np_X, np_y, str_outputFilePath = "", str_outputFileName = "RFRegressor.pkl", int_nJobs = 4):
    X = np_X
    y = np_y
    X_sparse = coo_matrix(X)
    X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
    
    estimator_rf = RandomForestRegressor(n_estimators=1000, n_jobs=int_nJobs)
    estimator_rf.fit(X, y)
    list_predict = estimator_rf.predict(X)
    float_pearson = stats.stats.pearsonr(y, list_predict)[0]
    float_spearman = stats.stats.spearmanr(y, list_predict)[0]
  
    joblib.dump(estimator_rf, os.path.join(str_outputFilePath, str_outputFileName))
    
    return (float_pearson + float_spearman) / 2

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def CrossGeneEpistasisLasso(str_inputFilePath_feature, str_inputFileName_phenotype, str_inputFileName_score = "", str_outputFilePath = "", int_kOfKFold = 2, int_nJobs = 4):
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
        with open(str_inputFilePath_feature + item, "r") as file_inputFile:
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
    
    ### get genotype file
    ### declare a dictionary for mapping snp and gene
    dict_geneMap ={}
    idx_genotype_rsid = 0
    np_genotype = np.empty([int_num_phenotype, int_num_genotype], dtype='int')
    for item in list_featureFileName:
        with open(str_inputFilePath_feature + item, "r") as file_inputFile:
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
    ### select degree 1 feature
    np_genotype_rsid_degree = np.array([str(x).count('_') for x in np_genotype_rsid])
    np_selectedIdx = np.array([x == 1 for x in np_genotype_rsid_degree])
    np_genotype_degree1 = np_genotype[:, np_selectedIdx]
    np_genotype_degree1_rsid = np_genotype_rsid[np_selectedIdx]
    
    ### remove degree 1 feature from dataset
    np_selectedIdx = np.array([x != 1 for x in np_genotype_rsid_degree])
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    
    ### generate cross gene interations
    if np_genotype_degree1.shape[1] > 0:
        np_genotype_crossGene_rsid, np_genotype_crossGene = CrossGeneInteractionGenerator(np_genotype_degree1_rsid, np_genotype_degree1)
        np_genotype = np.concatenate((np_genotype, np_genotype_crossGene), axis=1)
        np_genotype_rsid = np.concatenate((np_genotype_rsid, np_genotype_crossGene_rsid))
    
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
    np_randWeight = np.array(RandomizedLassoRegression(np_genotype, np_phenotype[:, -1].astype(int), int_nJobs))
    np_selectedIdx = np.array([x >= 0.25 for x in np_randWeight])
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
    
    ### use random forest classifier for ensemble
    float_AVG_S_P_rf_te = RFRegressor(np_genotype, np_phenotype[:, -1].astype(int), int_kOfKFold, int_nJobs)
    ### random forest classifier model persistence
    float_AVG_S_P_rf_tr = RFRegressorModelPersistence(np_genotype, np_phenotype[:, -1].astype(int), str_outputFilePath)
    
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
    with open(str_outputFilePath + "Result.csv", "w") as file_outputFile:
        file_outputFile.writelines("rsid,weight,student-t-test_log_p-value,genotype_frequency,geneSymbol,singleGeneScore" + "\n")
        for idx_feature in range(0, np_genotype_rsid.shape[0]):
            ### if this feature is single gene epistasis
            if np_genotype_rsid[idx_feature,] in dict_geneMap.keys():
                str_thisOutput = str(np_genotype_rsid[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_fRegression[idx_feature,]) + "," + str(np_genotypeFreq[idx_feature]) + "," + str(dict_geneMap[np_genotype_rsid[idx_feature,]]) + "," + str(dict_score[dict_geneMap[np_genotype_rsid[idx_feature,]]]) + "\n"
                file_outputFile.writelines(str_thisOutput)
            ### else this feature is cross gene epistasis
            else:
                str_thisOutput = str(np_genotype_rsid[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_fRegression[idx_feature,]) + "," + str(np_genotypeFreq[idx_feature]) + "," + str(dict_geneMap[np_genotype_rsid[idx_feature,].split(" ")[0]]) + "*" + str(dict_geneMap[np_genotype_rsid[idx_feature,].split(" ")[1]]) + ", " + "\n"
                file_outputFile.writelines(str_thisOutput)
            
    ### output feature
    with open(str_outputFilePath + "Feature.csv", "w") as file_outputFile:
        file_outputFile.writelines(",".join(np_genotype_rsid) + "\n")
        for idx_subject in range(0, np_genotype.shape[0]):
            file_outputFile.writelines(",".join(np_genotype[idx_subject, :].astype(str)) + "\n")

    print("step5: Detect cross gene epistasis. DONE! (Training score:" + "{0:.2f}".format(float_AVG_S_P_rf_tr) + "; " + str(int_kOfKFold) + "-fold Test Score:" + "{0:.2f}".format(float_AVG_S_P_rf_te) + ")")
    
    return float_AVG_S_P
    