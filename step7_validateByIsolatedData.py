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
from sklearn.externals import joblib
import sklearn.metrics as skMetric
import scipy.stats as stats

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def SplittingDataAsIsolatedData(str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = "", int_randomState = 0):
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
    list_genotype_info = []
    list_genotype = []
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_thisSnp = line.strip().split(" ")
            list_genotype_info.append(list_thisSnp[:5])
            list_genotype.append(list_thisSnp[5:])
    np_genotype_info = np.array(list_genotype_info)
    del list_genotype_info
    np_genotype = np.array(list_genotype)
    del list_genotype
    
    ### set random state
    random = np.random.RandomState(int_randomState)
    ### random sample
    np_random = random.choice(int_num_phenotype, int(int_num_phenotype/2), replace=False)
    np_random_complement = np.ones(int_num_phenotype, dtype=bool)
    np_random_complement[np_random] = False
    
    #-------------------------
    # output results
    #-------------------------
    ### output phenotype files
    with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_phenotype).replace(".csv", "_subset_1.csv")), "w") as file_outputFile:
        np_phenotype_selected = np_phenotype[np_random, :]
        for idx_phenotype in range(np_phenotype_selected.shape[0]):
            str_line = ",".join(np_phenotype_selected[idx_phenotype, :])
            file_outputFile.writelines(str_line + "\n")
    with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_phenotype).replace(".csv", "_subset_2.csv")), "w") as file_outputFile:
        np_phenotype_selected = np_phenotype[np_random_complement, :]
        for idx_phenotype in range(np_phenotype_selected.shape[0]):
            str_line = ",".join(np_phenotype_selected[idx_phenotype, :])
            file_outputFile.writelines(str_line + "\n")
    
    ### output genotype files
    with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).replace(".gen", "_subset_1.gen")), "w") as file_outputFile:
        np_genotype_selected = np_genotype[:, np.array([[x * 3, x * 3 + 1, x * 3 + 2] for x in np_random]).ravel()]
        for idx_genotype in range(int_num_genotype):
            str_line = " ".join(np_genotype_info[idx_genotype, :]) + " " + " ".join(np_genotype_selected[idx_genotype, :])
            file_outputFile.writelines(str_line + "\n")
    with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).replace(".gen", "_subset_2.gen")), "w") as file_outputFile:
        np_genotype_selected = np_genotype[:, np.array([[x, x, x] for x in np_random_complement]).ravel()]
        for idx_genotype in range(int_num_genotype):
            str_line = " ".join(np_genotype_info[idx_genotype, :]) + " " + " ".join(np_genotype_selected[idx_genotype, :])
            file_outputFile.writelines(str_line + "\n")

def IsolatedDataFeatureGenerator(str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype) + "/isolatedValidation/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)

    ### get all selected snp ids
    list_feature_rsid_all = []
    with open(str_inputFileName_feature, "r") as file_inputFile:
        ### grep the header
        list_rsids = file_inputFile.readline().strip().split(",")
        for rsid in list_rsids:
            list_feature_rsid_all.append(rsid)
    ### get unique selected snp ids
    dict_feature_rsid_unique = {}
    for item in list_feature_rsid_all:
        for subitem in item.replace("_AA", "").replace("_AB", "").replace("_BB", "").split(" "):
            if subitem not in dict_feature_rsid_unique:
                dict_feature_rsid_unique[subitem] = 1
    
    ### extract selected snp from genotype file
    list_inputFile_genotype = []
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_thisSnp = line.strip().split(" ")
            if list_thisSnp[1] in dict_feature_rsid_unique:
                list_inputFile_genotype.append(line)
    
    ### count lines of input files
    int_num_genotype = len(list_inputFile_genotype)
    int_num_phenotype = sum(1 for line in open(str_inputFileName_phenotype))
    
    ### get phenotype file
    list_phenotype = []
    with open(str_inputFileName_phenotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_phenotype.append(line.strip().split(","))
    np_phenotype = np.array(list_phenotype)
    del list_phenotype
    
    list_genotype_rsid = []
    ### declare a numpy array for one-hot-encoded genotype
    np_genotype = np.empty([int_num_phenotype, int_num_genotype * 3], dtype='int')
    idx_snp = 0
    for line in list_inputFile_genotype:
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
    
    ### generate feature
    np_feature = np.empty([int_num_phenotype, len(list_feature_rsid_all)], dtype='int')
    for idx_feature in range(len(list_feature_rsid_all)):
        list_feature_rsid = list_feature_rsid_all[idx_feature].split(" ")
        if len(list_feature_rsid) == 1:
            np_feature[:, idx_feature] = np_genotype[:, int(np.argwhere(np_genotype_rsid == list_feature_rsid[0]))]
        else:
            np_feature[:, idx_feature] = np.multiply(np_genotype[:, int(np.argwhere(np_genotype_rsid == list_feature_rsid[0]))], np_genotype[:, int(np.argwhere(np_genotype_rsid == list_feature_rsid[1]))])
    
    ### output feature
    with open(os.path.join(str_outputFilePath, "Feature.csv"), "w") as file_outputFile:
        file_outputFile.writelines(",".join(list_feature_rsid_all) + "\n")
        for idx_subject in range(int_num_phenotype):
            file_outputFile.writelines(",".join(np_feature[idx_subject, :].astype(str)) + "\n")
    
    return np_feature, np_phenotype

def ValidateByIsolatedDataClassifier(str_inputFileName_model, str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype) + "/isolatedValidation/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    estimator = joblib.load(str_inputFileName_model)
    np_genotype, np_phenotype = IsolatedDataFeatureGenerator(str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath)
    
    list_target = []
    list_predict = []
    list_label = estimator.predict(np_genotype)
    for idx_target, idx_prdict in zip(list(np_phenotype[:, -1].astype(int)), list_label):
        list_target.append(float(idx_target))
        list_predict.append(idx_prdict)
    float_f1Score = skMetric.f1_score(list_target, list_predict)
    
    print("step7: Validate by isolated data. DONE! (Test score:" + "{0:.2f}".format(float_f1Score) + ")")
    
    return float_f1Score

def ValidateByIsolatedDataRegressor(str_inputFileName_model, str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype) + "/isolatedValidation/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    estimator = joblib.load(str_inputFileName_model)
    np_genotype, np_phenotype = IsolatedDataFeatureGenerator(str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath)
    
    list_target = []
    list_predict = []
    list_label = estimator.predict(np_genotype)
    for idx_target, idx_prdict in zip(list(np_phenotype[:, -1].astype(int)), list_label):
        list_target.append(float(idx_target))
        list_predict.append(idx_prdict)
    float_pearson = stats.stats.pearsonr(list_target, list_predict)[0]
    float_spearman = stats.stats.spearmanr(list_target, list_predict)[0]
    float_AVG_S_P = (float_pearson + float_spearman) / 2
    
    print("step7: Validate by isolated data. DONE! (Test score:" + "{0:.2f}".format(float_AVG_S_P) + ")")

    return float_AVG_S_P

def ValidateByIsolatedDataCovariateClassifier(str_inputFileName_model, str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype) + "/isolatedValidation/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    estimator = joblib.load(str_inputFileName_model)
    np_genotype, np_phenotype = IsolatedDataFeatureGenerator(str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath)
    
    if np_phenotype.shape[1] < 2:
        print("step7: Error no other factors exist.")
        return 0.0
    
    ### concatenate genotype and other factors
    np_genotype = np.concatenate((np_genotype, np_phenotype[:, :-1]), axis=1).astype(float)
    
    list_target = []
    list_predict = []
    list_label = estimator.predict(np_genotype)
    for idx_target, idx_prdict in zip(list(np_phenotype[:, -1].astype(int)), list_label):
        list_target.append(float(idx_target))
        list_predict.append(idx_prdict)
    float_f1Score = skMetric.f1_score(list_target, list_predict)
    
    print("step7: Validate by isolated data(with covariates). DONE! (Test score:" + "{0:.2f}".format(float_f1Score) + ")")
    
    return float_f1Score

def ValidateByIsolatedDataCovariateRegressor(str_inputFileName_model, str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype) + "/isolatedValidation/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    estimator = joblib.load(str_inputFileName_model)
    np_genotype, np_phenotype = IsolatedDataFeatureGenerator(str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath)
    
    if np_phenotype.shape[1] < 2:
        print("step7: Error no other factors exist.")
        return 0.0
    
    ### concatenate genotype and other factors
    np_genotype = np.concatenate((np_genotype, np_phenotype[:, :-1]), axis=1).astype(float)
    
    list_target = []
    list_predict = []
    list_label = estimator.predict(np_genotype)
    for idx_target, idx_prdict in zip(list(np_phenotype[:, -1].astype(int)), list_label):
        list_target.append(float(idx_target))
        list_predict.append(idx_prdict)
    float_pearson = stats.stats.pearsonr(list_target, list_predict)[0]
    float_spearman = stats.stats.spearmanr(list_target, list_predict)[0]
    float_AVG_S_P = (float_pearson + float_spearman) / 2
    
    print("step7: Validate by isolated data(with covariates). DONE! (Test score:" + "{0:.2f}".format(float_AVG_S_P) + ")")

    return float_AVG_S_P