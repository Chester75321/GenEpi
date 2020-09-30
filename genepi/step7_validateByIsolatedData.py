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
import subprocess

from genepi.step5_crossGeneEpistasis_Logistic import PlotPolygenicScore

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def SplittingDataAsIsolatedData(str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = "", int_randomState = 0):
    ### set path of output file
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype)    
    
    ### count lines of input files
    int_num_genotype = sum(1 for line in open(str_inputFileName_genotype))
    int_num_phenotype = sum(1 for line in open(str_inputFileName_phenotype))

    ### set random state
    random = np.random.RandomState(int_randomState)
    ### random sample
    np_choice = random.choice(int_num_phenotype, int(int_num_phenotype * 0.9), replace=False)
    np_random = np.zeros(int_num_phenotype, dtype=bool)
    np_random[np_choice] = True
    np_random_complement = np.ones(int_num_phenotype, dtype=bool)
    np_random_complement[np_choice] = False
    
    list_random = [x for x in range(int_num_phenotype) if x in np_choice]
    list_random = [str(x*3 + 6) + "-" + str(x*3 + 6 + 2) for x in list_random]
    list_random_complement =  [x for x in range(int_num_phenotype) if x not in np_choice]
    list_random_complement = [str(x*3 + 6) + "-" + str(x*3 + 6 + 2) for x in list_random_complement]
    
    ### get phenotype file
    list_phenotype = []
    with open(str_inputFileName_phenotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_phenotype.append(line.strip().split(","))
    np_phenotype = np.array(list_phenotype)
    del list_phenotype

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
    
    ### get genotype file
    str_filename_subset = os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).replace(".gen", "_subset_1.gen"))
    is_cont = False
    str_choice = ""
    for x in range(int_num_phenotype):
        if is_cont == False and np_random[x]:
            is_cont = True
            str_choice += str(x*3 + 6) + "-"
        elif is_cont == True and not np_random[x]:
            is_cont = False
            str_choice += str((x-1)*3 + 6 + 2) + ","
    str_choice = str_choice[:-1] if str_choice[-1] == "," else str_choice
    
    str_command = 'cut -d " " -f1-5,'
    str_command += str_choice + " "
    str_command += str_inputFileName_genotype + " > "
    str_command += str_filename_subset
    subprocess.call(str_command, shell=True)
    
    is_cont = False
    str_choice = ""
    for x in range(int_num_phenotype):
        if is_cont == False and np_random_complement[x]:
            is_cont = True
            str_choice += str(x*3 + 6) + "-"
        elif is_cont == True and not np_random_complement[x]:
            is_cont = False
            str_choice += str((x-1)*3 + 6 + 2) + ","
    str_choice = str_choice[:-1] if str_choice[-1] == "," else str_choice

    str_filename_subset = os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).replace(".gen", "_subset_2.gen"))
    str_command = 'cut -d " " -f1-5,'
    str_command += str_choice + " "
    str_command += str_inputFileName_genotype + " > "
    str_command += str_filename_subset
    subprocess.call(str_command, shell=True)

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
        for subitem in item.replace("_AA", "").replace("_AB", "").replace("_BB", "").split("*"):
            if subitem not in dict_feature_rsid_unique:
                dict_feature_rsid_unique[subitem] = 1
    
    ### extract selected snp from genotype file
    list_inputFile_genotype = []
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_thisSnp = line.strip().split(" ")
            if list_thisSnp[1] in dict_feature_rsid_unique.keys():
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

    ### get genotype file
    list_genotype = [[] for x in range(int_num_phenotype)]
    list_genotype_rsid = []
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_thisSnp = line.strip().split(" ")
            if list_thisSnp[1] not in dict_feature_rsid_unique.keys():
                continue
            np_this_genotype = np.empty([int_num_phenotype, 3], dtype='int8')
            for idx_subject in range(0, int_num_phenotype):
                list_allelType = [0, 0, 0]
                list_allelType[np.argmax(list_thisSnp[idx_subject * 3 + 5 : idx_subject * 3 + 8])] = 1
                list_genotype[idx_subject].extend(list_allelType)
            list_genotype_rsid.append(list_thisSnp[1] + "_AA")
            list_genotype_rsid.append(list_thisSnp[1] + "_AB")
            list_genotype_rsid.append(list_thisSnp[1] + "_BB")
    np_genotype = np.array(list_genotype, dtype=np.int8)
    np_genotype_rsid = np.array(list_genotype_rsid)
    
    ### generate feature
    np_feature = np.empty([int_num_phenotype, len(list_feature_rsid_all)], dtype='int')
    for idx_feature in range(len(list_feature_rsid_all)):
        list_feature_rsid = list_feature_rsid_all[idx_feature].split("*")
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
    list_predict_proba = []
    list_label = estimator.predict(np_genotype)
    list_prob = estimator.predict_proba(np_genotype)
    for idx_y, idx_label, idx_prob in zip(list(np_phenotype[:, -1].astype(int)), list_label, list_prob):
        list_target.append(float(idx_y))
        list_predict.append(idx_label)
        list_predict_proba.append(idx_prob)
    float_f1Score = skMetric.f1_score(list_target, list_predict)
    dict_y = {"target": list_target, "predict": list_predict, "predict_proba": list_predict_proba}

    ### calculate statistic
    tn, fp, fn, tp = skMetric.confusion_matrix(dict_y["target"], dict_y["predict"]).ravel()
    float_specificity = tn / (tn + fp)
    float_sensitivity = tp / (fn + tp)

    fpr, tpr, _ = skMetric.roc_curve(dict_y["target"], np.array(dict_y["predict_proba"])[:,1])
    float_auc = skMetric.auc(fpr, tpr)

    PlotPolygenicScore(dict_y["target"], dict_y["predict"], dict_y["predict_proba"], str_outputFilePath, "ISO")
    
    print("step7: Validate by isolated data. DONE! (Isolated test score:" + "{0:.2f}".format(float_f1Score) + ")")
    print("AUC: " + "{0:.2f}".format(float_auc) + "; Specificity: " + "{0:.2f}".format(float_specificity) + "; Sensitivity: " + "{0:.2f}".format(float_sensitivity))

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
    
    print("step7: Validate by isolated data. DONE! (Isolated test score:" + "{0:.2f}".format(float_AVG_S_P) + ")")

    return float_AVG_S_P

def ValidateByIsolatedDataCovariateClassifier(str_inputFileName_model, str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype) + "/isolatedValidation/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    if not os.path.isfile(str_inputFileName_model):
        print("step7: Error no other factors exist.")
        return 0.0

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
    
    print("step7: Validate by isolated data (with covariates). DONE! (Isolated test score:" + "{0:.2f}".format(float_f1Score) + ")")
    
    return float_f1Score

def ValidateByIsolatedDataCovariateRegressor(str_inputFileName_model, str_inputFileName_feature, str_inputFileName_genotype, str_inputFileName_phenotype, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype) + "/isolatedValidation/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    if not os.path.isfile(str_inputFileName_model):
        print("step7: Error no other factors exist.")
        return 0.0
    
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
    
    print("step7: Validate by isolated data (with covariates). DONE! (Isolated test score:" + "{0:.2f}".format(float_AVG_S_P) + ")")

    return float_AVG_S_P