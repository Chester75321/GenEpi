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
from sklearn.feature_selection import chi2
from scipy.sparse import coo_matrix
from sklearn.utils import shuffle
from sklearn import linear_model
from sklearn.cross_validation import KFold
from sklearn import grid_search
import scipy.stats as stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
import step4_singleGeneEpistasis_Logistic

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def CrossGeneEpistasisLogistic(str_inputFilePath_feature, str_inputFileName_phenotype, str_inputFileName_score = "", str_outputFilePath = ""):
    str_inputFilePath_feature = "D:\\Phd\\Grade_05\\Alzheimer\\GenEpi\\genepi\\example\\snpSubsets\\singleGeneResult\\"
    str_inputFileName_phenotype = "D:\\Phd\\Grade_05\\Alzheimer\\GenEpi\\genepi\\example\\sample.csv"
    
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.abspath(os.path.join(str_inputFilePath_feature, os.pardir)) + "/crossGeneResult/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)

    #-------------------------
    # load data
    #-------------------------
    ### load score file
    dict_score = {}
    int_count_memErr = 0
    with open(input_scoreFileName, "r") as file_inputFile:
        file_inputFile.readline()
        for line in file_inputFile:
            list_thisScore = line.strip().split(",")
            if list_thisScore[1] == "MemErr":
                int_count_memErr = int_count_memErr + 1
            elif float(list_thisScore[1]) == 0.0:
                pass
            else:
                dict_score[list_thisScore[0]] = float(list_thisScore[1])
    
    ### get all the file names of feature file
    list_featureFileName = []
    for str_fileName in os.listdir(str_inputFilePath_feature):
        if "Feature.csv" in str_fileName:
            list_featureFileName.append(str_fileName)
    
    ### get genotype SNP id
    list_genotype_rsid = []
    for item in list_featureFileName:
        with open(str_inputFilePath_feature + item, "r") as file_inputFile:
            list_rsids = file_inputFile.readline().strip().split(",")
            for rsid in list_rsids:
                list_genotype_rsid.append(rsid)
    np_genotype_rsid = np.array(list_genotype_rsid)
    
    ### count lines in file
    int_num_genotype = len(np_genotype_rsid)
    int_num_phenotype = sum(1 for line in open(str_inputFileName_phenotype))
    
    print(int_num_genotype)
    print(int_num_phenotype)
    with open(str_outputFilePath + "SNPList.csv", "w") as file_outputFile:
        for idx_SNP in range(0, len(np_genotype_rsid)):
            file_outputFile.writelines(str(np_genotype_rsid[idx_SNP]) + "\n")
    
    ### get phenotype file
    list_phenotype = []
    with open(str_inputFileName_phenotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_phenotype.append(line.strip().split(","))
    np_phenotype = np.array(list_phenotype)
    del list_phenotype
    
    ### get genotype file
    dict_geneMap ={}
    int_idx_genotype_rsid = 0
    np_genotype = np.empty([int_num_phenotype, int_num_genotype], dtype='int')
    for item in list_featureFileName:
        with open(str_inputFilePath_feature + item, "r") as file_inputFile:
            list_rsids = file_inputFile.readline().strip().split(",")
            for rsid in list_rsids:
                dict_geneMap[rsid] = item.split("_")[0]
            int_idx_phenotype = 0
            for line in file_inputFile:
                np_genotype[int_idx_phenotype, int_idx_genotype_rsid:int_idx_genotype_rsid + len(list_rsids)] = np.array([float(x) for x in line.strip().split(",")], dtype='int')
                int_idx_phenotype = int_idx_phenotype + 1
            int_idx_genotype_rsid = int_idx_genotype_rsid + len(list_rsids)
    
    #-------------------------
    # data preprocessing
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
    
    if np_genotype_degree1.shape[1] > 0:
        np_genotype_crossGene_rsid, np_genotype_crossGene = CrossGeneInteractionGenerator(np_genotype_degree1_rsid, np_genotype_degree1, np_phenotype)
        np_genotype = np_genotype
        np_genotype = np.concatenate((np_genotype, np_genotype_crossGene), axis=1)
        np_genotype_rsid = np.concatenate((np_genotype_rsid, np_genotype_crossGene_rsid))
    
    #-------------------------
    # feature selection
    #-------------------------
    ### chi-square test selection
    np_chi2 = -np.log10(chi2(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
    np_selectedIdx = np.array([x > 5 for x in np_chi2])
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0
    
    ### random logistic feature selection
    np_randWeight = np.array(STEP04_PloyLogistic.RandomizedLogisticRegressionL1(np_genotype, np_phenotype[:, -1].astype(int)))
    np_selectedIdx = np.array([x >= 0.35 for x in np_randWeight])
    np_randWeight = np_randWeight[np_selectedIdx]
    np_genotype = np_genotype[:, np_selectedIdx]
    np_genotype_rsid = np_genotype_rsid[np_selectedIdx]
    if np_genotype_rsid.shape[0] == 0:
        return 0.0    
    
    #-------------------------
    # modeling
    #-------------------------
    ### model selection
    float_f1Score, np_weight = STEP04_PloyLogistic.LogisticRegressionL1(np_genotype, np_phenotype[:, -1].astype(int))
    #print("F1 Score: " + str(float_f1Score) + "\n")
    ### model persistence
    LogisticRegressionL1ModelPersistence(np_genotype, np_phenotype[:, -1].astype(int), str_outputFilePath)
    
    ### result analysis
    np_chi2 = -np.log10(chi2(np_genotype.astype(int), np_phenotype[:, -1].astype(int))[1])
    list_oddsRatio = []
    for idx_feature in range(0, np_genotype.shape[1]):
        np_contingency = STEP04_PloyLogistic.GenerateContingencyTable(np_genotype[:, idx_feature], np_phenotype[:, -1])
        oddsratio, pvalue = stats.fisher_exact(np_contingency)
        list_oddsRatio.append(oddsratio)
    
    ### calculate allele frequency
    np_alleleFreq = np.sum(np_genotype, axis=0).astype(float) / np_genotype.shape[0]
    
    ### calculate dependecy
    np_dependency = np.empty([np_genotype.shape[1], np_genotype.shape[1]], dtype='float')
    for idx_gene in range(0, np_genotype.shape[1]):
        np_thisChi2 = -np.log10(chi2(np_genotype[:, idx_gene:].astype(int), np_genotype[:, idx_gene].astype(int))[1])
        np_dependency[idx_gene, idx_gene:] = np_thisChi2
    np_dependency[np.tril_indices(np_genotype.shape[1])] = np.transpose(np_dependency)[np.tril_indices(np_genotype.shape[1])]
    
    #-------------------------
    # outputting results
    #-------------------------
    ### statistics results
    with open(str_outputFilePath + "All_R.csv", "w") as file_outputFile:
        file_outputFile.writelines("rsID,rand_weight,weight,chi2_p-value,odds_ratio,alleleFreq,geneSymbol,singleGeneScore" + "\n")
        for idx_feature in range(0, np_genotype_rsid.shape[0]):
            if np_genotype_rsid[idx_feature,] in dict_geneMap.keys():
                str_thisOutput = str(np_genotype_rsid[idx_feature,]) + "," + str(np_randWeight[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_chi2[idx_feature,]) + "," + str(list_oddsRatio[idx_feature]) + "," + str(np_alleleFreq[idx_feature]) + "," + str(dict_geneMap[np_genotype_rsid[idx_feature,]]) + "," + str(dict_score[dict_geneMap[np_genotype_rsid[idx_feature,]]]) + "\n"
                file_outputFile.writelines(str_thisOutput)
            else:
                str_thisOutput = str(np_genotype_rsid[idx_feature,]) + "," + str(np_randWeight[idx_feature,]) + "," + str(np_weight[idx_feature,]) + "," + str(np_chi2[idx_feature,]) + "," + str(list_oddsRatio[idx_feature]) + "," + str(np_alleleFreq[idx_feature]) + "," + str(dict_geneMap[np_genotype_rsid[idx_feature,].split(" ")[0]]) + "*" + str(dict_geneMap[np_genotype_rsid[idx_feature,].split(" ")[1]]) + ", " + "\n"
                file_outputFile.writelines(str_thisOutput)
            
    ### output feature
    with open(str_outputFilePath + "All_F.csv", "w") as file_outputFile:
        file_outputFile.writelines(",".join(np_genotype_rsid) + "\n")
        for idx_subject in range(0, np_genotype.shape[0]):
            file_outputFile.writelines(",".join(np_genotype[idx_subject, :].astype(str)) + "\n")

    return float_f1Score