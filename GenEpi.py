# -*- coding: utf-8 -*-
"""
Created on Apr 2019

@author: Chester (Yu-Chuan Chang)
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import argparse
import time
import os
import sys
import genepi

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def ArgumentsParser():
    ### define arguments
    str_description = ''
    'GenEpi is a package to uncover epistasis associated with phenotypes by a machine learning approach, '
    'developed by Yu-Chuan Chang at c4Lab of National Taiwan University. '
    'Github: https://github.com/Chester75321/GenEpi '
    'Reference: GenEpi: Gene-based Epistasis Discovery Using Machine Learning '
    '(https://www.biorxiv.org/content/10.1101/421719v5) '
    parser = argparse.ArgumentParser(prog='GenEpi', description=str_description)
    
    ### define arguments for I/O
    parser.add_argument("-g", required=True, help="filename of the input .gen file")
    parser.add_argument("-p", required=True, help="filename of the input phenotype")
    parser.add_argument("-s", required=False, help="self-defined genome regions")
    parser.add_argument("-o", required=False, help="output file path")
    
    ### define arguments for modeling
    parser.add_argument("-m", required=False, default="c", choices=["c", "r"], help="choose model type: c for classification; r for regression")
    parser.add_argument("-k", required=False, default=2, help="k of k-fold cross validation")
    parser.add_argument("-t", required=False, default=1, help="number of threads")
    
    ### define arguments for step1_downloadUCSCDB
    parser_group_1 = parser.add_argument_group("update UCSC database")
    parser_group_1.add_argument('--updatedb', action='store_true', default=False, help="enable this function")
    parser_group_1.add_argument("-b", required=False, default="hg19", choices=["hg19", "hg38"], help="human genome build")
    
    ### define arguments for step2_estimateLD
    parser_group_2 = parser.add_argument_group("compress data by LD block")
    parser_group_2.add_argument('--compressld', action='store_true', default=False, help="enable this function")
    parser_group_2.add_argument("-d", required = False, default=0.9, type=float, help="threshold for compression: D prime")
    parser_group_2.add_argument("-r", required = False, default=0.9, type=float, help="threshold for compression: R square")
 
    return parser

def InputChecking(str_inputFileName_genotype, str_inputFileName_phenotype):
    ### check file name exist
    if str_inputFileName_genotype is None:
        sys.exit("There is no input genotype file.")
    if str_inputFileName_phenotype is None:
        sys.exit("There is no input phenotype file.")
    
    ### count lines of input files
    int_num_genotype = sum(1 for line in open(str_inputFileName_genotype))
    int_num_phenotype = sum(1 for line in open(str_inputFileName_phenotype))
    
    ### count sample num. in genotype file
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
        list_line = file_inputFile.readline().strip().split(" ")
        int_num_genotype_sample = (len(list_line) - 5) / 3
        if int_num_genotype_sample != int_num_phenotype:
            sys.exit("The number of samples in genotype file does not match the number of samples in phenotype file.")
    
    return int_num_genotype, int_num_phenotype

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def main(args=None):
    ### obtain arguments from argument parser
    args = ArgumentsParser().parse_args(args)
    
    ### get arguments for I/O
    str_inputFileName_genotype = args.g
    str_inputFileName_phenotype = args.p
    str_inputFileName_regions = ""
    if args.s is not None:
        str_inputFileName_regions = args.s
    else:
        str_inputFileName_regions = "None"
    str_outputFilePath = ""
    if args.o is not None:
        str_outputFilePath = args.o
    else:
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype)
        
    if str_inputFileName_genotype == "example" and str_inputFileName_phenotype == "example":
        str_command = "cp " + os.path.join(os.path.dirname(genepi.__file__), "example", "sample.csv") + " " + str_outputFilePath
        os.system(str_command)
        str_command = "cp " + os.path.join(os.path.dirname(genepi.__file__), "example", "sample.gen") + " " + str_outputFilePath
        os.system(str_command)
        str_command = "GenEpi -g " + os.path.join(str_outputFilePath, "sample.gen") + " "
        str_command += "-p " + os.path.join(str_outputFilePath, "sample.csv") + " "
        str_command += "-o " + str_outputFilePath + " "
        str_command += "--updatedb --compressld"
        os.system(str_command)
        return
    
    with open(os.path.join(str_outputFilePath, "GenEpi_Log_" + time.strftime("%Y%m%d-%H%M", time.localtime()) + ".txt"), "w") as file_outputFile:
        ### create log
        file_outputFile.writelines("start analysis at: " + time.strftime("%Y%m%d-%H:%M:%S", time.localtime()) + "\n" + "\n")
    
        ### log arguments
        file_outputFile.writelines("Arguments in effect:" + "\n")
        file_outputFile.writelines("\t" + "-g (input genotype filename): " + str_inputFileName_genotype + "\n")
        file_outputFile.writelines("\t" + "-p (input phenotype filename): " + str_inputFileName_phenotype + "\n")
        file_outputFile.writelines("\t" + "-s (self-defined genome regions): " + str_inputFileName_regions + "\n")
        file_outputFile.writelines("\t" + "-o (output filepath): " + str_outputFilePath + "\n" + "\n")
        
        file_outputFile.writelines("\t" + "-m (model type): " + "Classification" if args.m=="c" else "Regression"  + "\n")
        file_outputFile.writelines("\t" + "-k (k-fold cross validation): " + str(args.k) + "\n")
        file_outputFile.writelines("\t" + "-t (number of threads): " + str(args.t) + "\n" + "\n")
        
        file_outputFile.writelines("\t" + "--updatedb (enable function of update UCSC database): " + str(args.updatedb) + "\n")
        file_outputFile.writelines("\t" + "-b (human genome build): " + args.b + "\n" + "\n")
        
        file_outputFile.writelines("\t" + "--compressld (enable function of LD data compression): " + str(args.compressld) + "\n")
        file_outputFile.writelines("\t" + "-d (D prime threshold): " + str(args.d) + "\n")
        file_outputFile.writelines("\t" + "-r (R square threshold): " + str(args.r) + "\n" + "\n")
        
        ### check input format
        int_num_genotype, int_num_phenotype = InputChecking(str_inputFileName_genotype, str_inputFileName_phenotype)
        file_outputFile.writelines("Number of variants: " + str(int_num_genotype) + "\n")
        file_outputFile.writelines("Number of samples: " + str(int_num_phenotype) + "\n")
        print("Number of variants: " + str(int_num_genotype))
        print("Number of samples: " + str(int_num_phenotype))
        
        ### step1_downloadUCSCDB
        if args.updatedb:
            genepi.DownloadUCSCDB(str_hgbuild=args.b)
    
        ### step2_estimateLD
        if args.compressld:
            genepi.EstimateLDBlock(str_inputFileName_genotype, str_outputFilePath=str_outputFilePath, float_threshold_DPrime=float(args.d), float_threshold_RSquare=float(args.r))
            str_inputFileName_genotype = os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype.replace(".gen", "_LDReduced.gen")))
        
        ### step3_splitByGene
        if str_inputFileName_regions == "None":
            genepi.SplitByGene(str_inputFileName_genotype, str_outputFilePath=os.path.join(str_outputFilePath, "snpSubsets"))
        else:
            genepi.SplitByGene(str_inputFileName_genotype, str_inputFileName_UCSCDB=str_inputFileName_regions, str_outputFilePath=os.path.join(str_outputFilePath, "snpSubsets"))
        
        if args.m=="c":
            ### step4_singleGeneEpistasis_Logistic (for case/control trial)
            genepi.BatchSingleGeneEpistasisLogistic(os.path.join(str_outputFilePath, "snpSubsets"), str_inputFileName_phenotype, int_kOfKFold=int(args.k), int_nJobs=int(args.t))
            ### step5_crossGeneEpistasis_Logistic (for case/control trial)
            float_score_train, float_score_test = genepi.CrossGeneEpistasisLogistic(os.path.join(str_outputFilePath, "singleGeneResult"), str_inputFileName_phenotype, int_kOfKFold=int(args.k), int_nJobs=int(args.t))
            file_outputFile.writelines("Overall genetic feature performance (F1 score)" + "\n")
            file_outputFile.writelines("Training: " + str(float_score_train) + "\n")
            file_outputFile.writelines("Testing (" + str(args.k) + "-fold CV): " + str(float_score_test) + "\n" + "\n")
            ### step6_ensembleWithCovariates (for case/control trial)
            float_score_train, float_score_test = genepi.EnsembleWithCovariatesClassifier(os.path.join(str_outputFilePath, "crossGeneResult", "Feature.csv"), str_inputFileName_phenotype, int_kOfKFold=int(args.k), int_nJobs=int(args.t))
            file_outputFile.writelines("Ensemble with co-variate performance (F1 score)" + "\n")
            file_outputFile.writelines("Training: " + str(float_score_train) + "\n")
            file_outputFile.writelines("Testing (" + str(args.k) + "-fold CV): " + str(float_score_test) + "\n" + "\n")
        else:
            ### step4_singleGeneEpistasis_Lasso (for quantitative trial)
            genepi.BatchSingleGeneEpistasisLasso(os.path.join(str_outputFilePath, "snpSubsets"), str_inputFileName_phenotype, int_kOfKFold=int(args.k), int_nJobs=int(args.t))
            ### step5_crossGeneEpistasis_Lasso (for quantitative trial)
            float_score_train, float_score_test = genepi.CrossGeneEpistasisLasso(os.path.join(str_outputFilePath, "singleGeneResult"), str_inputFileName_phenotype, int_kOfKFold=int(args.k), int_nJobs=int(args.t))
            file_outputFile.writelines("Overall genetic feature performance (Average of the Pearson and Spearman correlation)" + "\n")
            file_outputFile.writelines("Training: " + str(float_score_train) + "\n")
            file_outputFile.writelines("Testing (" + str(args.k) + "-fold CV): " + str(float_score_test) + "\n" + "\n")
            ### step6_ensembleWithCovariates (for quantitative trial)
            float_score_train, float_score_test = genepi.EnsembleWithCovariatesRegressor(os.path.join(str_outputFilePath, "crossGeneResult", "Feature.csv"), str_inputFileName_phenotype, int_kOfKFold=int(args.k), int_nJobs=int(args.t))
            file_outputFile.writelines("Ensemble with co-variate performance (Average of the Pearson and Spearman correlation)" + "\n")
            file_outputFile.writelines("Training: " + str(float_score_train) + "\n")
            file_outputFile.writelines("Testing (" + str(args.k) + "-fold CV): " + str(float_score_test) + "\n" + "\n")
        
        file_outputFile.writelines("end analysis at: " + time.strftime("%Y%m%d-%H:%M:%S", time.localtime()) + "\n")

if __name__ == "__main__":
    main()