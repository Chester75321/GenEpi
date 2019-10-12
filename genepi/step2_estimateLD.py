# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester (Yu-Chuan Chang)
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import os
import sys
import numpy as np

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def EstimateAlleleFrequency(gen_snp):
    """

    A function for estimating allele frequency of a single varaint

    Args:
        gen_snp (list): The genotypes of a variant of all samples
        
    Returns:
        (tuple): tuple containing:

            - float_frequency_A (float): The reference allele type frequency
            - float_frequency_B (float): The alternative allele type frequency
    
    """
    
    ### get all subject's genotype
    list_snp = gen_snp.split(" ")[5:]
    ### get the number of subjects
    int_num_subject = int(len(list_snp)/3)
    
    ### generate count table (AA, AB, BB)
    list_count = [0, 0, 0]
    for idx_subject in range(0, int_num_subject):
        idx_col = np.argmax(list_snp[idx_subject * 3: idx_subject * 3 + 3])
        list_count[idx_col] = list_count[idx_col] + 1
    
    ### calculate allele frequency
    ### frequency of A = AA + AB/2
    ### frequency of B = BB + AB/2
    float_frequency_A = float(list_count[0] + list_count[1] / 2) / int_num_subject
    float_frequency_B = float(list_count[2] + list_count[1] / 2) / int_num_subject
    
    return float_frequency_A, float_frequency_B

def EstimatePairwiseLD(gen_snp_1, gen_snp_2):
    """

    Lewontin (1964) linkage disequilibrium (LD) estimation.

    Args:
        gen_snp_1 (list): The genotypes of first variant of all samples
        gen_snp_2 (list): The genotypes of second variant of all samples

    Returns:
        (tuple): tuple containing:

            - float_D_prime (float): The DPrime of these two variants
            - float_R_square (float): The RSquare of these two variants
    
    """

    ### get all subject's genotype
    list_snp1 = gen_snp_1.split(" ")[5:]
    list_snp2 = gen_snp_2.split(" ")[5:]
    ### get the number of subjects
    int_num_subject = int(len(list_snp1)/3)
    
    ### generate contigency table
    ### row: SNP1_AA; SNP1_Aa; SNP1_aa
    ### col: SNP2_bb; SNP2_Bb; SNP2_bb
    np_contigency = np.zeros((3, 3), np.dtype(int))
    for idx_subject in range(0, int_num_subject):
        idx_row = np.argmax(np.array(list_snp1[idx_subject * 3: idx_subject * 3 + 3]))
        idx_col = np.argmax(np.array(list_snp2[idx_subject * 3: idx_subject * 3 + 3]))
        np_contigency[idx_row, idx_col] = np_contigency[idx_row, idx_col] + 1
    
    ### estimate single locus haplotyes
    ### snp1_A = (AABB + AABb + AAbb) + (AaBB + AaBb + Aabb)/2; snp1_a = snp1_A - 1
    float_probability_A = float(np.sum(np_contigency[0, :]) + float(np.sum(np_contigency[1, :])) / 2) / int_num_subject
    float_probability_a = 1 - float_probability_A
    ### snp2_B = (AABB + AaBB + aaBB) + (AABb + AaBb + aaBb)/2; snp2_b = snp2_B - 1
    float_probability_B = float(np.sum(np_contigency[:, 0]) + float(np.sum(np_contigency[:, 1])) / 2) / int_num_subject
    float_probability_b = 1 - float_probability_B
    
    ### set arbitrary probability of AB
    float_probability_AB = float_probability_A * float_probability_B
    
    try:
        ### EM algorithm
        for idx_loop in range(0, 10000):
            ### E(num_AB|prob_AB) = 2 * num_AABB + num_AABb + num_AaBB +
            ### (prob_AB * (1 + prob_AB - prob_A - prob_B) * num_AbBb) / 
            ### ((prob_A - prob_AB) * (prob_B - prob_AB) + prob_AB * (1 + prob_AB - prob_A - prob_B))
            float_num_AB_estimateByEM = 2 * float(np_contigency[0, 0]) + float(np_contigency[0, 1]) + float(np_contigency[1, 0]) + (float_probability_AB * (1 + float_probability_AB - float_probability_A - float_probability_B) * float(np_contigency[1, 1])) / ((float_probability_A - float_probability_AB) * (float_probability_B - float_probability_AB) + float_probability_AB * (1 + float_probability_AB - float_probability_A - float_probability_B))
            float_probability_AB_estimateByEM = float_num_AB_estimateByEM / (int_num_subject * 2)
            if abs(float_probability_AB_estimateByEM - float_probability_AB) < 0.0000001:
                break
            else:
                float_probability_AB = float_probability_AB_estimateByEM
        
        ### calculate D
        float_D = float_probability_AB - float_probability_A * float_probability_B
        ### calculate D prime
        if float_D >= 0:
            float_D_min = min([float_probability_A * (1 - float_probability_B), (1 - float_probability_A) * float_probability_B])
        else:
            float_D_min = max([-float_probability_A * float_probability_B, -(1 - float_probability_A) * (1 - float_probability_B)])
        float_D_prime = float_D / float_D_min
        ### calculate R square
        float_R_square = (float_D**2) / (float_probability_A * float_probability_a * float_probability_B * float_probability_b)
        
        return float_D_prime, float_R_square
    
    except ZeroDivisionError:
        return 1.0, 1.0

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def EstimateLDBlock(str_inputFileName_genotype, str_outputFilePath = "", float_threshold_DPrime = 0.8, float_threshold_RSquare = 0.8):
    """

    A function for implementing linkage disequilibrium (LD) dimension reduction. In genotype data, a variant often exhibits high dependency with its nearby variants because of LD. In the practical implantation, we prefer to group these dependent features to reduce the dimension of features. In other words, we can take the advantages of LD to reduce the dimensionality of genetic features. In this regard, this function adopted the same approach developed by Lewontin (1964) to estimate LD. We used D’ and r2 as the criteria to group highly dependent genetic features as blocks. In each block, we chose the features with the largest minor allele frequency to represent other features in the same block.

    Args:
        str_inputFileName_genotype (str): File name of input genotype data
        str_outputFilePath (str): File path of output file
        float_threshold_DPrime (float): The Dprime threshold for discriminating a LD block (default: 0.8)
        float_threshold_RSquare (float): The RSquare threshold for discriminating a LD block (default: 0.8)

    Returns:
        - Expected Success Response::

            "step2: Estimate LD. DONE!"
    
    """
    
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype)
    
    ### get the number of snp
    int_num_snp = sum(1 for line in open(str_inputFileName_genotype))
    
    ### read .gen file and estimate the LD block
    list_outputLDBlock = []
    with open(str_inputFileName_genotype, "r") as file_inputFile:
        with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).replace(".gen", "_LDReduced.gen")), "w") as file_outputFile:
            ### create dictionary for LD block
            ### key: rsID; value:[minor allele requency, raw genotypes data]
            dict_thisLDBlock = {}
            ### put first snp into dictionary
            line_previousSnp = file_inputFile.readline()
            list_previousSnp = line_previousSnp.strip().split(" ")
            dict_thisLDBlock[list_previousSnp[1]] = [min(EstimateAlleleFrequency(line_previousSnp)), line_previousSnp]
            
            ### scan all other snps
            int_count_snp = 1
            for line in file_inputFile:
                list_thisSnp = line.strip().split(" ")
                
                ### estimate pairwise LD for all of the snps in dictionary
                bool_flag_inLD = True
                for key in dict_thisLDBlock.keys():
                    float_DPrime, float_RSquare = EstimatePairwiseLD(dict_thisLDBlock[key][1], line)
                    if float_DPrime < float_threshold_DPrime or float_RSquare < float_threshold_RSquare:
                        bool_flag_inLD = False
                        break
                
                ### if this snp not in this LD block, then output and clear the content of dictionary
                if bool_flag_inLD == False:
                    ### find a snp with maximum minor allele frequency to be representative snp
                    str_representative_rsid = list(dict_thisLDBlock.keys())[0]
                    for key in dict_thisLDBlock.keys():
                        if dict_thisLDBlock[key][0] > dict_thisLDBlock[str_representative_rsid][0]:
                            str_representative_rsid = key
                    list_outputLDBlock.append(str_representative_rsid + ":" + ",".join(dict_thisLDBlock.keys()))
                    file_outputFile.writelines(dict_thisLDBlock[str_representative_rsid][1])
                    dict_thisLDBlock.clear()
                ### add this snp to current dictionary
                dict_thisLDBlock[list_thisSnp[1]] = [min(EstimateAlleleFrequency(line)), line]
                
                ### show progress
                int_count_snp = int_count_snp + 1
                str_print = "step2: Processing: " + "{0:.2f}".format(float(int_count_snp) / int_num_snp * 100) + "%"
                sys.stdout.write('%s\r' % str_print)
                sys.stdout.flush()
            
            ### output the final LD block in dictionary
            ### find a snp with maximum minor allele frequency to be representative snp
            str_representative_rsid = list(dict_thisLDBlock.keys())[0]
            for key in dict_thisLDBlock.keys():
                if dict_thisLDBlock[key][0] > dict_thisLDBlock[str_representative_rsid][0]:
                    str_representative_rsid = key
            list_outputLDBlock.append(str_representative_rsid + ":" + ",".join(dict_thisLDBlock.keys()))
            file_outputFile.writelines(dict_thisLDBlock[str_representative_rsid][1])
    
    ### output the file of LD block
    ### output file format: rsid_representative: rsid_1,rsid_2,rsid_3,...(the snps in the same LD block)
    with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).replace(".gen", ".LDBlock")), "w") as file_outputFile:
        for item in list_outputLDBlock:
            file_outputFile.writelines(item + "\n")
    
    print("step2: Estimate LD. DONE! \t\t\t\t")