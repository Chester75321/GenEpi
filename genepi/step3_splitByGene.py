# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester (Yu-Chuan Chang)
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import os
import numpy as np

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def SplitMegaGene(list_snpsOnGene, int_window, int_step, str_outputFilePath, str_outputFileName):
    """

    In order to extract genetic features for a gene, this function used the start and end positions of each gene from the local UCSC database to split the genetic features. Then, generate the .GEN files for each gene in the folder named snpSubsets.

    Args:
        list_snpsOnGene (list): A list contains SNPs on a gene
        int_window (int): The size of the sliding window
        int_step (int): The step of the sliding window
        str_outputFilePath (str): File path of output file
        str_outputFileName (str): File name of output file

    Returns:
        None
    
    """
    
    ### write to gen file if this gene is not mega gene
    int_total_window = int((len(list_snpsOnGene)-int_window)/int_step)
    if int_total_window <= 0:
        with open(os.path.join(str_outputFilePath, str_outputFileName + ".gen"), "w") as file_outputFile:
            for item in list_snpsOnGene:
                file_outputFile.writelines(item)

    ### write gen file of each window on current gene (output file name: geneSymbol_numOfSNPOnGene@windowNum.gen)
    else: 
        for idx_w in range(int_total_window):
            with open(os.path.join(str_outputFilePath, str_outputFileName.split("_")[0] + "@" + str(idx_w) + "_" + str_outputFileName.split("_")[1] + ".gen"), "w") as file_outputFile:
                for item in list_snpsOnGene[int_step*idx_w:int_step*idx_w+int_window]:
                    file_outputFile.writelines(item)
        
        ### write reminder SNPs to gen file
        with open(os.path.join(str_outputFilePath, str_outputFileName.split("_")[0] + "@" + str(int_total_window) + "_" + str_outputFileName.split("_")[1] + ".gen"), "w") as file_outputFile:
            for item in list_snpsOnGene[int_step*int_total_window:]:
                file_outputFile.writelines(item)
    
    return

def SplitByGene(str_inputFileName_genotype, str_inputFileName_UCSCDB = os.path.dirname(os.path.abspath(__file__)) + "/UCSCGenomeDatabase.txt", str_outputFilePath = ""):
    """

    In order to extract genetic features for a gene, this function used the start and end positions of each gene from the local UCSC database to split the genetic features. Then, generate the .GEN files for each gene in the folder named snpSubsets.

    Args:
        str_inputFileName_genotype (str): File name of input genotype data
        str_inputFileName_UCSCDB (str): File name of input genome regions
        str_outputFilePath (str): File path of output file

    Returns:
        - Expected Success Response::

            "step3: Split by gene. DONE!"
    
    Warnings:
        "Warning of step3: .gen file should be sorted by chromosome and position"
    
    """
    
    print("Warning of step3: .gen file should be sorted by chromosome and position")

    int_window = 1000
    int_step = 200
    
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.join(os.path.dirname(str_inputFileName_genotype), "snpSubsets")
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    ### load UCSC Genome Database
    list_UCSCGenomeDatabase = []
    with open(str_inputFileName_UCSCDB, "r") as file_inputFile:
        for line in file_inputFile:
            list_UCSCGenomeDatabase.append(line.strip().split(","))
    np_UCSCGenomeDatabase = np.array(list_UCSCGenomeDatabase)
    
    ### scan all snp
    with open(str_inputFileName_genotype, "r") as file_inputFile:
        idx_gene = 0
        list_snpsOnGene = []
        for line in file_inputFile:
            ### get information of each snp
            list_thisSnp = line.strip().split(" ")
            int_chromosome = int(list_thisSnp[0])
            int_position = int(list_thisSnp[2])

            ### current gene is in next chromosome
            if int_chromosome < int(np_UCSCGenomeDatabase[idx_gene, 0]):
                continue
            ### current snp of genotype data is in next chromosome
            elif int_chromosome > int(np_UCSCGenomeDatabase[idx_gene, 0]):
                if len(list_snpsOnGene) != 0:
                    #### write gen file of current gene (output file name: geneSymbol_numOfSNPOnGene.gen)
                    #str_outputFileName = str(np_UCSCGenomeDatabase[idx_gene, 4]) + "_" + str(len(list_snpsOnGene)) + ".gen"
                    #with open(os.path.join(str_outputFilePath, str_outputFileName), "w") as file_outputFile:
                    #    for item in list_snpsOnGene:
                    #        file_outputFile.writelines(item)
                    str_outputFileName = str(np_UCSCGenomeDatabase[idx_gene, 4]) + "_" + str(len(list_snpsOnGene))
                    SplitMegaGene(list_snpsOnGene, int_window, int_step, str_outputFilePath, str_outputFileName)
                list_snpsOnGene = []
                while int_chromosome > int(np_UCSCGenomeDatabase[idx_gene, 0]):
                    ### jump to next gene
                    idx_gene = idx_gene + 1
                    ### if no next gene then break
                    if idx_gene == np_UCSCGenomeDatabase.shape[0]:
                        break
                    ### current snp on next gene
                    if int(np_UCSCGenomeDatabase[idx_gene, 1]) <= int_position and int_position <= int(np_UCSCGenomeDatabase[idx_gene, 2]) and int_chromosome == int(np_UCSCGenomeDatabase[idx_gene, 0]):
                        list_snpsOnGene.append(line)
    
            ### chromosome numbers of current snp and gene are match
            else:
                ### current snp on current gene
                if int(np_UCSCGenomeDatabase[idx_gene, 1]) <= int_position and int_position <= int(np_UCSCGenomeDatabase[idx_gene, 2]):
                    list_snpsOnGene.append(line)
                ### snp position exceed this gene
                elif int_position > int(np_UCSCGenomeDatabase[idx_gene, 2]):
                    if len(list_snpsOnGene) != 0:
                        #### write gen file of current gene (output file name: geneSymbol_numOfSNPOnGene.gen)
                        #str_outputFileName = str(np_UCSCGenomeDatabase[idx_gene, 4]) + "_" + str(len(list_snpsOnGene)) + ".gen"
                        #with open(os.path.join(str_outputFilePath, str_outputFileName), "w") as file_outputFile:
                        #    for item in list_snpsOnGene:
                        #        file_outputFile.writelines(item)
                        str_outputFileName = str(np_UCSCGenomeDatabase[idx_gene, 4]) + "_" + str(len(list_snpsOnGene))
                        SplitMegaGene(list_snpsOnGene, int_window, int_step, str_outputFilePath, str_outputFileName)
                    list_snpsOnGene = []
                    while int_position > int(np_UCSCGenomeDatabase[idx_gene, 2]) and int_chromosome == int(np_UCSCGenomeDatabase[idx_gene, 0]):
                        ### jump to next gene
                        idx_gene = idx_gene + 1
                        ### if no next gene then break
                        if idx_gene == np_UCSCGenomeDatabase.shape[0]:
                            break
                        ### snp on next gene
                        if int(np_UCSCGenomeDatabase[idx_gene, 1]) <= int_position and int_position <= int(np_UCSCGenomeDatabase[idx_gene, 2]) and int_chromosome == int(np_UCSCGenomeDatabase[idx_gene, 0]):
                            list_snpsOnGene.append(line)

            ### if the index of gene out of the boundary of DB then break
            if idx_gene >= np_UCSCGenomeDatabase.shape[0]:
                break
    
    print("step3: Split by gene. DONE!")
