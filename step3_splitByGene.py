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
def SplitByGene(str_inputFileName_genotype, str_inputFileName_UCSCDB = os.path.dirname(os.path.abspath(__file__)) + "/UCSCGenomeDatabase.txt", str_outputFilePath = ""):
    print("Warning of step3: .gen file should be sorted by chromosome and position")
    
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
                    ### write gen file of current gene (output file name: geneSymbol_numOfSNPOnGene.gen)
                    str_outputFileName = str(np_UCSCGenomeDatabase[idx_gene, 4]) + "_" + str(len(list_snpsOnGene)) + ".gen"
                    with open(os.path.join(str_outputFilePath, str_outputFileName), "w") as file_outputFile:
                        for item in list_snpsOnGene:
                            file_outputFile.writelines(item)
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
                        ### write gen file of current gene (output file name: geneSymbol_numOfSNPOnGene.gen)
                        str_outputFileName = str(np_UCSCGenomeDatabase[idx_gene, 4]) + "_" + str(len(list_snpsOnGene)) + ".gen"
                        with open(os.path.join(str_outputFilePath, str_outputFileName), "w") as file_outputFile:
                            for item in list_snpsOnGene:
                                file_outputFile.writelines(item)
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