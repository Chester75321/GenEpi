# -*- coding: utf-8 -*-
"""
Created on Apr 2019

@author: Chester (Yu-Chuan Chang)
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import argparse
import os
import genepi
import numpy as np

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def ArgumentsParser():
    ### define arguments
    str_description = ''
    'This script is a tool for extracting single gene sequence from whole genome reference'
    parser = argparse.ArgumentParser(prog='creategeneref', description=str_description)
    
    ### define arguments for I/O
    parser.add_argument("-i", required=True, help="filename of the input txt file")
    parser.add_argument("-r", required=True, help="filename of the reference genome")    
    parser.add_argument("-o", required=False, help="output file path")
    
    ### define other optional arguments    
    parser.add_argument("-b", required=False, default="hg19", choices=["hg19", "hg38"], help="human genome build")
    
    return parser

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def main(args=None):
    ### obtain arguments from argument parser
    args = ArgumentsParser().parse_args(args)

    ### get arguments for I/O
    str_inputFileName_in = args.i
    str_inputFileName_ref = args.r
    str_outputFilePath = ""
    if args.o is not None:
        str_outputFilePath = args.o
    else:
        str_outputFilePath = os.path.dirname(str_inputFileName_ref)
        
    ### get gene symbol from input file
    list_gene = []
    with open(str_inputFileName_in, "r") as file_inputFile:
        for line in file_inputFile:
            list_gene.append(line.strip())
    
    ### switch genome build
    if args.b == "hg38":
        genepi.DownloadUCSCDB(str_hgbuild=args.b)
    
    ### load UCSC Genome Database
    list_UCSCGenomeDatabase = []
    with open(os.path.join(os.path.dirname(genepi.__file__), "UCSCGenomeDatabase.txt"), "r") as file_inputFile:
        for line in file_inputFile:
            list_UCSCGenomeDatabase.append(line.strip().split(","))
    np_UCSCGenomeDatabase = np.array(list_UCSCGenomeDatabase)
    
    ### get start & end of gene
    np_gene = np_UCSCGenomeDatabase[np.where(np_UCSCGenomeDatabase[:, 4] == list_gene[0]), :].ravel()
    if np_gene[3] == "+":
        np_gene[1] = int(np_gene[1]) + 1000
    else:
        np_gene[2] = int(np_gene[2]) - 1000
    
    ### output bed file
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp.bed"), "w") as file_outputFile:
        file_outputFile.writelines("chr" + str(np_gene[0]).replace("chr", "") + "\t" + str(np_gene[1]) + "\t" + str(np_gene[2]))
    
    ### use bedtools to extract single gene sequence
    str_command = "bedtools getfasta -fi " + str_inputFileName_ref + " -bed " + os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp.bed")
    list_fasta = str(os.popen(str_command).read()).strip().split("\n")
    list_fasta[0] = ">" + str(list_gene[0])
    
    ### output single gene fasta
    with open(os.path.join(str_outputFilePath, str(np_gene[4]) + ".fasta"), "w") as file_outputFile:
        for line in list_fasta:
            file_outputFile.writelines(line + "\n")
    
    ### use samtools to generate index
    str_command = "samtools faidx " + os.path.join(str_outputFilePath, str(np_gene[4]) + ".fasta")
    os.system(str_command)
    
    ### use GATK to generate dict
    str_command = "python3 /opt/gatk/gatk CreateSequenceDictionary -R " + os.path.join(str_outputFilePath, str(np_gene[4]) + ".fasta") + " -O " + os.path.join(str_outputFilePath, str(np_gene[4]) + ".dict")
    os.system(str_command)
    
if __name__ == "__main__":
    main()