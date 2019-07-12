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
import sys

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def ArgumentsParser():
    ### define arguments
    str_description = ''
    'This script is a preprossing tool of GenEpi for converting vcf file to oxford format'
    parser = argparse.ArgumentParser(prog='vcf2oxford', description=str_description)
    
    ### define arguments for I/O
    parser.add_argument("-v", required=True, help="filename of the input vcf file")
    parser.add_argument("-o", required=False, help="output file path")
    
    return parser

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def main(args=None):
    ### obtain arguments from argument parser
    args = ArgumentsParser().parse_args(args)

    ### get arguments for I/O
    str_inputFileName_vcf = args.v
    str_outputFilePath = ""
    if args.o is not None:
        str_outputFilePath = args.o
    else:
        str_outputFilePath = os.path.dirname(str_inputFileName_vcf)

    ### count lines of input files
    int_num_vcf = sum(1 for line in open(str_inputFileName_vcf))

    ### scan .vcf file
    list_sample_id = []
    list_gen_all = []
    int_count_vcf = 0
    with open(str_inputFileName_vcf, "r") as file_inputFile:
        for line in file_inputFile:
            int_count_vcf += 1
            ### skip headers
            if line[0] == "#":
                list_sample_id = line.strip().split("\t")
            else:
                list_line = line.strip().split("\t")
                ### skip low quality variants
                if list_line[6] != "LowQual":
                    list_gen = []
                    ### chromosome
                    list_gen.append(list_line[0].replace("chr", ""))
                    ### rsID: chr|position|alt|ref
                    if list_line[2] == ".":
                        list_line[2] = list_line[0].replace("chr", "") + "|" + list_line[1] + "|" + list_line[4] + "|" + list_line[3]
                    list_gen.append(list_line[2])
                    ### position
                    list_gen.append(list_line[1])
                    ### alternative allele
                    list_gen.append(list_line[4])
                    ### reference allele
                    list_gen.append(list_line[3])
                    
                    ### encode
                    for item in list_line[9:]:
                        list_item = item.split(":")
                        if list_item[0] == ".":
                            list_gen = list_gen + ["0", "0", "1"]
                        elif list_item[0] == "0/1":
                            list_gen = list_gen + ["0", "1", "0"]
                        elif list_item[0] == "1/1":
                            list_gen = list_gen + ["1", "0", "0"]
                        else:
                            list_gen = list_gen + ["0", "0", "0"]
                    list_gen_all.append(list_gen)
            str_print = "Preprocessing vcf2oxford: " + "{0:.2f}".format(float(int_count_vcf) / int_num_vcf * 100) + "%" + "\t\t"
            sys.stdout.write('%s\r' % str_print)
            sys.stdout.flush()

    ### ouput sample file
    list_sample_id = list_sample_id[9:]
    with open(os.path.join(str_outputFilePath, str_inputFileName_vcf.replace(".vcf", ".sample")), "w") as file_outputFile:
        file_outputFile.writelines("ID_1 ID_2 missing sex phenotype" + "\n")
        file_outputFile.writelines("0 0 0 D B" + "\n")
        for item in list_sample_id:    
            file_outputFile.writelines(item + " " + item + " " + "NA 0 NA" + "\n")

    ### output gen file        
    with open(os.path.join(str_outputFilePath, str_inputFileName_vcf.replace(".vcf", ".gen")), "w") as file_outputFile:
        for item in list_gen_all:
            file_outputFile.writelines(" ".join(item) + "\n")

if __name__ == "__main__":
    main()