# -*- coding: utf-8 -*-
"""
Created on Feb 2019

@author: Chester (Yu-Chuan Chang)
"""

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import sys
import os
import subprocess

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def CorrectMissingID(str_inputFileName_genotype, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype)

    ### count lines of input files
    int_num_genotype = sum(1 for line in open(str_inputFileName_genotype))
    
    ### read .gen file
    with open(str_inputFileName_genotype, "r") as file_inputFile:
        with open(os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).replace(".gen", "_addID.gen")), "w") as file_outputFile:
            int_count_snp = 0
            for line in file_inputFile:
                list_thisSnp = line.strip().split(" ")
                ### if current snp have no ID, assign temporary ID (chromosome + position) to it
                if list_thisSnp[1] == ".":
                    list_thisSnp[1] = 'C' + str(list_thisSnp[0]) + 'P' + str(list_thisSnp[2])
                    file_outputFile.writelines(' '.join(list_thisSnp) + "\n")
                else:
                    file_outputFile.writelines(line)
                
                ### show progress
                int_count_snp = int_count_snp + 1
                str_print = "step0: Processing: " + "{0:.2f}".format(float(int_count_snp) / int_num_genotype * 100) + "\t\t"
                sys.stdout.write('%s\r' % str_print)
                sys.stdout.flush()
    
    print("step0: Correct Missingn ID DONE! \t\t\t\t")

def CorrectMissingIDAWK(str_inputFileName_genotype, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype)
    
    str_command = "awk -F \" \" '{if ($2==\".\") $2=\"C\"$1\"P\"$3; print}' "
    str_command += str_inputFileName_genotype + " > "
    str_command += os.path.join(str_outputFilePath, os.path.basename(str_inputFileName_genotype).replace(".gen", "_addID.gen"))
    subprocess.call(str_command, shell=True)

if __name__ == '__main__':
    CorrectMissingIDAWK(sys.argv[1], sys.argv[2])