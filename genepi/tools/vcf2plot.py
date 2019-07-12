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
import numpy as np
import scipy.stats as sp
import pandas as pd
import matplotlib.pyplot as plt

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def ArgumentsParser():
    ### define arguments
    str_description = ''
    'This script is a visualization tool for displaying the potential rare variants in vcf file'
    parser = argparse.ArgumentParser(prog='vcf2plot', description=str_description)
    
    ### define arguments for I/O
    parser.add_argument("-v", required=True, help="filename of the input vcf file")
    parser.add_argument("-p", required=True, help="filename of the input phenotype")    
    parser.add_argument("-o", required=False, help="output file path")
    
    ### define other optional arguments
    parser.add_argument("-m", required=False, default="r", choices=["c", "r"], help="choose phenotype type: c for classification; r for regression")
    parser.add_argument("-t", required=False, default=0.0, help="the threshold for categorize")
    
    return parser

def DrawRegressionPlot(df_variant, df_data, float_threshold, str_outputFilePath):
    ### draw scatter plot for real number phenotype
    plt.figure(figsize=(15, 40))
    plt.scatter(df_data["phenotype"], df_data["order"], color=df_data["color"], s=40, marker="o", alpha=0.5)
    plt.ylim(-1 , 100)
    plt.yticks(df_data["order"], df_data["variant_ID_x"])
    for idx_y in df_data["order"]:
        plt.axhline(y=idx_y, color="#CCCCCC", linewidth=1)
    plt.axvline(x=float_threshold, color="#CCCCCC", linewidth=1)
    plt.savefig(os.path.join(str_outputFilePath, "vcf2plot.jpg"), bbox_inches='tight', dpi=100)

def DrawClassificationPlot(df_variant, df_data, str_outputFilePath):
    ### count the number of mutations in each variant
    df_data = df_data.groupby(["order", "variant_ID_x", "phenotype"])["allele_type"].sum().unstack(fill_value=0).reset_index()
    ### draw barplot plot for dualistic phenotype
    float_barwidth = 0.35
    float_opacity = 0.8
    fig, ax = plt.subplots(1, 1, figsize=(15, 40))
    plt.barh(df_data["order"] + float_barwidth * 0.5, df_data.iloc[:,-1], float_barwidth, alpha=float_opacity, color="#CC2900", label='Case')
    plt.barh(df_data["order"] + float_barwidth * 1.5, df_data.iloc[:,-2], float_barwidth, alpha=float_opacity, color="#29A329", label='Control')
    plt.ylim(-1 , 100)
    plt.yticks(df_data["order"] + float_barwidth, df_data["variant_ID_x"])
    for idx_y in df_data["order"]:
        plt.axhline(y=idx_y + float_barwidth * 1, color="#CCCCCC", linewidth=1)
    plt.legend() 
    plt.savefig(os.path.join(str_outputFilePath, "vcf2plot.jpg"), bbox_inches='tight', dpi=100)

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def main(args=None):
    ### obtain arguments from argument parser
    args = ArgumentsParser().parse_args(args)

    ### get arguments for I/O
    str_inputFileName_vcf = args.v
    str_inputFileName_phenotype = args.p
    str_outputFilePath = ""
    if args.o is not None:
        str_outputFilePath = args.o
    else:
        str_outputFilePath = os.path.dirname(str_inputFileName_vcf)
    float_threshold = float(args.t)    

    ### get phenotype of each sample
    dict_pheno = {}
    with open(str_inputFileName_phenotype, "r") as file_inputFile:
        ### skip two headers
        file_inputFile.readline()
        for line in file_inputFile:
            list_line = line.strip().split("\t")
            list_line[1] = "{:010d}".format(int(list_line[1]))
            dict_pheno[list_line[1]] = list_line[-1]
    
    ### scan vcf file
    # [idx, variant ID, odds ratio, p-value, mutation in control, mutation in case]
    np_variant = []
    # [idx, variant ID, sample ID, allele type, phenotype]
    np_data = []
    idx_variant = 0
    with open(str_inputFileName_vcf, "r") as file_inputFile:
        for line in file_inputFile:
            ### skip headers
            if line[0] == "#":
                list_sample_id = line.strip().split("\t")
            else:
                ### create contigency table for fisher-exact test
                np_contingency = np.array([[0, 0], [0, 0]])
                list_line = line.strip().split("\t")
                ### skip low quality variants
                if list_line[6] != "LowQual":
                    ### grep genotype of each sample
                    for idx_sample in range(len(list_line[9:])):
                        str_sample = list_sample_id[idx_sample + 9].split("_")[1]
                        if str_sample in dict_pheno.keys():
                            list_item = list_line[idx_sample + 9].split(":")
                            ### if the sample be classified as class I
                            if float(dict_pheno[str_sample]) > float_threshold:
                                if list_item[0] == ".":
                                    np_contingency[0, 1] += 2
                                elif list_item[0] == "0/1":
                                    np_contingency[0, 1] += 1
                                    np_contingency[1, 1] += 1
                                    np_data.append([idx_variant, str(list_line[0] + ":" + list_line[1]), str_sample, 1, dict_pheno[str_sample]])                        
                                elif list_item[0] == "1/1":
                                    np_contingency[1, 1] += 2
                                    np_data.append([idx_variant, str(list_line[0] + ":" + list_line[1]), str_sample, 2, dict_pheno[str_sample]])
                                else:
                                    np_contingency[0, 1] += 2
                            ### the sample is belong class II
                            else:
                                if list_item[0] == ".":
                                    np_contingency[0, 0] += 2
                                elif list_item[0] == "0/1":
                                    np_contingency[0, 0] += 1
                                    np_contingency[1, 0] += 1
                                    np_data.append([idx_variant, str(list_line[0] + ":" + list_line[1]), str_sample, 1, dict_pheno[str_sample]])                        
                                elif list_item[0] == "1/1":
                                    np_contingency[1, 0] += 2
                                    np_data.append([idx_variant, str(list_line[0] + ":" + list_line[1]), str_sample, 2, dict_pheno[str_sample]])
                                else:
                                    np_contingency[0, 0] += 2
                    ### execute fisher-exact test
                    oddsratio, pvalue = sp.fisher_exact(np_contingency)
                    np_variant.append([idx_variant, str(list_line[0] + ":" + list_line[1]), oddsratio, pvalue, np_contingency[1, 0], np_contingency[1, 1]])
                    idx_variant += 1
    
    ### prepare variant's info for plot
    df_variant = pd.DataFrame(np_variant)
    df_variant.columns = ["idx", "variant_ID", "odds_ratio", "p-value", "mutation_in_control", "mutation_in_case"]
    df_variant = df_variant.assign(mutation_sum = df_variant["mutation_in_control"] + df_variant["mutation_in_case"])
    df_variant = df_variant.sort_values(by=["odds_ratio", "mutation_sum"])
    df_variant = df_variant.dropna(subset=["odds_ratio"])
    ### obtain top 100 significant variants
    if df_variant.shape[0] > 100:
        df_temp = df_variant[df_variant["odds_ratio"] < 1].sort_values(by=["odds_ratio", "mutation_sum"], ascending=[True, False])
        df_temp = pd.concat([df_temp, df_variant[df_variant["odds_ratio"] > 1]])
        df_variant = df_temp
        df_temp = df_variant[:50]
        df_temp = pd.concat([df_temp, df_variant[-50:]])
        df_variant = df_temp
    ### set order of variants for plot
    df_variant = df_variant.assign(order = range(df_variant.shape[0]))
    df_variant["label"]= df_variant["odds_ratio"].map('{:,.4f}'.format) + "-" + df_variant["mutation_in_control"].map(str) + ":" + df_variant["mutation_in_case"].map(str)
    
    ### prepare mutations's info for plot
    df_data = pd.DataFrame(np_data)
    df_data.columns = ["idx", "variant_ID", "sample_ID", "allele_type", "phenotype"]
    df_data["phenotype"] = df_data["phenotype"].astype("float")
    df_data = pd.merge(df_data, df_variant, how="inner", on="idx")
    ### set color of mutations for plot
    df_data["color"] = np.where(df_data["phenotype"] > float_threshold, "#CC2900", "#29A329")
    
    ### draw plot
    if args.m=="c":
        DrawClassificationPlot(df_variant, df_data, str_outputFilePath)
    else:
        DrawRegressionPlot(df_variant, df_data, float_threshold, str_outputFilePath)

if __name__ == "__main__":
    main()