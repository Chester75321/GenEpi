# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester Chang
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
from .step1_downloadUCSCDB import DownloadUCSCDB
from .step2_estimateLD import EstimateLDBlock
from .step3_splitByGene import SplitByGene
from .step4_singleGeneEpistasis_Logistic import SingleGeneEpistasisLogistic
from .step4_singleGeneEpistasis_Logistic import BatchSingleGeneEpistasisLogistic