# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester Chang
"""

import os
import genepi

### step1_downloadUCSCDB
genepi.DownloadUCSCDB(str_hgbuild="hg19")

### step2_estimateLD
genepi.EstimateLDBlock(os.path.dirname(os.path.abspath(__file__)) + "/sample.gen", float_threshold_DPrime=0.9, float_threshold_RSquare=0.9)