# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester (Yu-Chuan Chang)
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import os
import pymysql

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
# database schema: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
# human genome build could be: hg18, hg19, hg38, etc.
def DownloadUCSCDB(str_outputFilePath = os.path.dirname(os.path.abspath(__file__)), str_hgbuild = "hg19"):
    """

    To obtain the gene information such as official gene symbols and genomic coordinates, this function is for retrieving kgXref and knownGene data table from the UCSC human genome annotation database

    Args:
        str_outputFilePath (str): File path of output database
        str_hgbuild (str): Genome build (eg. "hg19")

    Returns:
        - Expected Success Response::

            "step1: Down load UCSC Database. DONE!"
    
    """

    ### create connection
    conv = {pymysql.constants.FIELD_TYPE.LONG: int}
    conn = pymysql.Connect(host = "genome-mysql.cse.ucsc.edu", user = "genome", passwd = "",db = str_hgbuild,  local_infile = 1, conv = conv)
    
    ### execute sql command
    str_sqlCommand = "SELECT chr, CASE WHEN strand='+' THEN txStart-1000 ELSE txStart END AS txStart, CASE WHEN strand='-' THEN txEnd+1000 ELSE txEnd END AS txEnd, strand, geneSymbol FROM "
    str_sqlCommand = str_sqlCommand + "(SELECT REPLACE(chr, 'chr', '') AS chr, txStart, txEnd, strand, geneSymbol, MAX(ABS(txEnd-txStart)) FROM ( "
    str_sqlCommand = str_sqlCommand + "SELECT knownGene.chrom AS chr, knownGene.txStart AS txStart, knownGene.txEnd AS txEnd, knownGene.strand AS strand, kgXref.geneSymbol AS geneSymbol FROM kgXref INNER JOIN knownGene ON kgXref.kgID=knownGene.name WHERE LEFT(kgXref.mRNA, 2) IN ('NR', 'NM')) AS L1 "
    str_sqlCommand = str_sqlCommand + "GROUP BY geneSymbol) AS L2 WHERE LEFT(chr, 1) NOT IN ('X', 'Y', 'M', 'U') AND chr NOT LIKE '%\_%' ORDER BY CAST(chr AS UNSIGNED), txStart"
    cur = conn.cursor()
    cur.execute(str_sqlCommand)
    db_result = cur.fetchall()
    cur.close()
    conn.close()

    ### output database
    with open(str_outputFilePath + "/UCSCGenomeDatabase.txt", "w") as file_outputFile:
        for item in db_result:
            file_outputFile.writelines(",".join(item) + "\n")
    
    print("step1: Down load UCSC Database. DONE!")
