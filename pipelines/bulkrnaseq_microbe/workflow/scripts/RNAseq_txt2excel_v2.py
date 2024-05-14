
# -*- coding: utf-8 -*-
# Created on Tue Jun  5 15:56:48 2018
# @author: simone oberhaensli



import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--de', dest='de', help="One or more tables with DESeq2 results as produced by pipeline", nargs="*")
args = parser.parse_args()


def txt2exc(file):
    #-- read DESeq output file into dataframe --#
    df=pd.read_table(file,  sep='\t', index_col=None)
	# insert an empty column to have a space between DE results and normalised counts
    df['spacer'] = np.nan # remove if spacer is not needed
    df=df.fillna('') # remove if spacer is not needed
    
    #-- Do cosmetics --#
    # extract column names
    cols = list(df)
    
    # round *.norm columns
    roundlist=[s for s in cols if '.norm' in s]
    for index in roundlist:
        df[[index]]=df[[index]].apply(lambda x: pd.Series.round(x, 2))
    
    # remove unnecessary columns
    cols.remove('heatmap.id')
    cols.remove('baseMean')

    # identify indices of columns with *rlog (number depends on number of replicates)
    indicesRlog = [i for i, s in enumerate(cols) if 'rlog' in s]
    # identify indices of columns with raw counts (number depends on number of replicates)
    NrSamples=len(indicesRlog)
    # indices of raw counts plus base.mean
    indicesSamples = range(6,7+NrSamples)
    # combine indices of unnecessary columns
    allInd=indicesRlog+indicesSamples
    # delete all unnecessary column names from list
    for index in sorted(allInd, reverse=True):
        del cols[index]
    # move spacer to separate normalized counts from the rest    
    cols.insert(6, cols.pop(cols.index('spacer'))) # remove if spacer is not needed
    # move some column names to the front and insert them in front of spacer and norm. counts
    l=['log2FoldChange','lfcSE','stat','pvalue','padj']
    c=5
    for i in l:
        c=c+1
        cols.insert(c, cols.pop(cols.index(i)))
    
    
    #--sort df according to the labels in the list "cols"--#
    dff = df.loc[:, cols]
    
    #-- export dataframe to excel --##    
    # define file name
    filename=str(file)[:-18] 
    filename=filename + '.xlsx'
    writer = pd.ExcelWriter(filename, engine='xlsxwriter')
    dff.to_excel(writer, sheet_name='Sheet1',index=False) # change sheet name
    writer.save()
   
    
    


#-- Read in de files and process them --#
for currentFile in args.de:  
    txt2exc(currentFile)    

