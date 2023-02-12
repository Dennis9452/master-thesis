# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 14:00:40 2017

@author: user
"""



def Concatenating_gene_expression_profile(USER_FOLDER, Case_ID, mRNA_file_name, lncRNA_file_name):
    import os,numpy
    import pandas as pd

    #PATH = os.path(UPLOAD_FOLDER)
    
    mRNA_file = pd.read_table(os.path.join(USER_FOLDER, mRNA_file_name), sep='\t', header=0, index_col=0, encoding='utf-8')
    lncRNA_file = pd.read_table(os.path.join(USER_FOLDER, lncRNA_file_name), sep='\t', header=0, index_col=0, encoding='utf-8')
    
    sample_list = sorted(list(set(mRNA_file.columns.values.tolist())&set(lncRNA_file.columns.values.tolist())))
    if sample_list == 0:
        return 1
    
    frames = [mRNA_file.loc[:,sample_list], lncRNA_file.loc[:,sample_list]]
    result = pd.concat(frames)
    
    fname = 'Gene_Expression_profile_raw_data.txt'
    result.to_csv(os.path.join(USER_FOLDER, Case_ID+'_'+fname), sep='\t', encoding='utf-8')
 
    return 0
