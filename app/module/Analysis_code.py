# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 18:50:23 2018

@author: user
"""
import matplotlib
matplotlib.use('Agg')
import pandas as pd, numpy as np, matplotlib.font_manager as fm, matplotlib.pyplot as plt, seaborn as sns
from itertools import product, combinations
from scipy.stats import spearmanr
from collections import Counter
from scipy.stats import linregress
from matplotlib import gridspec
import sys

def correlation(PATH, Case_ID, mRNA_file_name, lncRNA_file_name, miRNA_file_name):
    cutoff = np.log2(0.7)
    for (rna, rna_file_name) in [('mRNA', mRNA_file_name), ('miRNA', miRNA_file_name)]:
        data1 = pd.read_table('/'.join([PATH, lncRNA_file_name]), sep='\t', header=0, index_col=0, encoding='utf-8')
#        data1 = np.log2(data1)
        Median = data1.quantile(q=0.5, axis='columns')
        data1 = data1.loc[Median.index.values[np.where(Median > cutoff)[0]].tolist(), :]
        data2 = pd.read_table('/'.join([PATH, rna_file_name]), sep='\t', header=0, index_col=0, encoding='utf-8')
#        data2 = np.log2(data2)
        if 'mRNA' in rna:
            Median = data2.quantile(q=0.5, axis='columns')
            data2 = data2.loc[Median.index.values[np.where(Median > cutoff)[0]].tolist(), :]
        del Median
        data = pd.concat([data1, data2], axis='index', keys=['lncRNA', rna], names=['group_key', 'ID'])
        del data1, data2
        data_index = data.index.get_level_values('ID').tolist()
        n = data.shape[1]
        idx = data.index.get_level_values('group_key') == 'lncRNA'
        lncRNA_idx = idx.nonzero()[0]
        idx = data.index.get_level_values('group_key') == rna
        rna_idx = idx.nonzero()[0]
        del idx
        pairs_idx = pd.DataFrame([[x, y] for x, y in product(lncRNA_idx, rna_idx)], columns=['lncRNA_idx', '_'.join([rna, 'idx'])])
        pairs_ID = pd.DataFrame([[data_index[x], data_index[y]] for x, y in product(lncRNA_idx, rna_idx)], columns=['lncRNA_ID', '_'.join([rna, 'ID'])])
        del data_index, lncRNA_idx, rna_idx
        fname = '_'.join([Case_ID, '2'.join(['lncRNA', rna]), 'ID.txt'])
        pairs_ID.to_csv(str('/'.join([PATH, fname])), sep='\t', na_rep='NA', index=False, encoding='utf-8')
        del pairs_ID
        correlation = spearmanr(data, axis=1)[0]
        del data
        rho = pd.DataFrame(correlation[pairs_idx['lncRNA_idx'], pairs_idx['_'.join([rna, 'idx'])]], columns=['SCC'])
        del correlation
        rho['z-score'] = np.sqrt((n-3)/1.06)*(0.5*np.log((1+rho['SCC'])/(1-rho['SCC'])))
        fname = '_'.join([Case_ID, '2'.join(['lncRNA', rna]), 'correlation.txt'])
        rho.to_csv(str('/'.join([PATH, fname])), sep='\t', na_rep='NA', index=False, encoding='utf-8')
        del n, pairs_idx, rho
        
    print('correlation finish')

def bipartite_network(PATH, Case_ID):
    for (rna_x, rna_y) in [('lncRNA', 'mRNA'), ('lncRNA', 'miRNA')]:
       fname = '_'.join([Case_ID, '2'.join([rna_x, rna_y]), 'correlation.txt'])
       data = pd.read_table('/'.join([PATH, fname]), sep='\t', header=0, index_col=None, usecols=['z-score'], encoding='utf-8')
       fname = '_'.join([Case_ID, '2'.join([rna_x, rna_y]), 'ID.txt'])
       try:
           ID_list = pd.read_table('/'.join([PATH, fname]), sep='\t', header=0, index_col=None, encoding='utf-8') 
       except Exception as error:
           print (str(error))
           return str(error)
       for threshold in np.arange(2, 10.1, 0.25):
            Interaction_nodeX = {}
            Interaction_nodeY = {}
            Edges = 0
            for [nodeX, nodeY] in ID_list.loc[data.index.values[np.where(data['z-score'] > threshold)[0]].tolist(), :].values:
                if nodeX in Interaction_nodeX:
                    Interaction_nodeX[nodeX] += 1
                else:
                    Interaction_nodeX[nodeX] = 1
                if nodeY in Interaction_nodeY:
                    Interaction_nodeY[nodeY] += 1
                else:
                    Interaction_nodeY[nodeY] = 1
                Edges += 1
            Degree = np.array(list(Interaction_nodeX.values())+list(Interaction_nodeY.values()))
            n = len(Degree)
            degree, probability = [], []
            for key, val in Counter(Degree).items():
                degree.append(np.log10(key))
                probability.append(np.log10(val/n))
            r_squared = linregress(degree, probability)[2]**2
            del Interaction_nodeX, Interaction_nodeY, Edges, nodeX, nodeY, Degree, n, degree, probability
            if r_squared > 0.8:
                break
       temp = ID_list.loc[data.index.values[np.where(data['z-score'] > threshold)[0]].tolist(), :]
       del data, ID_list
       fname = '_'.join([Case_ID, '2'.join([rna_x, rna_y]), 'bipartite_network.txt'])
       temp.to_csv(str('/'.join([PATH, fname])), sep='\t', na_rep='NA', index=False, encoding='utf-8')
       del fname, temp
    print('biparite finish')

def calculate_association_index(Interaction_A, Interaction_B, N_Y_nodes):
    setA = set(Interaction_A)
    setB = set(Interaction_B)
    NabIntersection = np.abs(len(setA&setB))
    NabUnion = np.abs(len(setA|setB))
    Ny = N_Y_nodes
    Na = np.abs(len(setA))
    Nb = np.abs(len(setB))
    minab = min(Na,Nb)
    Jaccard = NabIntersection/NabUnion
    Simpson = NabIntersection/minab
    Geometric = (NabIntersection**2)/(Na*Nb)
    Cosine = NabIntersection/((Na*Nb)**0.5)
    PCC = ((NabIntersection*Ny)-(Na*Nb))/((Na*Nb*(Ny-Na)*(Ny-Nb))**0.5)
    s = '-'
    if PCC >= 0:
        s = 'p'
    else:
        s = 'n'
    
    return {'Jaccard':Jaccard, 'Simpson':Simpson, 'Geometric':Geometric, 'Cosine':Cosine, 'PCC':np.abs(PCC), 'positive/negative':s}
    
def association_indices(PATH, Case_ID):
    Bipartite_Node = {'lncRNA':('mRNA', 'miRNA')}
    for typeX, (typeY, typeZ) in Bipartite_Node.items():
        networkY = '2'.join([typeX, typeY])
        Interaction_XY = {}
        temp = []
        fname = '_'.join([Case_ID, networkY, 'bipartite_network'])+'.txt'
        with open('/'.join([PATH, fname]), mode='r', encoding='utf-8') as fr:
            next(fr)
            for line in fr:
                X = line.strip().split('\t')[0]
                Y = line.strip().split('\t')[1]
                if X not in Interaction_XY:
                    Interaction_XY.update({X:[Y]})
                else:
                    Interaction_XY[X].append(Y)
                temp.append(Y)
        Total_typeY_nodes = len(set(temp))
        networkZ = '2'.join([typeX, typeZ])
        Interaction_XZ = {}
        temp = []
        fname = '_'.join([Case_ID, networkZ, 'bipartite_network'])+'.txt'
        with open('/'.join([PATH, fname]), mode='r', encoding='utf-8') as fr:
            next(fr)
            for line in fr:
                X = line.strip().split('\t')[0]
                Z = line.strip().split('\t')[1]
                if X not in Interaction_XZ:
                    Interaction_XZ.update({X:[Z]})
                else:
                    Interaction_XZ[X].append(Z)
                temp.append(Z)
        Total_typeZ_nodes = len(set(temp))
        temp = None
        list_typeX_nodes = list(set(Interaction_XY.keys())&set(Interaction_XZ.keys()))
        fname_Jaccard_Y = '_'.join([Case_ID, networkY, 'association_index_Jaccard'])+'.txt'
        fname_Simpson_Y = '_'.join([Case_ID, networkY, 'association_index_Simpson'])+'.txt'
        fname_Geometric_Y = '_'.join([Case_ID, networkY, 'association_index_Geometric'])+'.txt'
        fname_Cosine_Y = '_'.join([Case_ID, networkY, 'association_index_Cosine'])+'.txt'
        fname_PCC_Y = '_'.join([Case_ID, networkY, 'association_index_PCC'])+'.txt'
        fname_Jaccard_Z = '_'.join([Case_ID, networkZ, 'association_index_Jaccard'])+'.txt'
        fname_Simpson_Z = '_'.join([Case_ID, networkZ, 'association_index_Simpson'])+'.txt'
        fname_Geometric_Z = '_'.join([Case_ID, networkZ, 'association_index_Geometric'])+'.txt'
        fname_Cosine_Z = '_'.join([Case_ID, networkZ, 'association_index_Cosine'])+'.txt'
        fname_PCC_Z = '_'.join([Case_ID, networkZ, 'association_index_PCC'])+'.txt'
        with open('/'.join([PATH, fname_Jaccard_Y]), mode='w', encoding='utf-8') as fw_Jaccard_Y, open('/'.join([PATH, fname_Simpson_Y]), mode='w', encoding='utf-8') as fw_Simpson_Y, open('/'.join([PATH, fname_Geometric_Y]), mode='w', encoding='utf-8') as fw_Geometric_Y, open('/'.join([PATH, fname_Cosine_Y]), mode='w', encoding='utf-8') as fw_Cosine_Y, open('/'.join([PATH, fname_PCC_Y]), mode='w', encoding='utf-8') as fw_PCC_Y, open('/'.join([PATH, fname_Jaccard_Z]), mode='w', encoding='utf-8') as fw_Jaccard_Z, open('/'.join([PATH, fname_Simpson_Z]), mode='w', encoding='utf-8') as fw_Simpson_Z, open('/'.join([PATH, fname_Geometric_Z]), mode='w', encoding='utf-8') as fw_Geometric_Z, open('/'.join([PATH, fname_Cosine_Z]), mode='w', encoding='utf-8') as fw_Cosine_Z, open('/'.join([PATH, fname_PCC_Z]), mode='w', encoding='utf-8') as fw_PCC_Z:
            fw_Jaccard_Y.write('%s\t%s\t%s\n'%('nodeA', 'nodeB', 'Jaccard'))
            fw_Simpson_Y.write('%s\t%s\t%s\n'%('nodeA', 'nodeB', 'Simpson'))
            fw_Geometric_Y.write('%s\t%s\t%s\n'%('nodeA', 'nodeB', 'Geometric'))
            fw_Cosine_Y.write('%s\t%s\t%s\n'%('nodeA', 'nodeB', 'Cosine'))
            fw_PCC_Y.write('%s\t%s\t%s\t%s\n'%('nodeA', 'nodeB', 'PCC', 'positive/negative'))
            fw_Jaccard_Z.write('%s\t%s\t%s\n'%('nodeA', 'nodeB', 'Jaccard'))
            fw_Simpson_Z.write('%s\t%s\t%s\n'%('nodeA', 'nodeB', 'Simpson'))
            fw_Geometric_Z.write('%s\t%s\t%s\n'%('nodeA', 'nodeB', 'Geometric'))
            fw_Cosine_Z.write('%s\t%s\t%s\n'%('nodeA', 'nodeB', 'Cosine'))
            fw_PCC_Z.write('%s\t%s\t%s\t%s\n'%('nodeA', 'nodeB', 'PCC', 'positive/negative'))
            for (A_name, B_name) in combinations(list_typeX_nodes, 2):
                AI = calculate_association_index(Interaction_XY.get(A_name), Interaction_XY.get(B_name), Total_typeY_nodes)
                fw_Jaccard_Y.write('%s\t%s\t%f\n'%(A_name, B_name, AI.get('Jaccard')))
                fw_Simpson_Y.write('%s\t%s\t%f\n'%(A_name, B_name, AI.get('Simpson')))
                fw_Geometric_Y.write('%s\t%s\t%f\n'%(A_name, B_name, AI.get('Geometric')))
                fw_Cosine_Y.write('%s\t%s\t%f\n'%(A_name, B_name, AI.get('Cosine')))
                fw_PCC_Y.write('%s\t%s\t%f\t%s\n'%(A_name, B_name, AI.get('PCC'), AI.get('positive/negative')))
                AI = calculate_association_index(Interaction_XZ.get(A_name), Interaction_XZ.get(B_name), Total_typeZ_nodes)
                fw_Jaccard_Z.write('%s\t%s\t%f\n'%(A_name, B_name, AI.get('Jaccard')))
                fw_Simpson_Z.write('%s\t%s\t%f\n'%(A_name, B_name, AI.get('Simpson')))
                fw_Geometric_Z.write('%s\t%s\t%f\n'%(A_name, B_name, AI.get('Geometric')))
                fw_Cosine_Z.write('%s\t%s\t%f\n'%(A_name, B_name, AI.get('Cosine')))
                fw_PCC_Z.write('%s\t%s\t%f\t%s\n'%(A_name, B_name, AI.get('PCC'), AI.get('positive/negative')))
    print('association indice finish')
def association_index_histmatrix(PATH, Case_ID):
    BINS = 10
    Type_nodeX = 'lncRNA'
    Type_nodeY_1 = 'mRNA'
    Type_nodeY_2 = 'miRNA'
    for method in ['Jaccard', 'Simpson', 'Geometric', 'Cosine', 'PCC']:
        hist1 = {}
        hist2 = {}
        fname = '_'.join([Case_ID, '2'.join([Type_nodeX, Type_nodeY_1]), 'association_index', method])+'.txt'
        Data = pd.read_table('/'.join([PATH, fname]), sep='\t', header=0)
        Data = Data.loc[:, method].abs()
        Data = Data.iloc[list(np.nonzero(-np.isnan(Data))[0])]
        N = Data.shape[0]
        Bin_size = [int(np.fix(N/BINS)) for idx in range(BINS)]
        Bin_size[BINS-1] = int(np.fix(N/BINS)) + int(N-np.sum(Bin_size))
        argsort_data = np.argsort(Data)
        del Data
        idx = 0
        for bin_idx in range(BINS):
            for flag in range(Bin_size[bin_idx]):
                hist1.update({argsort_data.iloc[idx]:bin_idx})
                idx += 1
        del Bin_size, argsort_data
        fname = '_'.join([Case_ID, '2'.join([Type_nodeX, Type_nodeY_2]), 'association_index', method])+'.txt'
        Data = pd.read_table('/'.join([PATH, fname]), sep='\t', header=0)
        Data = Data.loc[:, method].abs()
        Data = Data.iloc[list(np.nonzero(-np.isnan(Data))[0])]
        N = Data.shape[0]
        Bin_size = [int(np.fix(N/BINS)) for idx in range(BINS)]
        Bin_size[BINS-1] = int(np.fix(N/BINS)) + int(N-np.sum(Bin_size))
        argsort_data = np.argsort(Data)
        del Data
        idx = 0
        for bin_idx in range(BINS):
            for flag in range(Bin_size[bin_idx]):
                hist2.update({argsort_data.iloc[idx]:bin_idx})
                idx += 1
        del Bin_size, argsort_data
        M = pd.DataFrame([list(map(int, np.zeros(BINS))) for idx in range(BINS)])
        for idx in hist1.keys():
            x_bin_idx = hist1[idx]
            y_bin_idx = hist2[idx]
            M.iloc[y_bin_idx, x_bin_idx] += 1
        fname = '_'.join([Case_ID, Type_nodeX, 'association_index_histmatrix', method])+'.txt'
        M.to_csv(str('/'.join([PATH, fname])), sep='\t', header=False, index=False)
        del M, hist1, hist2
    print('association_histmatrix finish')

def graph_association_index_histmatrix(PATH, Case_ID):
    sns.set(style='white', rc={'axes.edgecolor': '.2', 'axes.linewidth':1.5})
    fontpath = '/usr/share/fonts/truetype/msttcorefonts/Verdana.ttf'
    Fontprop = {'legend':fm.FontProperties(fname=fontpath, size=12), 
                'ticklabels':fm.FontProperties(fname=fontpath, size=10), 
                'label':fm.FontProperties(fname=fontpath, size=12), 
                'title':fm.FontProperties(fname=fontpath, size=14), 
                'cbar_ticklabels':fm.FontProperties(fname=fontpath, size=8)}
    DPI = 600
    CMAP = 'Spectral_r'
    xTickLabels = ['','','','','5','','','','','10']
    yTickLabels = ['10','','','','','5','','','','']
    
    Method = ['Jaccard', 'Simpson', 'Geometric', 'Cosine', 'PCC']
    for method in Method:
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111)
        cbar_ax = fig.add_axes([0.87, 0.25, 0.025, 0.5])
        fname = '_'.join([Case_ID, 'lncRNA_association_index_histmatrix', method])+'.txt'
        M = pd.read_table('/'.join([PATH, fname]), sep='\t', header=None, index_col=None)
        MAX_value = M.max().max()
        idx = sorted(range(M.shape[0]), reverse=True)
        M = M.iloc[idx,:]
        hm = sns.heatmap(M, vmin=0, vmax=MAX_value, cmap=CMAP, linewidths=.2, 
                         linecolor='white', cbar=True, square=True, ax=ax, 
                         xticklabels=xTickLabels, yticklabels=yTickLabels, 
                         cbar_ax=cbar_ax)
        plt.setp(cbar_ax.get_yticklabels(), fontproperties=Fontprop['cbar_ticklabels'])
        for item in hm.get_yticklabels():
            item.set_rotation(0)
        plt.setp(ax.get_xticklabels(), horizontalalignment='center', 
                 verticalalignment='top', fontproperties=Fontprop['ticklabels'])
        plt.setp(ax.get_yticklabels(), horizontalalignment='right', 
                 verticalalignment='center', fontproperties=Fontprop['ticklabels'])
        ax.set_ylabel('Similarity according to\nmiRNA', fontproperties=Fontprop['label'])
        ax.set_xlabel('Similarity according to mRNA', fontproperties=Fontprop['label'])
        fname = '_'.join([Case_ID, method, 'lncRNA_association_index_histmatrix_heatmap.png'])
        fig.savefig('/'.join([PATH, fname]), format='png', dpi=DPI)
        fig.show(False)
        plt.close()
    print('graph_finish')

def main(cmdArgument):
    print('test')
    CASE_ID = cmdArgument[1]
    PATH = cmdArgument[2]
    mRNA_file_name = CASE_ID + '_mRNA_file.txt'
    lncRNA_file_name = CASE_ID + '_lncRNA_file.txt'
    miRNA_file_name = CASE_ID + '_miRNA_file.txt'
    print('start function',CASE_ID,PATH)
    correlation(PATH, CASE_ID, mRNA_file_name, lncRNA_file_name, miRNA_file_name)
    bipartite_network(PATH, CASE_ID)
    association_indices(PATH, CASE_ID)
    association_index_histmatrix(PATH, CASE_ID)
    graph_association_index_histmatrix(PATH, CASE_ID)
    print('finish')
if __name__ == '__main__':
    main(sys.argv)
    
