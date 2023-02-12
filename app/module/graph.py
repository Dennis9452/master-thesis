# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 16:53:08 2017

@author: user
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from pathlib import Path
import os

def graph():
    path = Path('/usr/share/fonts/truetype/msttcorefonts/Verdana.ttf')
    fontprop = fm.FontProperties(fname=path, size=18)
    Fontsize = {'legend':10, 'ticklabels':10, 'label':12, 'title':14}
    sns.set(style='whitegrid', rc={'axes.edgecolor': '.2', 'axes.linewidth':1.5, 
    'grid.linestyle': ':', 'grid.linewidth':1})
    
    fig = plt.figure(figsize=(12, 9), dpi=100)
    ax = fig.add_subplot(111)
    # Create a random dataset across several variables
    rs = np.random.RandomState(0)
    n, p = 40, 8
    d = rs.normal(0, 2, (n, p))
    d += np.log(np.arange(1, p + 1)) * -5 + 10
    
    # Use cubehelix to get a custom sequential palette
    pal = sns.cubehelix_palette(p, rot=-.5, dark=.3)
    
    # Show each distribution with both violins and points
    sns.violinplot(data=d, palette=pal, inner="points")
    
    plt.setp(ax.get_xticklabels(), horizontalalignment='center', 
             verticalalignment='top', fontsize=Fontsize['ticklabels'], fontproperties=fontprop)
    plt.setp(ax.get_yticklabels(), horizontalalignment='right', 
             verticalalignment='center', fontsize=Fontsize['ticklabels'], fontproperties=fontprop)
    ax.set_xlabel('Spearman rank-order correlation coefficient', fontsize=Fontsize['label'], fontproperties=fontprop)
    ax.set_ylabel('Probability Density', fontsize=Fontsize['label'], fontproperties=fontprop)
    #os.path.join(UPLOAD_FOLDER, Case_ID+fname)os.path.join(UPLOAD_FOLDER, Case_ID+fname) 
    UPLOAD_FOLDER = '/var/www/helloworldapp/app/uploads'
    fname = 'temp.png'
    fig.savefig(os.path.join(UPLOAD_FOLDER, fname), format='png', dpi=500)
    fig.show(False)
    plt.close()
    return 0
