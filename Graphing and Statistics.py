# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:07:40 2024

@author: Jacob D. Holt

Graphing and statistics for acuumulation and recobmination of DNA on chitin
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

"""
https://stackoverflow.com/questions/36153410/how-to-create-a-swarm-plot-with-matplotlib
"""
def simple_beeswarm2(y, nbins=10, width=.25):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """
    y = np.asarray(y)
    if nbins is None:
        nbins = len(y) // 6

    # Get upper bounds of bins
    x = np.zeros(len(y))
    ylo = np.min(y)
    yhi = np.max(y)
    dy = (yhi - ylo) / nbins
    ybins = np.linspace(ylo + dy, yhi - dy, nbins - 1)

    # Divide indices into bins
    i = np.arange(len(y))
    ibs = [0] * nbins
    ybs = [0] * nbins
    nmax = 0
    for j, ybin in enumerate(ybins):
        f = y <= ybin
        ibs[j], ybs[j] = i[f], y[f]
        nmax = max(nmax, len(ibs[j]))
        f = ~f
        i, y = i[f], y[f]
    ibs[-1], ybs[-1] = i, y
    nmax = max(nmax, len(ibs[-1]))

    # Assign x indices
    try:
        dx = width / (nmax // 2)
    except:
        dx = 0.25
    for i, y in zip(ibs, ybs):
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(y)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x

"""
for log plot
transform the data (0 becomes 10**-5)
maximum biovolume of reporter strain (for both pilU and WT) at 4 d is ~10^5, -> limit of detection is set to 10^-5
to put the data on a log scale 0 values for the GFP channel, the transformaiton reporter channel, are set to 1x10^-5
"""
def transform(dict_in):
    for key in list(dict_in.keys()):
        holder = []
        for i in dict_in[key]:
            if i == 0:
                holder.append(10**-5)
            else:
                holder.append(i)
        dict_in[key] = holder
    return dict_in
#%%
"""
Figure 1C
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure1C')

plot_dict = {'No_Label':df['Label-_DNA+'].dropna(), 
             'Label':df['Label+_DNA-'].dropna(),
             'DNA':df['Label+_DNA+'].dropna()}

print(mannwhitneyu(plot_dict['DNA'], plot_dict['Label']))
print(mannwhitneyu(plot_dict['DNA'], plot_dict['No_Label']))

plt.figure(figsize=(2.9,3))
plt.boxplot(plot_dict.values(), widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2,3], labels=plot_dict.keys(), fontsize=12)

plt.scatter(simple_beeswarm2(plot_dict['No_Label'])+1, plot_dict['No_Label'], alpha=0.5, color='black', s=75)
plt.scatter(simple_beeswarm2(plot_dict['Label'])+2, plot_dict['Label'], alpha=0.5, color='red', s=75)
plt.scatter(simple_beeswarm2(plot_dict['DNA'])+3, plot_dict['DNA'], alpha=0.5, color='red', s=75)

plt.ylim(0,2)
plt.ylabel(r'Label Intensity (A.U.)/area', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.savefig('Figure1C.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

#%%
"""
Figure 1D
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure1D')

plot_dict = {'DNA 1 d':df['DNA+_24h'].dropna(), 
             'DNA 2 d':df['DNA+_48h'].dropna(),
             'Blank 1 d':df['DNA-_24h'].dropna(),
             'Blank 2 d':df['DNA-_48h'].dropna()}

print(mannwhitneyu(plot_dict['DNA 1 d'], plot_dict['Blank 1 d']))
print(mannwhitneyu(plot_dict['DNA 2 d'], plot_dict['Blank 2 d']))

plt.figure(figsize=(3,3))
plt.boxplot(plot_dict.values(), widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2, 3, 4], labels=plot_dict.keys(), fontsize=20, rotation=45)

plt.scatter((simple_beeswarm2(plot_dict['DNA 1 d']))+1, plot_dict['DNA 1 d'], alpha=0.5, color='red', s=75)
plt.scatter((simple_beeswarm2(plot_dict['DNA 2 d']))+2, plot_dict['DNA 2 d'], alpha=0.5, color='red', s=75)
plt.scatter((simple_beeswarm2(plot_dict['Blank 1 d']))+3, plot_dict['Blank 1 d'], alpha=0.5, color='black', s=75, marker='D')
plt.scatter((simple_beeswarm2(plot_dict['Blank 2 d']))+4, plot_dict['Blank 2 d'], alpha=0.5, color='black', s=75, marker='D')

#the limit of detection for the qubit based off of our dilution factor
plt.hlines(y=0.05, xmin=0.5, xmax=4.5, linewidth=1.25, linestyles='--', color='r')

plt.yticks(np.arange(0,0.22,.04),fontsize=12)
plt.ylim(0,.2)
plt.ylabel(r'DNA ng/uL', fontsize=12)
plt.xlabel(r'Treatment', fontsize=12)
plt.savefig('Figure1D.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

#%%
"""
Figure 1F
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure1F')

plot_dict = {'DNA':df['DNA+_48h'].dropna(), 
             'Blank':df['DNA-_48h'].dropna()}

print(mannwhitneyu(plot_dict['DNA'], plot_dict['Blank']))

plot_dict = transform(plot_dict)
plt.figure(figsize=(1.6,3))
plt.boxplot(plot_dict.values(),widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2], labels=plot_dict.keys(), fontsize=12, rotation=45)

plt.scatter(simple_beeswarm2(plot_dict['DNA'])+1, plot_dict['DNA'], alpha=0.5, color='red', s=75)
plt.scatter(simple_beeswarm2(plot_dict['Blank'])+2, plot_dict['Blank'], alpha=0.5, color='black', s=75, marker='D')
plt.hlines(y=10**-5, xmin=0.5, xmax=2.5, linewidth=1, linestyles='--', color='r')

plt.yticks(np.arange(0,0.1,0.02),fontsize=12)
plt.ylim(10**-6, 10**0)
plt.yscale('log')
plt.ylabel(r'Transformation Frequency', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.savefig('Figure1F.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2A
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure2A')

plot_dict = {'WT':df['WT'].dropna(), 
             'pilU':df['pilU'].dropna()}

print(mannwhitneyu(plot_dict['WT'], plot_dict['pilU']))

plt.figure(figsize=(2,3))
plt.boxplot(plot_dict.values(), widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2], labels=plot_dict.keys(), fontsize=12)

plt.scatter((simple_beeswarm2(plot_dict['WT'], nbins=1)+1), plot_dict['WT'], alpha=0.5, color='black', s=75)
plt.scatter((simple_beeswarm2(plot_dict['pilU'])+2), plot_dict['pilU'], alpha=0.5, color='blue', s=75)
plt.axhline(y=1/2500, color='r', linestyle='--', linewidth=1)

plt.ylim(10**-6, 10**0)
plt.yscale('log')
plt.ylabel(r'Transformation Efficiency', fontsize=12)
plt.xlabel(r'Strain', fontsize = 12)

plt.savefig('Figure2A.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2B
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure2B')

plot_dict = {'WT':df['WT'].dropna(), 
             'pilU':df['pilU'].dropna()}

print(mannwhitneyu(plot_dict['WT'], plot_dict['pilU']))

plot_dict = transform(plot_dict)
plt.figure(figsize=(2,3))
plt.boxplot(plot_dict.values(), widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2], labels=plot_dict.keys(), fontsize=12)

plt.scatter((simple_beeswarm2(plot_dict['WT']))+1, plot_dict['WT'], alpha=0.5, c='black', s=75)
plt.scatter((simple_beeswarm2(plot_dict['pilU']))+2, plot_dict['pilU'], alpha=0.5, c='blue', s=75)

plt.hlines(y=10**-5, xmin=0.5, xmax=2.5, linewidth=1, linestyles='--', color='r')

plt.xticks(rotation=45)
plt.ylim(10**-6, 10**0)
plt.yscale('log')
plt.ylabel(r'Transformation Efficiency', fontsize=12)
plt.xlabel(r'Strain', fontsize = 12)
plt.savefig('Figure2B.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()