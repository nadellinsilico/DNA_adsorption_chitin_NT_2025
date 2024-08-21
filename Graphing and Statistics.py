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
#%%
"""
Figure 1C
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure1C')

plot_dict = {'No_Label':df['Label-_DNA+'].dropna(), 
             'Label':df['Label+_DNA-'].dropna(),
             'DNA':df['Label+_DNA+'].dropna()}

plt.figure(figsize=(2.9,3))
plt.boxplot(plot_dict.values(), widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2,3], labels=plot_dict.keys(), fontsize=12)

plt.scatter(np.ones(len(plot_dict['No_Label'])), plot_dict['No_Label'], alpha=0.5, color='black', s=75)
plt.scatter(np.ones(len(plot_dict['Label']))*2, plot_dict['Label'], alpha=0.5, color='red', s=75)
plt.scatter(np.ones(len(plot_dict['DNA']))*3, plot_dict['DNA'], alpha=0.5, color='red', s=75)

plt.ylim(0,2)
plt.ylabel(r'Label Intensity (A.U.)/area', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.savefig('Figure_1C.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['DNA'], plot_dict['Label']))
print(mannwhitneyu(plot_dict['No_Label'], plot_dict['Label']))
#%%
"""
Figure 1D
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure1D')

plot_dict = {'DNA 1 d':df['DNA+_24h'].dropna(), 
             'DNA 2 d':df['DNA+_48h'].dropna(),
             'Blank 1 d':df['DNA-_24h'].dropna(),
             'Blank 2 d':df['DNA-_48h'].dropna()}

plt.figure(figsize=(3,3))
plt.boxplot(plot_dict.values(), widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2, 3, 4], labels=plot_dict.keys(), fontsize=20, rotation=45)

plt.scatter(np.ones(len(plot_dict['DNA 1 d'])), plot_dict['DNA 1 d'], alpha=0.5, color='red', s=75)
plt.scatter(np.ones(len(plot_dict['DNA 2 d']))*2, plot_dict['DNA 2 d'], alpha=0.5, color='red', s=75)
plt.scatter(np.ones(len(plot_dict['Blank 1 d']))*3, plot_dict['Blank 1 d'], alpha=0.5, color='black', s=75, marker='D')
plt.scatter(np.ones(len(plot_dict['Blank 2 d']))*4, plot_dict['Blank 2 d'], alpha=0.5, color='black', s=75, marker='D')

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

print(mannwhitneyu(plot_dict['DNA 1 d'], plot_dict['Blank 1 d']))
print(mannwhitneyu(plot_dict['DNA 2 d'], plot_dict['Blank 2 d']))


#%%
"""
Figure 1F
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure1F')

plot_dict = {'DNA':df['DNA+_48h'].dropna(), 
             'Blank':df['DNA-_48h'].dropna()}

plt.figure(figsize=(1.6,3))
plt.boxplot(plot_dict.values(),widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2], labels=plot_dict.keys(), fontsize=12, rotation=45)

plt.scatter(np.ones(len(plot_dict['DNA'])), plot_dict['DNA'], alpha=0.5, color='red', s=75)
plt.scatter(np.ones(len(plot_dict['Blank']))*2, plot_dict['Blank'], alpha=0.5, color='black', s=75, marker='D')

#average biovolume of reporter strain (for both pilU and WT) at 4 d is ~10^4, -> limit of detection is set to 10^-4
#to put the data on a log scale 0 values for the GFP channel, the transformaiton reporter channel, were set to 1x10^-4
plt.hlines(y=10**-4, xmin=0.5, xmax=2.5, linewidth=1, linestyles='--', color='r')

plt.yticks(np.arange(0,0.1,0.02),fontsize=12)
plt.ylim(10**-5, 10**0)
plt.yscale('log')
plt.ylabel(r'Transformation Frequency', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.savefig('Figure1F.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

#untransform the data (10**-4 becomes 0)
for key in list(plot_dict.keys()):
    holder = []
    for i in plot_dict[key]:
        if i == 10**-4:
            holder.append(0)
        else:
            holder.append(i)
    plot_dict[key] = holder

print(mannwhitneyu(plot_dict['DNA'], plot_dict['Blank']))

#%%
"""
Figure 2A
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure2A')

plot_dict = {'WT':df['WT'].dropna(), 
             'pilU':df['pilU'].dropna()}

plt.figure(figsize=(2,3))
plt.boxplot(plot_dict.values(), widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2], labels=plot_dict.keys(), fontsize=12)

plt.scatter(np.ones(len(plot_dict['WT'])), plot_dict['WT'], alpha=0.5, color='black', s=75)
plt.scatter(np.ones(len(plot_dict['pilU']))*2, plot_dict['pilU'], alpha=0.5, color='blue', s=75)

#average biovolume of reporter strain (for both pilU and WT) at 4 d is ~10^4, -> limit of detection is set to 10^-4
#to put the data on a log scale 0 values for the GFP channel, the transformaiton reporter channel, were set to 1x10^-4
plt.hlines(y=10**-4, xmin=0.5, xmax=2.5, linewidth=1, linestyles='--', color='r')

plt.xticks(rotation=45)
plt.ylim(10**-5, 10**0)
plt.yscale('log')
plt.ylabel(r'Transformation Efficiency', fontsize=12)
plt.xlabel(r'Strain', fontsize = 12)
plt.savefig('Figure2A.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

#untransform the data (10**-4 becomes 0)
for key in list(plot_dict.keys()):
    holder = []
    for i in plot_dict[key]:
        if i == 10**-4:
            holder.append(0)
        else:
            holder.append(i)
    plot_dict[key] = holder

test_res = mannwhitneyu(plot_dict['WT'], plot_dict['pilU'])
print(test_res)

#%%
"""
Figure 2B
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Figure2B')

plot_dict = {'WT':df['WT'].dropna(), 
             'pilU':df['pilU'].dropna()}

plt.figure(figsize=(2,3))
plt.boxplot(plot_dict.values(), widths=0.5, showfliers=False)
plt.xticks(ticks=[1,2], labels=plot_dict.keys(), fontsize=12)

plt.scatter(np.ones(len(plot_dict['WT'])), plot_dict['WT'], alpha=0.5, color='black', s=75)
plt.scatter(np.ones(len(plot_dict['pilU']))*2, plot_dict['pilU'], alpha=0.5, color='blue', s=75)
plt.axhline(y=1/2500, color='r', linestyle='--', linewidth=1)

plt.ylim(10**-5, 10**0)
plt.yscale('log')
plt.ylabel(r'Transformation Efficiency', fontsize=12)
plt.xlabel(r'Strain', fontsize = 12)

plt.savefig('Figure2B.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['WT'], plot_dict['pilU']))