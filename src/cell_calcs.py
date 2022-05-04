# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""



import pandas as pd
import numpy as np
from matplotlib import *
import matplotlib.pyplot as plt
from scipy.integrate import *
from scipy import *
from pylab import *



############################

#  Data Import from csv   

############################


df_all = pd.read_csv("../data/NH4_add.csv")
df_all['rep1'] = df_all['rep1'].fillna(value = 0.0) #filling Nans with 0.0 in 'rep1' column 
df_all['rep2'] = df_all['rep2'].fillna(value = 0.0 )#filling Nans with 0.0 in 'rep2' column 
df_all = df_all.dropna(axis = 1)     #taking NaN columns off the end of df but had to fill rep 1 and 2 Nans first
df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'

####################################

# Slicing Data

####################################

ros_uni = df_all['treatment'].unique()
treatments = ros_uni*1e-3   #1e-3 to take treatment nM (M is in L, not mL...so have to )

#print(treatments)
Qn = (9.6e-15*(1/(14.0))*1e+9)
P = 1e4

nc = np.array([])
for t in treatments: 
    new_cells = t/Qn
    nc = np.append(nc,new_cells)

print(list(zip(treatments,nc)))

'''   
total  = (nc)+1e4

ms = np.array([]) #made cells 
for t in treatments: 
    count = treatments.index(t)
    #print(count)
    df = dc['df_'+str(t)]
    times = df['times']
    data = avgs['avg_df_'+ str(t)]
    yerr_graph = yerrs['yerr_df_'+ str(t)]
    m = max(data)
    #print(m)
    mc = np.append(ms,m)
    #print(mc)


'''
