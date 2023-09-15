#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 23:13:57 2023
location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/src

author: DKM


goal: loop modeling portion of this code


working on: 
"""



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
from scipy.integrate import odeint
import ODElib
import random as rd




df_all = pd.read_csv("../data/low_N_growth_BCC_compiled.csv")

df_all = df_all.rename({'Time(days)':'time'}, axis=1)    #'renaming column to make it callable by 'times'


df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df = df_all
df['log_sigma'] = df['log_sigma'].clip(lower=0) 
df = df.rename(columns={"log_abundance": "abundance"})

inits = pd.read_csv("../data/inits/low_N_weeks_inits.csv")




#assays are weeks0-8
#graph Pro without N assition each time to get parameterizartiion of SN of media as well as K1 and k2 of Pro


