'''

name:   modeled_low.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/src

author: DKM


goal: import and visualise and model NH4 addition experient using MIT9215

'''

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

df_all = pd.read_csv("/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/data/NH4_add.csv")

#print(df_all)

df_all['rep1'] = df_all['rep1'].fillna(value = 0.0) #filling Nans with 0.0 in 'rep1' column 
df_all['rep2'] = df_all['rep2'].fillna(value = 0.0 )#filling Nans with 0.0 in 'rep2' column 


df_all = df_all.dropna(axis = 1)     #taking NaN columns off the end of df but had to fill rep 1 and 2 Nans first

#print(df_all)

df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'


####################################

# Slicing Data

####################################

#print(df_all['treatment'].value_counts())   #finding how many uniquie counts (different treatments) there are and the number of each

df_40 = df_all[df_all["treatment"].isin([40])]



rep_cols = ['rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6']     # columns of just replicate assay abundance values
avg_40 = df_40[rep_cols].mean(axis=1) #takes mena value across rep1-6 column first each row 


##################################

# np data arrays

#################################

step = 1
ndays = 40
times = np.linspace(0,ndays,int(ndays/step))


################################

# Model


################################



#dPdt = (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration) (Population size) 
#dSdt = (Supply of nutriet) - (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration)

#dPdt = k2 * P * S /( (k2/k1) + S)       # k2 = Vmax  K1 = affinity for nutirent (alpha) 
#dSdt =  -P*( k2*S)/((k2/k1)+S)




##################################

#   Integration  

##################################

P = 1e4
S = 7e6 
k1= 0.00000078

k2 = 0.6       # seems to control steepness of slope

nrdelta = 0.02    #nutrient replete delta
nddelta = 0.1    #nutrient deplete delta

SsEuler40 = np.array([])
PsEuler40 = np.array([])

for t in times:
        PsEuler40 = np.append(PsEuler40,P)
        SsEuler40 = np.append(SsEuler40,S)
        if (S>1e-4):
            delta = nrdelta
        else:
            delta = nddelta
        dPdt = k2 * P * S /( (k2/k1) + S) - delta*P
        dSdt =-P*( k2*S)/((k2/k1)+S)
        if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                S = 0.00000000000000000004
        else:
                S = S + dSdt*step
        P = P + dPdt*step
              
        
        

#####################################

#  Graphing

#####################################


fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('N and P during NH addition trials')
ax1.plot(times,PsEuler40, label = "Prochlorococcus Biomass over time")
ax1.plot(df_40['times'],avg_40,linestyle = 'None',  marker='o', label = '40 NH4 added')
ax1.set(xlabel='Time (day $^(-1)$', ylabel='number of cells(10^_)')
ax2.plot(times, SsEuler40, label = "Nutrient Concentration over time")
ax2.set(xlabel='Time (day $^(-1)$', ylabel='Nutrient concentration(10^_)')

ax1.semilogy()
ax2.semilogy()
'''
plt.plot(times, PsEuler,label = "Phytoplankton Biomass")
plt.xlabel('Time (day $^{-1}$)') 
plt.ylabel('number of cells')

'''

plt.show()




print("Done")
