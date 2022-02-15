'''

name:   modeled_high.py 

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

df_400000 = df_all[df_all["treatment"].isin([400000])]



rep_cols = ['rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6']     # columns of just replicate assay abundance values
avg_400000 = df_400000[rep_cols].mean(axis=1) #takes mean value across rep1-6 column for each row 
yerr = df_400000[rep_cols].std(axis=1)

######################################

#  Graphing Data                 

#####################################
'''

fig , ax = plt.subplots(sharex=True, sharey=True) 
plt.scatter(x = df_400000['times'], y = [avg_400000], label = '400000 NH4 added')
yerr = df_400000[rep_cols].std(axis=1)

#plt.errorbar(x = df_0['times'], y = [avg_0], yerr = 'df_0[reps_cols].std(axis=1)')    #trying to get error bars to print on graph

plt.semilogy()

plt.legend()
plt.title('Replete N for MIT9215 Pro', fontsize = '22')
plt.xlabel('Time',fontsize = '18')
plt.legend(prop={"size":14})
plt.ylabel('Cell Abundance (cell/ml)',fontsize = '18')
plt.xticks(fontsize = 14) 
plt.yticks(fontsize = 14)



#plt.show()
'''

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
S = 2e8 
k1= 0.00000008

k2 = 0.6      # seems to control steepness of slope


SsEuler = np.array([])
PsEuler = np.array([])


for t in times:
        PsEuler = np.append(PsEuler,P)
        SsEuler = np.append(SsEuler,S)
        dPdt = k2 * P * S /( (k2/k1) + S)
        dSdt =-P*( k2*S)/((k2/k1)+S)
        if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                S = 0.00000000000000000004
        else:
                S = S + dSdt*step
        P = P + dPdt*step
                #S = S + dSdt*step
        #print( dPdt,t,P)
        #print(np.log(dSdt), t, S)
        #print(' \n')



#####################################

#  Graphing

#####################################


fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,7))
fig1.suptitle('400000 NH4 trial modeled data', fontweight='bold', fontsize=25)

#cell abundance subplot
ax1.plot(times,PsEuler, color = 'purple', label = "modeled Pro abundance, high N") #plotting model on same ax plot as data
ax1.plot(df_400000['times'],avg_400000,linestyle = 'None',  marker='o', color = 'green' , label = ' + 400000 NH4 treatment') #plotting data on ax 1 of fig1
ax1.errorbar(df_400000['times'], avg_400000, yerr=yerr,fmt='none', color = 'green' ) #error bars for avg400000 data. yerr
ax1.set(xlabel='Time (day $^-1$)', ylabel='number of cells(10^_)')
ax1.set_title('Prochlorococcus Biomass over time', fontsize=20)
ax1.legend(loc='lower right', fontsize=12)
   
#nutrient subplot
ax2.plot(times, SsEuler, color = 'purple', label = "high N model")
ax2.set(xlabel='Time (day $^-1)$', ylabel='Nutrient concentration(10^_)')
ax2.set_title('NH4 concentrations over time',fontsize=20)
ax2.legend(loc='lower left', fontsize=12)

ax1.semilogy()
ax2.semilogy()


plt.show()



print('k1')
print(k1)
print('k2')
print(k2) 



print("Done")
