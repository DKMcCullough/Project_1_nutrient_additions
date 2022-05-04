'''

name:   control_modeled.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/src

author: DKM


get a real fix on zero treatment for NH4 data set of  MIT9215

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


df_all = pd.read_csv("../data/NH4_add.csv")
df_all['rep1'] = df_all['rep1'].fillna(value = 0.0) #filling Nans with 0.0 in 'rep1' column 
df_all['rep2'] = df_all['rep2'].fillna(value = 0.0 )#filling Nans with 0.0 in 'rep2' column 
df_all = df_all.dropna(axis = 1)     #taking NaN columns off the end of df but had to fill rep 1 and 2 Nans first

#make a for loop to fill all nans for all lines with 'rep' in them

df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'

####################################

# Slicing Data

####################################


treatments = [0,40,400,4000,40000,400000]  

dc = dict()

####slicing df into treatments and saving into dictionary (dc)######
for i in treatments:
    df_i = df_all[df_all["treatment"].isin([i])]
    name = ('df_' + str(i))
    dc.update({name : df_i})  #update dictionary of dfs with each loop itteration. 


avgs = dict()
yerrs = dict()

#calculating avg point and yerr for each treatment

for i in dc :
    df_i = dc[i]   
    rep_df_i = df_i[['rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6']]
    avg_i = rep_df_i.mean(axis=1)  
    avgs.update({'avg_'+i : avg_i })
    yerr_i = rep_df_i.std(axis=1)   
    yerrs.update({'yerr_'+i : yerr_i })


################################

# Model


################################
'''
#dPdt = (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration) (Population size) 
#dSdt = (Supply of nutriet) - (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration)
#dPdt = k2 * P * S /( (k2/k1) + S)    - delta*P     # k2 = Vmax  K1 = affinity for nutirent (alpha) 
#dSdt =  -P*( k2*S)/((k2/k1)+S) + (nutrinet concentration)*Cell quota
#Qn = (9.4e-15*(1/14.0)*1e+9)  #Nitrogen Quota for Pro 
#9.4*10^-15 g N per cell   #Nitrogen fg quota of Pro cells from Bertillison et al 2003
#14 g per Mol N     #N g to M converstion from periodic table
#10^9 ng per g      #conversion from grams to n grams
#nutrient replete or deplete delta dependant on if S ~ 0.0
'''

##################################

# np data arrays

#################################

step = 0.1
ndays = 35
times = np.linspace(0,ndays,int(ndays/step))
Qn = (9.4e-15*(1/(14.0))*1e+9)   #Nitrogen Quota for Pro from Bertillison? 

P = 1e4
S_base = 3.0    #nM N             #(treatment  = S

##################################

#   Integration  

##################################


#0 NH4 added


S = (S_base + 0 )    #nM N per ml for units      #    0.164 micromolar rediual N from Calfee_et_al 2022

k1= 0.9

k2 = 0.5       # seems to control steepness of slope

nrdelta = 0.002      #nutrient replete delta
nddelta = 0.115       #nutrient deplete delta


SsEuler0 = np.array([])
PsEuler0 = np.array([])


for t in times:
        PsEuler0 = np.append(PsEuler0,P)
        SsEuler0 = np.append(SsEuler0,S)
        if (S>1e-4):
            delta = nrdelta
        else:
            delta = nddelta
        dPdt = k2 * P * S /( (k2/k1) + S) - delta*P   # <--  should this also have Qn in it
        dSdt = -P*(( k2*S)/((k2/k1)+S))*Qn         # <--  should Qn be in fractrion on both sides next to S?
        if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                S = 4e-47
        else:
                S = S + dSdt*step
        P = P + dPdt*step


####################################

#graphing 

####################################


fig, (ax1) = plt.subplots(figsize=(10,7))

ax1.plot(times, PsEuler0, color = 'orange' , label  = 'model')   #model


df_0 = dc['df_0']              #slicing data for graphing ease
times = df_0['times']
data = avgs['avg_df_0']
yerr_graph = yerrs['yerr_df_0']


ax1.plot(times, data, linestyle = 'None', marker= 's',markersize = 12, label = ('0 nM NH4'), color = 'orange')  #color = colors(i))
ax1.errorbar(times, data, yerr = yerr_graph, fmt='none',color = 'orange')   
#ax1.plot(times,data,linestyle='-', linewidth=0.25, color='black', marker = 'None')


ax1.set(xlabel= 'Time (days)', ylabel='Biomass (cells  ml$^{-1}$)', yscale = "log")
ax1.set_title('Prochlorococcus Biomass over time', fontsize=20)
ax1.legend(loc='upper left',prop={'size': 12}, fontsize=16)

ax1.tick_params(axis = 'both', labelsize=20)
ax1.xaxis.label.set_size(18)
ax1.yaxis.label.set_size(18)

#pnames = ('alpha','vmax','nrdelta','nddelta')
#pvalues = (k1,k2,nrdelta,nddelta)

#ax2 = plt.table([pvalues], colLabels=(pnames) ,loc='lower center')
#ax2.set_fontsize(14)

plt.show()

print(S_base)







