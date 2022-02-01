'''

name:   All_graphing.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/src

author: DKM


goal: import and visualise all treatments of NH4 addition experients using MIT9215

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


#Ben_data_import = pd.read_excel('../Data/Ben_data_import_fixed.xlsx', header=0)

df_all = pd.read_csv("/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/data/NH4_add.csv")

#print(df_all)

df_all['rep1'] = df_all['rep1'].fillna(value = 0.0) #filling Nans with 09.0 in 'rep1' column 
df_all['rep2'] = df_all['rep2'].fillna(value = 0.0 )#filling Nans with 09.0 in 'rep2' column 



df_all = df_all.dropna(axis = 1)     #taking NaN columns off the end of df but had to fill rep 1 and 2 Nans first

print(df_all)

df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'




####################################

# Slicing Data

####################################

print(df_all['treatment'].value_counts())   #finding how many uniquie counts (different treatments) there are and the number of each

#df_smol = df_all[df_all["treatment"].isin([0,40])]
df_0 = df_all[df_all["treatment"].isin([0])]
df_40 = df_all[df_all["treatment"].isin([40])]
df_400 = df_all[df_all["treatment"].isin([400])]
df_4000 = df_all[df_all["treatment"].isin([4000])]
df_40000 =  df_all[df_all["treatment"].isin([40000])]
df_400000 = df_all[df_all["treatment"].isin([400000])]




rep_cols = ['rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6']     # columns of just replicate assay abundance values
avg_0 = df_0[rep_cols].mean(axis=1) #takes mena value across rep1-6 column first each row
avg_40 = df_40[rep_cols].mean(axis=1)
avg_400 = df_400[rep_cols].mean(axis=1)
avg_4000 = df_4000[rep_cols].mean(axis=1)
avg_40000 = df_40000[rep_cols].mean(axis=1)
avg_400000 = df_400000[rep_cols].mean(axis=1) 




######################################

#  Graphing Data                 

#####################################

'''
plt.figure()           #graphed all treatments in same rep column each time.
plt.scatter(df_all['times'],df_all['rep1'], label = 'rep1')
plt.scatter(df_all['times'],df_all['rep2'], label = 'rep2')
plt.scatter(df_all['times'],df_all['rep3'], label = 'rep3')
plt.scatter(df_all['times'],df_all['rep4'], label = 'rep4')
plt.scatter(df_all['times'],df_all['rep5'], label = 'rep5')
plt.scatter(df_all['times'],df_all['rep6'], label = 'rep6')

'''

fig , ax = plt.subplots(sharex=True, sharey=True) 
#plt.scatter(x = df_0['times'], y = [df_0['rep1'], label = '0 NH4 added')
plt.scatter(x = df_0['times'], y = [avg_0], label = '0 NH4 added')
plt.scatter(x = df_40['times'], y = [avg_40], label = '40 NH4 added')
plt.scatter(x = df_400['times'], y = [avg_400], label = '400 NH4 added')
plt.scatter(x = df_4000['times'], y = [avg_4000], label = '4000 NH4 added')
plt.scatter(x = df_40000['times'], y = [avg_40000], label = '40000 NH4 added')
plt.scatter(x = df_400000['times'], y = [avg_400000], label = '400000 NH4 added')

#plt.errorbar(x = df_0['times'], y = [avg_0], yerr = 'df_0[reps_cols].std(axis=1)')    #trying to get error bars to print on graph

plt.semilogy()

plt.legend()
plt.title('NH4 additions to MIT9215 Pro', fontsize = '22')
plt.xlabel('Time',fontsize = '18')
plt.legend(prop={"size":14})
plt.ylabel('Cell Abundance (flow)',fontsize = '18')
plt.xticks(fontsize = 14) 
plt.yticks(fontsize = 14)



plt.show()




print("Done") 

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.integrate import *
from scipy import *
from pylab import *
from matplotlib import *



##################################

# np data arrays

#################################

step = 1
ndays = 80
times = np.linspace(0,ndays,int(ndays/step))

# pd data frames 

#df = 



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

P = 1e3
S = 1e6 
k1_big  = 0.0000006

k2 = 1.2       # seems to control steepness of slope


big_SsEuler = np.array([])
big_PsEuler = np.array([])

#big alpha

for t in times:
        big_PsEuler = np.append(big_PsEuler,P)
        big_SsEuler = np.append(big_SsEuler,S)
        dPdt = k2 * P * S /( (k2/k1_big) + S)
        dSdt =-P*( k2*S)/((k2/k1_big)+S)
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


fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Alpha (k1) control on dynamics')
ax1.plot(times,big_PsEuler, label = "Phytoplankton Biomass over time")
ax1.set(xlabel='Time (day $^(-1)$', ylabel='number of cells(10^_)')
#ax1.xlabel('Time (day $^(-1)$')
#ax1.ylabel('number of cells')
ax2.plot(times, big_SsEuler, label = "Nutrient Concentration over time")
ax2.set(xlabel='Time (day $^(-1)$', ylabel='Nutrient concentration(10^_)')
#ax1.xlabel('Time (day $^(-1)$')
#ax1.ylabel('Nutrient concentration')

ax1.semilogy()
ax2.semilogy()
'''
plt.plot(times, PsEuler,label = "Phytoplankton Biomass")
plt.xlabel('Time (day $^{-1}$)') 
plt.ylabel('number of cells')

'''

plt.show()



"""

still working on getting the S to be hella high  
