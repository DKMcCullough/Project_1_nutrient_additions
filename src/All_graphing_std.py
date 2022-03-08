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


#error calc

yerr_0 = df_0[rep_cols].std(axis=1)
yerr_40 = df_40[rep_cols].std(axis=1)
yerr_400 = df_400[rep_cols].std(axis=1)
yerr_4000 = df_4000[rep_cols].std(axis=1)
yerr_40000 = df_40000[rep_cols].std(axis=1)
yerr_400000 = df_400000[rep_cols].std(axis=1)



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
plt.scatter(x = df_0['times'], y = [avg_0], color = 'm', label = '0 NH4 ')
plt.scatter(x = df_40['times'], y = [avg_40], color = 'red',  label = '40 NH4 ')
plt.scatter(x = df_400['times'], y = [avg_400],  color = 'green' , label = '400 NH4 ')
plt.scatter(x = df_4000['times'], y = [avg_4000], color = 'c' ,  label = '4000 NH4 ')
plt.scatter(x = df_40000['times'], y = [avg_40000], color = 'b', label = '40000 NH4 ')
plt.scatter(x = df_400000['times'], y = [avg_400000],color = 'k',  label = '400000 NH4 ')

plt.errorbar(df_0['times'], avg_0, yerr=yerr_0,fmt='none', color = 'm')
plt.errorbar(df_40['times'], avg_40, yerr=yerr_40,fmt='none', color = 'red')
plt.errorbar(df_400['times'], avg_400,yerr=yerr_400,fmt='none', color = 'green' )
plt.errorbar(df_4000['times'], avg_4000, yerr=yerr_4000,fmt='none', color = 'c' )
plt.errorbar(df_40000['times'], avg_40000, yerr=yerr_40000,fmt='none', color = 'b')
plt.errorbar(df_400000['times'], avg_400000, yerr=yerr_400000,fmt='none', color = 'k' )


plt.semilogy()


plt.title('NH4 additions to MIT9215 Pro', fontsize = '22')
plt.xlabel('Time',fontsize = '18')
plt.legend(loc='upper left',prop={'size': 10}, fontsize=12)
plt.ylabel('Cell Abundance (flow)',fontsize = '18')
plt.xticks(fontsize = 14) 
plt.yticks(fontsize = 14)



plt.show()




print("Done") 


