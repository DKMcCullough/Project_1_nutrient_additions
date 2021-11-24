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



######################################

#  Graphing Data 

#####################################

plt.figure()
plt.scatter(df_all['times'],df_all['rep1'], label = 'rep1')
plt.scatter(df_all['times'],df_all['rep2'], label = 'rep2')
plt.scatter(df_all['times'],df_all['rep3'], label = 'rep3')
plt.scatter(df_all['times'],df_all['rep4'], label = 'rep4')
plt.scatter(df_all['times'],df_all['rep5'], label = 'rep5')
plt.scatter(df_all['times'],df_all['rep6'], label = 'rep6')




#df_40 = df_all.loc[(df_all['treatment'] == "40.0"),::]    #not working 


'''

time = df_main.iloc[:,1]
rep1 = df_main.iloc[:,2]
rep2 = df_main.iloc[:,3]

plt.semilogy()
'''

plt.legend()
plt.title('NH4 additions to MIT9215 Pro', fontsize = '22')
plt.xlabel('Time',fontsize = '18')
plt.legend(prop={"size":14})
plt.ylabel('Cell Abundance (flow)',fontsize = '18')
plt.xticks(fontsize = 14) 
plt.yticks(fontsize = 14)



plt.show()







