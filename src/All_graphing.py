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

df_all = pd.read_csv( "/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/data/NH4_add.csv")

print(df_all)



df_main = pd.read_csv("/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_6_transfer_datasets/data/Pro_transfer_no_N_fixed.csv")

#print(df_main.loc[time()])


#slicing into arrays since column names aren't callable rn :/ 

time = df_main.iloc[:,1]
rep1 = df_main.iloc[:,2]
rep2 = df_main.iloc[:,3]
rep3 = df_main.iloc[:,4]
rep4 = df_main.iloc[:,5]
rep5 = df_main.iloc[:,6]
rep6 = df_main.iloc[:,7]

#Graphing

plt.figure()

plt.plot(time,rep1,label = 'rep1')
plt.plot(time,rep2,label = 'rep2')
plt.plot(time,rep3,label = 'rep3')
plt.plot(time,rep4,label = 'rep4')
plt.plot(time,rep5,label = 'rep5')
plt.plot(time,rep6,label = 'rep6')


plt.semilogy()


plt.legend()
plt.title('MIT9215 Transfers', fontsize = '22')
plt.xlabel('Time',fontsize = '18')
plt.legend(prop={"size":14})
plt.ylabel('Cell Abundance (flow)',fontsize = '18')
plt.xticks(fontsize = 14) 
plt.yticks(fontsize = 14)



plt.show()







