'''

name:   All_modeled.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/src

author: DKM


goal: import and all treatments of NH4 addition experients using MIT9215

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
df_all['rep1'] = df_all['rep1'].fillna(value = 0.0) #filling Nans with 09.0 in 'rep1' column 
df_all['rep2'] = df_all['rep2'].fillna(value = 0.0 )#filling Nans with 09.0 in 'rep2' column 
df_all = df_all.dropna(axis = 1)     #taking NaN columns off the end of df but had to fill rep 1 and 2 Nans first

#print(df_all)

df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'




####################################

# Slicing Data

####################################

#print(df_all['treatment'].value_counts())   #finding how many uniquie counts (different treatments) there are and the number of each

treatments = [0,40,400,4000,40000,400000]  
    
df_0 = df_all[df_all["treatment"].isin([0])]
df_40 = df_all[df_all["treatment"].isin([40])]
df_400 = df_all[df_all["treatment"].isin([400])]
df_4000 = df_all[df_all["treatment"].isin([4000])]
df_40000 =  df_all[df_all["treatment"].isin([40000])]
df_400000 = df_all[df_all["treatment"].isin([400000])]


rep_cols = ['rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6']     # columns of just replicate assay abundance values
avg_0 = df_0[rep_cols].mean(axis=1) #takes mena value across rep1-6 column for each row
avg_40 = df_40[rep_cols].mean(axis=1)
avg_400 = df_400[rep_cols].mean(axis=1)
avg_4000 = df_4000[rep_cols].mean(axis=1)
avg_40000 = df_40000[rep_cols].mean(axis=1)
avg_400000 = df_400000[rep_cols].mean(axis=1) 


yerr_0 = df_0[rep_cols].std(axis=1)
yerr_40 = df_40[rep_cols].std(axis=1)
yerr_400 = df_400[rep_cols].std(axis=1)
yerr_4000 = df_4000[rep_cols].std(axis=1)
yerr_40000 = df_40000[rep_cols].std(axis=1)
yerr_400000 = df_400000[rep_cols].std(axis=1)



##################################

# np data arrays

#################################

step = 1
ndays = 37
times = np.linspace(0,ndays,int(ndays/step))

P = 1e4
################################

# Model


################################

'''
#dPdt = (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration) (Population size) 
#dSdt = (Supply of nutriet) - (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration)

#dPdt = k2 * P * S /( (k2/k1) + S)    - delta*P     # k2 = Vmax  K1 = affinity for nutirent (alpha) 
#dSdt =  -P*( k2*S)/((k2/k1)+S)

#nutrient replete or deplete delta dependant on if S ~ 0.0
'''

##################################

#   Integration  

##################################


#0 NH4 added

P = 1e4
S = (0 + 1e6)    #    treatment to add in this case is 0
k1= 8e-7

k2 = 0.6      # seems to control steepness of slope

nrdelta = 0.002    #nutrient replete delta
nddelta = 0.01    #nutrient deplete delta


SsEuler0 = np.array([])
PsEuler0 = np.array([])


for t in times:
        PsEuler0 = np.append(PsEuler0,P)
        SsEuler0 = np.append(SsEuler0,S)
        if (S>1e-4):
            delta = nrdelta
        else:
            delta = nddelta
        dPdt = k2 * P * S /( (k2/k1) + S) - delta*P
        dSdt =-P*( k2*S)/((k2/k1)+S)
        if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                S = 4e-4
        else:
                S = S + dSdt*step
        P = P + dPdt*step
              
        
#40 

P = 1e4
S = (40 +2e6) 
k1= 8e-6

k2 = 0.6      # seems to control steepness of slope

nrdelta = 0.002    #nutrient replete delta
nddelta = 0.08   #nutrient deplete delta


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
                S = 4e-4
        else:
                S = S + dSdt*step
        P = P + dPdt*step
              
        
        

#400

P = 1e4
S = (4.0e2 + 5.8e6)   #4e2 is the added NH4 in the treatment
k1= 2.0e-7

k2 = 1.1      # seems to control steepness of slope

nrdelta = 0.03     #nutrient replete delta
nddelta = 0.19    #nutrient deplete delta


SsEuler400 = np.array([])
PsEuler400 = np.array([])


for t in times:
        PsEuler400 = np.append(PsEuler400,P)
        SsEuler400 = np.append(SsEuler400,S)
        if (S>1e-4):
            delta = nrdelta
        else:
            delta = nddelta
        dPdt = k2 * P * S /( (k2/k1) + S) - delta*P
        dSdt =-P*( k2*S)/((k2/k1)+S)
        if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                S = 4e-4
        else:
                S = S + dSdt*step
        P = P + dPdt*step
              
        
        
#4000

P = 1e4
S = (4.0e3 + 1.3e7)      #4e3 is the added treatement
k1= 2.0e-7

k2 = 0.9         # seems to control steepness of slope

nrdelta = 0.06    #nutrient replete delta
nddelta = 0.23    #nutrient deplete delta


SsEuler4000 = np.array([])
PsEuler4000 = np.array([])


for t in times:
        PsEuler4000 = np.append(PsEuler4000,P)
        SsEuler4000 = np.append(SsEuler4000,S)
        if (S>1e-4):
            delta = nrdelta
        else:
            delta = nddelta
        dPdt = k2 * P * S /( (k2/k1) + S) - delta*P
        dSdt =-P*( k2*S)/((k2/k1)+S)
        if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                S = 4e-4
        else:
                S = S + dSdt*step
        P = P + dPdt*step
              
        
        



#40000

P = 1e4
S = (4e4 + 2e3)#2e7 would be correct for assay set up
k1= 5e-7

k2 = 0.6      # seems to control steepness of slope

nrdelta = 0.002    #nutrient replete delta
nddelta = 0.2    #nutrient deplete delta


SsEuler40000 = np.array([])
PsEuler40000 = np.array([])


for t in times:
        PsEuler40000 = np.append(PsEuler40000,P)
        SsEuler40000 = np.append(SsEuler40000,S)
        if (S>1e-4):
            delta = nrdelta
        else:
            delta = nddelta
        dPdt = k2 * P * S /( (k2/k1) + S) - delta*P
        dSdt =-P*( k2*S)/((k2/k1)+S)
        if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                S = 4e-4
        else:
                S = S + dSdt*step
        P = P + dPdt*step
        
        
        
#400000

P = 1e4
S = (4e5+ 2e3) 
k1= 8e-8

k2 = 0.6      # seems to control steepness of slope

nrdelta = 0.002    #nutrient replete delta
nddelta = 0.01    #nutrient deplete delta


SsEuler400000 = np.array([])
PsEuler400000 = np.array([])


for t in times:
        PsEuler400000 = np.append(PsEuler400000,P)
        SsEuler400000 = np.append(SsEuler400000,S)
        if (S>1e-4):
            delta = nrdelta
        else:
            delta = nddelta
        dPdt = k2 * P * S /( (k2/k1) + S) - delta*P
        dSdt =-P*( k2*S)/((k2/k1)+S)
        if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                S = 4e-4
        else:
                S = S + dSdt*step
        P = P + dPdt*step
              
        
        

####################################

#graphing 

####################################


fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,7))
fig1.suptitle('NH4 trials modeled with data', fontweight='bold', fontsize=25)


#cell abundance subplot
#model
#ax1.plot(times,PsEuler0,color = 'm' , label  = 'zero NH4 added')
#ax1.plot(times,PsEuler40,color = 'r' , label = ' + 40 NH4 treatment')
ax1.plot(times,PsEuler400,color = 'green' , label = ' + 400 NH4 treatment')
ax1.plot(times,PsEuler4000, color = 'c', label = ' + 4000 NH4 added treatment')
#ax1.plot(times,PsEuler40000, color = 'b' , label = ' + 40000 NH4 treatment')
#ax1.plot(times,PsEuler400000 , color = 'k' , label = ' + 400000 NH4 treatment')

#data
#ax1.plot(df_0['times'], avg_0, linestyle = 'None',  marker='o', color = 'm' , label  = 'zero NH4 added')
#ax1.plot(df_40['times'], avg_40, linestyle = 'None',  marker='o', color = 'r' , label = ' + 40 NH4 treatment')
ax1.plot(df_400['times'], avg_400, linestyle = 'None',  marker='o', color = 'green' , label = ' + 400 NH4 treatment')
ax1.plot(df_4000['times'], avg_4000, linestyle = 'None',  marker='o', color = 'c' , label = ' + 4000 NH4 added treatment')
#ax1.plot(df_40000['times'], avg_40000, linestyle = 'None',  marker='o', color = 'b' , label = ' + 40000 NH4 treatment')
#ax1.plot(df_400000['times'], avg_400000, linestyle = 'None',  marker='o', color = 'k' , label = ' + 400000 NH4 treatment')

#errorbars
#ax1.errorbar(df_0['times'], avg_0, yerr=yerr_0,fmt='none', color = 'm')
#ax1.errorbar(df_40['times'], avg_40, yerr=yerr_40,fmt='none', color = 'r')
ax1.errorbar(df_400['times'], avg_400,yerr=yerr_400,fmt='none', color = 'green' )
ax1.errorbar(df_4000['times'], avg_4000, yerr=yerr_4000,fmt='none', color = 'c' )
#ax1.errorbar(df_40000['times'], avg_40000, yerr=yerr_40000,fmt='none', color = 'b')
#ax1.errorbar(df_400000['times'], avg_400000, yerr=yerr_400000,fmt='none', color = 'k' )



ax1.set(xlabel='Time (day $^-1$)', ylabel='number of cells (10^_)')
ax1.set_title('Prochlorococcus Biomass over time', fontsize=20)
ax1.legend(loc='lower right', fontsize=12)

   
#nutrient subplot
#ax2.plot(times,SsEuler0,color = 'm' , label  = 'zero NH4 added')
#ax2.plot(times,SsEuler40,color = 'r' , label = ' + 40 NH4 treatment')
ax2.plot(times,SsEuler400,color = 'green' , label = ' + 400 NH4 treatment')
ax2.plot(times,SsEuler4000, color = 'c' , label = ' + 4000 NH4 added treatment')
#ax2.plot(times,SsEuler40000, color = 'b' , label = ' + 40000 NH4 treatment')
#ax2.plot(times,SsEuler400000 , color = 'k' , label = ' + 400000 NH4 treatment')


ax1.semilogy()
ax2.semilogy()


ax2.set(xlabel='Time (day $^-1)$', ylabel='Nutrient concentration(10^_)')
ax2.set_title('NH4 concentrations over time',fontsize=20)
ax2.legend(loc='upper right', fontsize=12)



#plt.legend(prop={"size":14})

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
#plt.tick_params(axis='both', which = 'both', length = 2,  labelsize=16,)

#plt.add_subplot(3,1,1)


plt.show()




#To do 
'''

array of times, data/model, error array name, label for each treatment.

Then we can use treatment array to have all graphing things and can run a loop.  

Run treatment array through slicing if loop to cleam up code.

Individual plots for each treatment with printout of parameters as a sub subplot? 




'''
