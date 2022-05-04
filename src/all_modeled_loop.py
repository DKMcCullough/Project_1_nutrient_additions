'''

name:   All_modeled_loop.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/src

author: DKM


goal: loop modeling portion of this code


working on: Get params into loops so you can have differfent kdams for differentt lines


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
df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'

####################################

# Slicing Data

####################################

ros_uni = df_all['treatment'].unique()
treatments = ros_uni*1e-3   #1e-3 to take treatment nM (M is in L, not mL...so have to )

#np.array([0,40,400,4000,40000,400000])*1e-3   #1e-3 to take treatment nM (M is in L, not mL...so have to )


dc = dict()

####slicing df into treatments and saving into dictionary (dc)######
for count in range(treatments.shape[0]): 
    #print(count)
    t = treatments[count]
    df_t = df_all[df_all["treatment"].isin([t])]   #is this grabbing all the corrrect data???I can't tell.
    print(df_t)
    name = ('df_' + str(t))
    dc.update({name : df_t})  #update dictionary of dfs with each loop itteration. 


avgs = dict()
yerrs = dict()

#calculating avg point and yerr for each treatment

for d in dc :
    df_d = dc[d]   
    rep_df_d = df_d[['rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6']]
    avg_d = rep_df_d.mean(axis=1)  
    avgs.update({'avg_'+d : avg_d })
    yerr_d = rep_df_d.std(axis=1)   
    yerrs.update({'yerr_'+d : yerr_d })


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

step = 0.002
ndays = 35
times = np.linspace(0,ndays,int(ndays/step))

Qn = (9.6e-15*(1/(14.0))*1e+9)   #Nitrogen Quota for Pro from Bertillison? 
P = 1e4
k1s = [0.9, 0.9, 0.9, 0.9, 0.9,0.9]  #[0.02,]     
k2s = [0.5, 0.5, 0.5, 0.5, 0.5,0.5]   #[0.32]
kdams = [0.002,0.002, 0.002, 0.002, 0.002, 0.002]    #[0.007]
kddams = [0.115, 0.115, 0.115, 0.115, 0.115,0.115]    #[ 0.03,]

params = list(zip(k1s,k2s,kdams,kddams))
S_base = 3.0 
colors = ('orange', 'r', 'green', 'c', 'purple', 'k') #make into a set in loop? for c in colors, color = count(c)?????
markers = ('s','v','P','o','*','d')
##################################

#   Integration  

##################################

fig1,ax1  = plt.subplots(figsize=(10,7))
fig2,ax2  = plt.subplots(figsize=(12,7))
#fig3,ax3 = plt.subplots(figsize=(12,7))

    #nM N per ml for units      #    0.164 micromolar rediual N from Calfee_et_al 2022

mc = np.array([])

for count in range(treatments.shape[0]): 
    print(count)
    t = treatments[count]
    k1 = params[count][0]
    k2 = params[count][1]
    kdam = params[count][2]
    kddam = params[count][3]
    ks = (k2/k1)   #set max value for this that is know from lit? (between 0.01 and 0.015 for N metabolism in )
    SsEuler = np.array([])
    PsEuler = np.array([])
    P = 1e4
    S = (S_base + t )
    print(kdam,kddam)
    for t in times:
        PsEuler = np.append(PsEuler,P)
        SsEuler = np.append(SsEuler,S)
        if (S>1e-4):
            delta = kdam
        else:
            delta = kddam 
        dPdt = k2 * P * S /( (k2/k1) + S) - delta*P
        dSdt = -P*(( k2*S)/((k2/k1)+S))*Qn
        #if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
         #       S = 4e-47
        #else:
        S = S + dSdt*step
        P = P + dPdt*step
    ax1.plot(times,(PsEuler), linestyle = 'dashed', color = colors[count]) 
    ax2.plot(times,(SsEuler), linestyle = 'dashed', color = colors[count])



####################################

#graphing data

####################################

for count in range(treatments.shape[0]): 
    #print(count)
    t = treatments[count]
    df = dc['df_'+str(t)]
    times = df['times']
    data = avgs['avg_df_'+ str(t)]
    yerr_graph = yerrs['yerr_df_'+ str(t)]
    #m = max(data)
    print(count)
    ax1.plot(times, data, linestyle = 'None', marker= markers[count],  markersize= 12, label = (str(t) +' nM NH4'), color = colors[count])  
    ax1.plot(times,data,linestyle='-', linewidth=0.25, color='black', marker = 'None')
    ax1.errorbar(times, data, yerr = yerr_graph, fmt='none', color = colors[count])   

#Axes.set_ylim(self, ymin=1e-20, ymax=1e20)
ax1.set_ybound( lower=10**-20, upper=10*20)
ax1.set(xlabel= 'Time (days)', ylabel='Biomass (cells  ml$^{-1}$)', yscale = "log")
ax1.set_title('Prochlorococcus Biomass over time', fontsize=25)
ax1.legend(loc='upper left',prop={'size': 12}, fontsize=22)
ax1.set_xlabel(xlabel= 'Time (days)', fontsize=20)
ax1.set_ylabel(ylabel= 'Biomass (cells  ml$^{-1}$)', fontsize=20)
plt.semilogy()

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)


ax2.set(xlabel='Time (days)', ylabel=' Nutrient Concentration (nmol ml$^{-1}$)', yscale = "log")
ax2.set_title('NH4 concentrations over time',fontsize=20)
#ax2.legend(loc='lower left',prop={'size': 10}, fontsize=12)
#ax2.yaxis.set_label_position("right")
#ax2.yaxis.tick_right()


plt.show()

#finding number of cells possible with each N treatment
'''
nc = np.array([])
for t in treatments: 
    new_cells = t*Qn
    nc = np.append(nc,new_cells)
    print(nc)

total  = (nc)+1e4




#Qn = 6.857142857142856e-07

#Qn*4e5 ->  0.27428571428571424

a = (Qn*4e5)

#(1e4)*a   -> 2742.8571428571427

b = ((1e4)*a)

b*1e9 ->  2742857142857.1426

cells_possible  = b*1e9
'''
#  Find total cells made in each trial from data 
#from this we could find the actual amount of N used. 
'''
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


print('\n' + 'done')




#to do 
#nM in L not ml so need to change addendments for maths 
#HOOH vs kdam 