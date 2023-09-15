'''
vmax_differences.py
created by DKM

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_4_k_derevations/src

messing with K1 nd K2 values while S is heading towards infinity to better understand how K1 and K2 control the curve. 

k2 = Velocity of reaction = possibly Vmax = dissociation constant in a singular enzymatic reaction
k2/k1 = Km 


'''


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
k1_smol = 0.00000015      # I am strugglign to get a good sense of what this does to the graph
k1_med  = 0.0000003
k1_big  = 0.0000006

k2_smol = 1.2       # seems to control steepness of slope
k2_med  = 1.5
k2_big  = 1.7


smol_SsEuler = np.array([])
smol_PsEuler = np.array([])
med_SsEuler = np.array([])
med_PsEuler = np.array([])
big_SsEuler = np.array([])
big_PsEuler = np.array([])



#Smol

for t in times:
	smol_PsEuler = np.append(smol_PsEuler,P)
	smol_SsEuler = np.append(smol_SsEuler,S)
	dPdt = k2_smol * P * S /( (k2_smol/k1_smol) + S)
	dSdt =-P*( k2_smol*S)/((k2_smol/k1_smol)+S)
	if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
		S = 0.00000000000000000004
	else:
		S = S + dSdt*step 
	P = P + dPdt*step
		#S = S + dSdt*step
	#print( dPdt,t,P)
	#print(np.log(dSdt), t, S)
	#print(' \n')



#med

for t in times:
        med_PsEuler = np.append(med_PsEuler,P)
        med_SsEuler = np.append(med_SsEuler,S)
        dPdt = k2_med * P * S /( (k2_med/k1_med) + S)
        dSdt =-P*( k2_med*S)/((k2_med/k1_med)+S)
        if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                S = 0.00000000000000000004
        else:
                S = S + dSdt*step
        P = P + dPdt*step
                #S = S + dSdt*step
        #print( dPdt,t,P)
        #print(np.log(dSdt), t, S)
        #print(' \n')


#big alpha

for t in times:
        big_PsEuler = np.append(big_PsEuler,P)
        big_SsEuler = np.append(big_SsEuler,S)
        dPdt = k2_big * P * S /( (k2_big/k1_big) + S)
        dSdt =-P*( k2_big*S)/((k2_big/k1_big)+S)
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
ax1.plot(times,smol_PsEuler, label = "Phytoplankton Biomass over time")
ax1.plot(times,med_PsEuler, label = "Phytoplankton Biomass over time")
ax1.plot(times,big_PsEuler, label = "Phytoplankton Biomass over time")
ax1.set(xlabel='Time (day $^(-1)$', ylabel='number of cells(10^_)')
#ax1.xlabel('Time (day $^(-1)$')
#ax1.ylabel('number of cells')
ax2.plot(times, smol_SsEuler, label = "Nutrient Concentration over time")
ax2.plot(times, med_SsEuler, label = "Nutrient Concentration over time")
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

"""
