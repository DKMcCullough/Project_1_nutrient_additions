'''
Rstar with data using Vmax
created by DKM
Location: /Users/dkm/Documents/Talmy_research/scripts/Producers_and_Resources/
'''


from scipy import *
from pylab import *
from numpy import *
from matplotlib import *
import sys
import xlrd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 


No_N_imp_practice_full = pd.read_excel("../../AMP_No_Nitrogen_import_practice.xlsx", header = 0)    #had to put 2 dots and 1 / to go back a single fo    lder destination  
No_N_imp_practice_Tech = pd.read_excel('../../AMP_No_Nitrogen_import_practice.xlsx','Vol 1 AMP No N. Tech. Reps')

print(No_N_imp_practice_full.head())  #printing just the first 5 rows of each column in  No_N_imp file 
time = No_N_imp_practice_Tech ['Time(days)']
No_N_A1 = No_N_imp_practice_Tech ['NoN-A.1']
print(time) #printing all of time column





rcParams["legend.loc"]


#create figures
#f,ax = subplots()

#model time array (sampling for n days with delta being the # of days between each sample?)
delt,ndays = 1.0/24.0,42.0
times = linspace(0,ndays,int(ndays/delt)) #model time series (start at 0, end with ndays, number of steps =  ndays/delta)

#model parameters
alpha =0.000001  #Producer's (P) affinity for Resource (R)
delta =0.12   #loss rate (1/day) of Producer (P)
Vmax = 0.58   #max velocity of enzyme being used for consumption of R

#steady states
Rstar = (delta*Vmax)/(Vmax-(alpha*delta))

#initial conditions
R = 6.9e+6
P = No_N_A1[0]

#biomass / water arrays
Rs,Ps, = r_[[]],r_[[]]    #creates 2 empty arrays named Rs and Ps
 

for t in times:
        Rs = append(Rs,R) #eight spaces over used for indention      # adds value to R found in each itteration to the end of the array Rs
        Ps = append(Ps,P)      #addes the new values of P to the end of the array Ps each time this is itterated
        dRst = -(R*Vmax/(R+(Vmax/alpha)))*P
        dPst = ((R*Vmax)/(R+(Vmax/alpha)))*P - delta*P
        R = max(R + dRst*delt,1e-4) #numerical integration (euler's method)
        P = max(P + dPst*delt,1e-4) 
        

#graphing
fig, axis1 = plt.subplots()   #making multiple plots that will all show at the same time when show() command is given
axis1.set_xlabel('Time (days)')  #labeling the x axis of the first subplot
axis1.set_ylabel('Cell Density (cells per mL)') #labeling the y axis of the first subplot
semilogy()
axis1.plot(times,Ps,c='g',label='Producers') #plotting Ps on first subplot
plt.plot(time,No_N_A1, marker='x', c='b', label='No N A1') #plotting points with a line going through them; print    ing on same axes as the scatter plot graph

#axis 1 subplot legend
legend = plt.legend(loc='best')   #, bbox_to_anchor=(40,1e+7))
legend.draw_frame(False)

axis2 = axis1.twinx()   #making it so the two subplots sit on top of one antoher by twining the x axis

axis2.set_ylabel('Nutrient Concentration') #labeling the y axis of the second subplot
semilogy()
#axis2.axhline(Rstar,c='r',ls='--',label='Rstar')
axis2.plot(times,Rs,c='r',label='Resource')

#axis 2 subplot legend
legend = plt.legend(loc='lower right')    #, bbox_to_anchor=(0.5,1.5))
legend.draw_frame(False)

plt.title('No N A1 Raw Data on Vmax Model')   #labels entire figure with title

show()

print('Done!!!')

