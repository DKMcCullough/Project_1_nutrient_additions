'''
Rstar with data and alpha with preleminary data added
created by @DKM
Location: /Users/dkm/Documents/Talmy_research/scripts/Producers_and_Rosources/
'''


from scipy import *
from pylab import *
from matplotlib import *
import sys
import xlrd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 


No_N_imp_practice_full = pd.read_excel("../../AMP_No_Nitrogen_import_practice.xlsx", header = 0)    #had to put 2 dots and 1 / to go back a single folder destination  
No_N_imp_practice_Tech = pd.read_excel('../../AMP_No_Nitrogen_import_practice.xlsx','Vol 1 AMP No N. Tech. Reps')

print(No_N_imp_practice_full.head())  #printing just the first 5 rows of each column in  No_N_imp file 

time = No_N_imp_practice_Tech ['Time(days)']
No_N_A1 = No_N_imp_practice_Tech ['NoN-A.1']
print(time) #printing all of time column 


# model time array
delt,ndays = 1.0/24.0,40.0
times = linspace(0,ndays,int(ndays/delt)) # model time-series

# model parameters
alpha = 0.00000009    #Producer's affinity for the resource
delta = 0.12     #0.000007      #loss rate (1/day) of producer 

# steady-states 
Rstar = delta/alpha

# initial conditions
R = 5e+6   #Had to up the starting nutrient pool concentration from (1e+5) by a fair amount in order to get the P to get the a high enough peak to match the data
P = No_N_A1[0]

# biomass / water arrays
Rs,Ps = r_[[]],r_[[]]
 #print (Rs)
for t in times:
	Rs = append(Rs,R)
	Ps = append(Ps,P)
	   #print (P,R)
	dRsdt = -alpha*R*P
	dPsdt = alpha*R*P - delta*P
	R = max(R + dRsdt*delt,1e-4) # numerical integration ('Euler integration')
	P = max(P + dPsdt*delt,1e-4)

#graphing
fig, axis1 = plt.subplots()   #making multiple plots that will all show at the same time when show() command is given
axis1.set_xlabel('Time (days)')  #labeling the x axis of the first subplot
axis1.set_ylabel('Cell Density (cells per mL)') #labeling the y axis of the first subplot
semilogy()
axis1.plot(times,Ps,c='g',label='Producers') #plotting Ps on first subplot
plt.plot(time,No_N_A1, marker='x', c='b', label='No N A1') #plotting points with a line going through them; print    ing on same axes as the scatter plot graph

#axis 1 subplot legend
legend = plt.legend(loc='upper right')   #, bbox_to_anchor=(40,1e+7))
legend.draw_frame(False)

axis2 = axis1.twinx()   #making it so the two subplots sit on top of one antoher by twining the x axis

axis2.set_ylabel('Nutrient Concentration') #labeling the y axis of the second subplot
semilogy()
#axis2.axhline(Rstar,c='r',ls='--',label='Rstar')
axis2.plot(times,Rs,c='r',label='Resource')

#axis 2 subplot legend
legend = plt.legend(loc='lower right')    #, bbox_to_anchor=(0.5,1.5))
legend.draw_frame(False)

plt.title('No N A1 Raw Data on d/a Model')   #labels entire figure with titl

show()

print ('Done!!!')
