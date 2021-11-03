'''
40nM_N_A1
Fitting Ben's grpwth data to model via met alg.
created by DKM
Location: /Users/dkm/Documents/Talmy_research/Zinser\ and\ Ben/Scripts/
'''


from scipy.integrate import *
from scipy import *
from pylab import *
from numpy import *
from matplotlib import *
import sys
import xlrd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 



############################################

#   Importing the data 

############################################


#No_N_imp_practice_full = pd.read_excel("../../AMP_No_Nitrogen_import_practice.xlsx", header = 0)    #had to put 2 dots and 1 / to go back a single folder destination  
Ben_data_import = pd.read_excel('../Data/Ben_data_import_fixed.xlsx', header=0)

#No_N_imp_practice_Tech = pd.read_excel('../../AMP_No_Nitrogen_import_practice.xlsx','Vol 1 AMP No N. Tech. Reps')

40nM_N_A1 = Ben_data_import ['40nM - A.1']      #nameing the technical replicate Producer cell density values as 40nM_N_A1
40nM_N_A1 = array((Ben_data_import['40nM - A.1'])) #turining 40nM_N_A1 into an array
Ps = 40nM_N_A1 #Cells per mL of Pro from data array; renaming it Ps
Q = (9.5*1.0e-15)*14*1.0e06       #Nitrogen CEll quota from literature (9.5fg/cell)...multiply by molar wieght (14)  and 10^6 to get micrograms
Rs = r_[[((((max(40nM_N_A1))-40nM_N_A1[0])*(Q)),0)]] #making an array for Rs values, just the amount of R calculated from taking the peak cell growth and finding change in cell density and finding resource needed for that cell increase based on cell Nitrogen quota



# Standard deviation of data (should pull from data sets after I start using more of the file

Rsd = np.std(Rs)
Psd = np.std(Ps)




# Sample times

dtimes = array(Ben_data_import ['Time(days)'])   #data times; put into an array from time column in dataframe  
ptimes = dtimes      #data times; put into an array from time column in dataframe  
rtimes = dtimes[[0,12]]   #r data only at first and last time point so rtiems is just 2 of the positions of dtimes
#should rtimes be different? should they correspond to the full time series even if we have no resource concentration data points for any of those?  



############################################

#   Integration function

###########################################
Rsd = np.std(Rs)
Psd = np.std(Ps)




# function to be integrated

def f(u,t, alpha, delta, Vmax):      #making a function called 'f' that will run the given equations and spit out arrays of dRst and dPst when the values for R,P,t,alpha,delta,and Vmax are given?
        R,P = u[0],u[1]
        dRst = -(R*Vmax/(R+(Vmax/alpha)))*P*Q
        dPst = ((R*Vmax)/(R+(Vmax/alpha)))*P - delta*P
        #Rstar = (delta*Vmax)/(Vmax-(alpha*delta))  #should Rstar be here or no?n
        #Rs,Ps, = r_[[]],r_[[]]
        return concatenate([r_[[dRst]],r_[[dPst]]])


#   calling odeint

def integrate(params,inits,rtimes,ptimes,forshow=False,delt=1.0 / 24.0):
        ndays = 42.0
        mtimes = linspace(0,ndays,int(ndays/delt)) #mtimes = model times. telling m times to be an even number of time points starting at 0 and ending at ndays taking ndays/delt number of steps (can't can 'delt' because it is 1/24th of a day...which rounds to 0...so it gets confused if you do that...in linspace.
        alpha, delta, Vmax = exp(params[0]),exp(params[1]), exp(params[2])  #exponentiated so as to search space correctly? part of chi calculation? or for a different reason?  
        u = odeint(f,inits,mtimes,args=(alpha, delta, Vmax))
        if forshow==False:
                Rinds = r_[[where(abs(a-mtimes)==min(abs(a-mtimes)))[0][0] for a in rtimes]]  #saying that Rinds is looking for the time slices in mtimes (which is the time array that has nday/delt steps even taken inbetween 0 and nedays)...(a just means every value in an array...so searching all possibilities in the array)  that are closest to the time points taken in the data. Use absoluet value so the closes time could be a little behind or ahead of actual time point. This just tells you the slot number inthe array in mtimes, not the actual time slice value...so you then multiply the array slot number by the 'u' array that holds the values of the P or R model values at said time slice slot that was chosen. (below) ...this lets you come up with the nutritent concentration value or producer cell number/ml that would be given by the model at that time slice value.   
                Pinds = r_[[where(abs(a-mtimes)==min(abs(a-mtimes)))[0][0] for a in ptimes]]
                #print(u.T[0])
                Rnt = u.T[0][Rinds]
                Pnt = u.T[1][Pinds]

        else:
                Rnt = u.T[0]
                Pnt = u.T[1]

        return Rnt,Pnt



#################################################### 

#  generic arrays and optimization paramerts

####################################################


stds = zeros(3) + .05 # this controls how large the random steps are in the parameter search (see below)

opt = r_[[1,1,1]] # this allows you to control whether each parameter is imposed or fittedi

names = ['alpha', 'delta', 'Vmax'] # name each parameter array - useful for printing the output


nits = 1000 # number of iterations

pits = 100  # frequency with which to print results to the user about the progress of the algorithm

burnin = 100 # burnin period - don't 'store' the results for the first number of iterations as these are likely to be very far from the optimal ssolutions


##################################

# MODEL FITTING: Phytoplankton abundance 

##################################


## set up first guess params, MHA parameters, etc. 

alpha = 0.00000029/Q  #Producer's (P) affinity for Resource (R)
delta = 0.08   #loss rate (1/day) of Producer (P)
Vmax = 0.9   #max velocity of enzyme being used for consumption of R



# put in arrays for ease manipulating

params = r_[[alpha,delta,Vmax]] # put inside a single array ready for the algorithm

params = log(params) # we do the parameter search in log space (I will explain later)

npars = params.shape[0] # number of parameters being searched through

40nM_N_A1 = array((Ben_data_import['40nM - A.1']))
Ps = 40nM_N_A1
Rs = r_[[((((max(40nM_N_A1))-40nM_N_A1[0])*(Q)),0)]]     #should the endpoint be somthing other than 0? 


# initial conditions
  #inits = r_[[Rs, Ps]]
inits = r_[[(((max(40nM_N_A1))-40nM_N_A1[0])*(Q)), 40nM_N_A1[0]]]     #what initial conditions do I put here??? For R and P? 
#R = 6.9e+6 <--- was a guess...replace with ((P cell density at peak of growth) - (P cell density initital))*Nitrogen Quota for cells from a paper 
#P = 40nM_N_A1[0]   #makes the starting value of P be the first data point collected for P in Excel data file

mtimes = linspace(0,42.0,int(42.0*24.0))   #put this down here bc you have to graph model R outputs by mtomes to get the model to actually be time resolved on the graph


# first run just to get error

Rnt,Pnt = integrate(params,inits,rtimes,ptimes,forshow=True)    #shou

#chi =   sum((Rnt - Rs) ** 2 / (Rsd ** 2)) +  sum((Pnt - Ps) ** 2 / (Psd ** 2)) 
#this is the equation for Error summ of swuares      #how to make chi for my parameters? 


#print('Chi = ' ,chi)





#######################
#graphing
############################


#delt = 1.0/24.0
#mtimes = linspace(0,ndays,int(ndays/delt))



f1,(ax1,ax2)  = plt.subplots(1,2)
ax2.plot(mtimes,Pnt,c='g',label='Model Fit')
ax2.scatter(ptimes,Ps,marker='x',c='g',label='Producers')
ax2.set_xlabel('Time (Days)')
ax2.set_ylabel('Cell Density (Cells per mL)')
ax1.plot(mtimes,Rnt,c='r',label='Model Fit')
ax1.scatter(rtimes,Rs,marker='x',c='r',label='Resources')
ax1.set_xlabel('Time (Days)')
ax1.set_ylabel('Nitrogen Concentration (micromoles per mL)')



ax1.legend(loc='best')
ax2.legend(loc='best')
show()






#distribution arrays and aceptance ratios - these are contaioners to be added to

ar = 0.0

ars = r_[[]]

alphas, deltas, Vmaxes = r_[[]],r_[[]],r_[[]]

pall = [alphas, deltas, Vmaxes]  #what is pall???  .just 'parameters,all?


#now actually do the fitting 

for it in arange(1,nits,1):  #making an evenly spaces arrangement of values from 1 to nits (number of itterations) on the step scale of 1.  

        parsnew = params + opt*normal(0,stds,npars)  #this is where we randomly change the parameter values   #parsnew = parameters new

        R, P = (((max(40nM_N_A1))-40nM_N_A1[0])*(Q)), 40nM_N_A1[0]  #have to reassign initial values for R and P because we are now opperating in a loop

        inits = r_[[R,P]]
        Rnt,Pnt = integrate(parsnew,inits,rtimes,ptimes)   #call the integrate fucntion
        chinew = sum((Rnt - Rs) ** 2 / (Rsd **2)) + sum ((Pnt - Ps) ** 2 / (Psd ** 2))   #calculare the error
        if exp(chi-chinew) > rand(): #key step

                chi = chinew

                params = parsnew #new paramenters can be a bit 'wrong' ...tjis is what lets us search multiple cariations of the 'right' answer

                if it > burnin: #only letting you store variables if they occur after the brunin period. 
                        pall = append(pall,params[:,None],1)  #I don;t understand what this stuff in the [] is doing there   #pall = all parameters
                ar = ar + 1.0  # number of accepted results (NOT the acceptance ratio until you diveide by total number of interations)
        if (it % pits == 0):   #pits = number of p itterations
                print(it,chi,ar/pits)   # printing the iteration number, the chi value, and the acceptance ration which is the total number of accepted values (ar) divided by the total number of iterations (pits)
                ars = append(ars,ar/pits) #making an array of the aceptance ratios
                ar = 0.0 



#print output to the screen 

print(' ')

print ('Optimal Parameters for Phytoplankton Abundance')
pars = r_[[mean(p) for p in pall]]   #array of the average of each parameter in 'pall' (all parameters) array
for (p,l) in zip(pars,names):        #putting array of parameter values and parameter names together
        print(l,'=',exp(p))   #printing the name of the parameter and the parameter value given from algorithum


print(' ')

print('Standard Deviations')
for (p,l) in zip(pall,names):
        print(l+'std','=',std(exp(p)))

print(' ')




#redefine times for nicer looking plots

delt = 1.0/24.0

ftimes = linspace(0,amax(rtimes)/24.0, ((amax(ptimes) / delt)+1)*24.0) #should I be using rtimes or ptimes?    #Making a times array specifically for the figure so that the curve looks better

n = ftimes.shape[0]




# run again just for nicer looking plots (more even timesteps)

R,P = ((((max(40nM_N_A1))-40nM_N_A1[0])*(Q)), 40nM_N_A1[0])

inits = r_[[R,P]]

print('Standard deviations')

for (p,l) in zip(pall,names):

        print(l+'std','=',std(exp(p)))

print(' ')

Rnt, Pnt = integrate(pars,inits,ftimes,ftimes,forshow=True,delt=delt)    #can't graph this stuff when forshow = true because Pnt and Rnt end up being much smaller than ftimes (size of each array is different) so you can't graph then against each other...idk why changing 'forshow' argumnet to 'true' does this. :/ 




# plot

#ax1[0].errorbar(times,R/1e+6,yerr=Rsd/1e+6,c='maroon',fmt='o',label='R')

#ax1[1].errorbar(times,P/1e+8,yerr=Psd/1e+8,c='g',fmt='o',label='P')

#ax1[0].plot(times,Rnt/1e+6,c='maroon',lw=1.5,label='model fit')

#ax1[1].plot(times,Pnt/1e+8,c='g',lw=1.5,label='model fit')



#plotting parameters? 

ax2[0].hist(exp(pall[0]),label='alpha',color='maroon')
ax2[1].hist(exp(pall[1]),label='delta',color='maroon')
ax2[2].hist(exp(pall[2]),label='Vmax',color='maroon')

"""

#ax2[0].set_xlabel(r'Host growth rate (day $^{-1}$)',fontsize=fs)

#ax2[1].set_xlabel(r'Clearance rate (pl day$^{-1}$)',fontsize=fs)

#ax2[2].set_xlabel(r'lysis rate (day $^{-1}$)',fontsize=fs)

#ax2[3].set_xlabel(r'Burst size',fontsize=fs)



#ax2[0].set_ylabel(r'Frequency',fontsize=fs)

#ax2[1].set_ylabel(r'Frequency',fontsize=fs)

#ax2[2].set_ylabel(r'Frequency',fontsize=fs)

#ax2[3].set_ylabel(r'Frequency',fontsize=fs)



ax1[1].set_ylim([0,2])



l1 = ax1[0].legend(loc='upper right',prop={'size':14})

l1.draw_frame(False)

l3 = ax2[0].legend(loc='upper center',prop={'size':14})

l3.draw_frame(False)



ax1[0].text(0.07,0.9,'a',ha='center',va='center',color='k',transform=ax1[0].transAxes)

ax1[1].text(0.07,0.9,'b',ha='center',va='center',color='k',transform=ax1[1].transAxes)



#f1.savefig('temp_dynamics')

#f2.savefig('temp_params')

"""

show()

print('Done')



































