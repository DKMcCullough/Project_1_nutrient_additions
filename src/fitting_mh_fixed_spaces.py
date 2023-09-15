from scipy.integrate import *

from numpy import *

from pylab import *

from scipy import *

from averagedata import average_data

seed(0) # seed random number generator to have reproducible results

#####################################################

# integration functions

#####################################################



# function to be integrated

def f(u,t,phi,muh,lambd,beta,loss):

        S,I,V = u[0],u[1],u[2] 

        dSdt = muh*S - phi*S*V

        dIdt = phi*S*V-lambd*I

        dVdt = beta*lambd*I - phi*(S+I)*V - loss*V

        return concatenate([r_[[dSdt]],r_[[dIdt]],r_[[dVdt]]])



# calling odeint

def integrate(params,inits,htime,vtime,forshow=False,delt=900.0 / 86400.0):

        days = amax(htime)

        times = arange(0,days,delt)

        phi,muh,lambd,beta,loss = exp(params[0]),exp(params[1]),exp(params[2]),exp(params[3]),0.1 # parameters

        u = odeint(f,inits,times,args=(phi,muh,lambd,beta,loss)).T

        if forshow==False:

                hinds = r_[[where(abs(a-times)==min(abs(a-times)))[0][0] for a in htime]] # get time indices where host abundances were measured

                vinds = r_[[where(abs(a-times)==min(abs(a-times)))[0][0] for a in vtime]] # same for viruses

                hnt = u[0][hinds] + u[1][hinds] # host density

                vnt = u[2][vinds] # virus density

        else:

                hnt = u[0] + u[1]

                vnt = u[2]

        return hnt,vnt



#####################################################

# load in the data

#####################################################



# sample times

htimes = r_[[0.83, 4.52, 8.67, 16.74, 24.63, 32.70, 40.39, 48.47, 56.11, 60.44, 64.37, 72.21, 83.85, 92.32]]/24.0

vtimes = r_[[0.83, 4.52, 8.67, 16.74, 24.63, 32.70, 40.39, 48.47, 56.11, 60.44, 64.37, 72.21, 83.85, 92.32]]/24.0



# hosts

# htemp - 25 degrees C

# ltemp - 20 degrees C

hhtemp = r_[[4.9371, NaN, NaN, NaN, 5.1084, 5.0396, 4.4043, 4.2937, 4.106, 4.1491, 4.0944, NaN, NaN, 4.0139]]

hltemp = r_[[4.8594, 4.8441, 4.9695, 4.8373, 4.8692, 4.7605, 4.7689, 4.793, 4.3327, 4.2002, 4.0521, 3.6699, 3.7099, 3.6481]]



hhtemp = ma.masked_invalid(hhtemp) # mask nans

hltemp = ma.masked_invalid(hltemp)



hhtemp = 10**hhtemp # convert from log10

hltemp = 10**hltemp



hhtempsd = ma.mean(hhtemp)*0.1 # these data did not come with standard deviations so i made them up

hltempsd = ma.mean(hltemp)*0.1



# virus

vhtemp = r_[[7.0508, 5.9731, 6.1026, 5.8457, 5.7997, 6.1993, 6.2119, 6.752, 7.5733, 7.7027, 8.0078, 7.8681, 7.5649, 7.5191]]

vltemp = r_[[7.2483, NaN, NaN, NaN,6.0053, 7.486, 8.0018, 7.9091, 8.0051, 8.4369, 8.2288, NaN, NaN, 7.6049]]



vhtemp = ma.masked_invalid(vhtemp) # mask nans

vltemp = ma.masked_invalid(vltemp)



vhtemp = 10**vhtemp # convert from log10

vltemp = 10**vltemp



vhtempsd = ma.mean(vhtemp)*0.1 # made up standard deviations for the virus as well

vltempsd = ma.mean(vltemp)*0.1



#####################################################

# set up figures

#####################################################

f1,ax1 = subplots(1,2,figsize=[9.5,4.0])

f2,ax2 = subplots(4,1,figsize=[6,15])

f1.subplots_adjust(bottom=0.13,wspace=0.3,hspace=0.3)

f2.subplots_adjust(hspace=0.45)



#######################################################

# generic arrays and optimization parameters

#######################################################



stds = zeros(4) + .05 # this controls how large the random steps are in the parameter search (see below)

opt = r_[[1,1,1,1]] # this allows you to control whether each parameter is imposed or fitted

names = ['phi','muh','lam','beta'] # name each parameter array - useful for printing the output



nits = 1000 # number of iterations

pits = 100  # frequency with which to print results to the user about the progress of the algorithm

burnin = 100 # burnin period - don't 'store' the results for the first number of iterations as these are likely to be very far from the optimal ssolutions



##################################

# MODEL FITTING: high temp

##################################



## set up first guess params, MHA parameters, etc. 

phi = 6.57946111322e-08

muh = 1.51155534756

lam = 6.37271315481

beta = 82.9215651713



# put in arrays for ease manipulating

params = r_[[phi,muh,lam,beta]] # put inside a single array ready for the algorithm

params = log(params) # we do the parameter search in log space (I will explain later)

npars = params.shape[0] # number of parameters being searched through



# initial conditions

inits = r_[[hhtemp[0],0,vhtemp[0]]]



# first run just to get error

hnt,vint = integrate(params,inits,htimes,vtimes)

chi =   sum((hnt - hhtemp) ** 2 / (hhtempsd ** 2)) + sum((vnt - vhtemp) ** 2 / (vhtempsd ** 2))        



# distribution arrays and acceptance ratios - these are containers to be added to

ar = 0.0

ars = r_[[]]

phis,muhs,lams,betas = r_[[]],r_[[]],r_[[]],r_[[]]

pall = [phis,muhs,lams,betas]



# now actually do the fitting

for it in arange(1,nits,1):

        parsnew = params + opt*normal(0,stds,npars) # this is where we randomly change the parameter values 

        sus,inf,vir = hhtemp[0],0,vhtemp[0] # have to reassign initial conditions because it's in a loop

        inits = r_[[sus,inf,vir]] # put initial conditions in an array

        hnt,vnt = integrate(parsnew,inits,htimes,vtimes) # call the integration function

        chinew = sum((hnt - hhtemp) ** 2 / (hhtempsd ** 2)) + sum((vnt - vhtemp) ** 2 / (vhtempsd ** 2)) # calculate the error

        if exp(chi-chinew) > rand(): # KEY STEP

                chi = chinew 

                params = parsnew #  new parameters can be a little bit 'wrong'

        if it > burnin: # only store the parameters if you've gone through the burnin period

                pall = append(pall,params[:,None],1)

                ar = ar + 1.0 # acceptance ratio - I can explain this another time

        if (it % pits == 0):

                print(it,chi,ar/pits)

                ars = append(ars,ar/pits)

                ar = 0.0



# print output to screen

print('Optimal parameters for high temp conditions')

pars = r_[[ mean(p) for p in pall]]

for (p,l) in zip(pars,names):

        print(l,'=',exp(p))



print(' ')

print('Standard deviations')

for (p,l) in zip(pall,names):

        print(l+'std','=',std(exp(p)))

print(' ')



# redefine times for nicer looking plots

delt = 900.0 / 86400.0

ftimes = linspace(0,amax(htimes)/24.0,(amax(htimes) / delt)+1)*24.0

n = ftimes.shape[0]



# run again just for nicer looking plots (more even timesteps)

sus,inf,vir = hhtemp[0],0,vhtemp[0]

inits = r_[[sus,inf,vir]]

hnt,vnt = integrate(pars,inits,ftimes,ftimes,forshow=True,delt=delt)



# plot

ax1[0].errorbar(htimes,hhtemp/1e+6,yerr=hhtempsd/1e+6,c='maroon',fmt='o',label='25 C')

ax1[1].errorbar(vtimes,vhtemp/1e+8,yerr=vhtempsd/1e+8,c='maroon',fmt='o',label='25 C')

ax1[0].plot(ftimes,hnt/1e+6,c='maroon',lw=1.5,label='model fit')

ax1[1].plot(ftimes,vnt/1e+8,c='maroon',lw=1.5,label='model fit')



# high temp

ax2[0].hist(exp(pall[1]),label='25 C',color='maroon')

ax2[1].hist(exp(pall[0])*1e+9,label='25 C',color='maroon')

ax2[2].hist(exp(pall[2]),label='25 C',color='maroon')

ax2[3].hist(exp(pall[3]),label='25 C',color='maroon')



##################################

# MODEL FITTING: low temperature - 20 degrees celcius

##################################



# first guess params

phi = 2.92725900355e-08

muh = 1.26863572907

lam = 6.18514303652

beta = 447.232752177



# initial conditions

inits = r_[[hltemp[0],0,vltemp[0]]]



# first run just to get error

hnt,vnt = integrate(params,inits,htimes,vtimes)

chi =   sum((hnt - hltemp) ** 2 / (hltempsd ** 2)) + sum((vnt - vltemp) ** 2 / (vltempsd ** 2))        



# distribution arrays and acceptance ratios

ar = 0.0

ars = r_[[]]

phis,muhs,lams,betas = r_[[]],r_[[]],r_[[]],r_[[]]

pall = [phis,muhs,lams,betas]



# run mha algorithm

for it in arange(1,nits,1):

        parsnew = params + opt*normal(0,stds,npars)

        sus,inf,vir = hltemp[0],0,vltemp[0]

        inits = r_[[sus,inf,vir]]

        hnt,vnt = integrate(parsnew,inits,htimes,vtimes)

        chinew = sum((hnt - hltemp) ** 2 / (hltempsd ** 2)) + sum((vnt - vltemp) ** 2 / (vltempsd ** 2))       

        if exp(chi-chinew) > rand():

                chi = chinew

                params = parsnew

                if it > burnin:

                        pall = append(pall,params[:,None],1)

        ar = ar + 1.0

        if (it % pits == 0):

                print(it,chi,ar/pits)

                ars = append(ars,ar/pits)

                ar = 0.0



print('Optimal parameters for temp limited conditions')

pars = r_[[ mean(p) for p in pall]]

for (p,l) in zip(pars,names):

        print(l,'=',exp(p))



print(' ')

print('Standard deviations')

for (p,l) in zip(pall,names):

        print(l+'std','=',std(exp(p)))

print(' ')



# redefine times for nicer looking plots

delt = 900.0 / 86400.0

ftimes = linspace(0,amax(htimes),(amax(htimes) / delt)+1)

n = ftimes.shape[0]



# run again

sus,inf,vir = hltemp[0],0,vltemp[0]

inits = r_[[sus,inf,vir]]

hnt,vnt = integrate(pars,inits,ftimes,ftimes,forshow=True,delt=900.0 / 86400.0)



# plot

ax1[0].errorbar(htimes,hltemp/1e+6,yerr=hltempsd/1e+6,c='mediumseagreen',fmt='o',label='20 C')

ax1[1].errorbar(vtimes,vltemp/1e+8,yerr=vltempsd/1e+8,c='mediumseagreen',fmt='o',label='20 C')

ax1[0].plot(ftimes,hnt/1e+6,c='mediumseagreen',lw=1.5,label='model fit')

ax1[1].plot(ftimes,vnt/1e+8,c='mediumseagreen',lw=1.5,label='model fit')



# high temp

ax2[0].hist(exp(pall[1]),label='20 C',color='mediumseagreen')

ax2[1].hist(exp(pall[0])*1e+9,label='20 C',color='mediumseagreen')

ax2[2].hist(exp(pall[2]),label='20 C',color='mediumseagreen')

ax2[3].hist(exp(pall[3]),label='20 C',color='mediumseagreen')



# other axes bits

fs = 18

ax1[0].set_ylabel(r'Cells ($\times$10$^6$ml$^{-1}$)',fontsize=fs)

ax1[1].set_ylabel(r'Particles ($\times$10$^{8}$ ml$^{-1}$)',fontsize=fs)

ax1[0].set_xlabel('Time (days)',fontsize=fs)

ax1[1].set_xlabel('Time (days)',fontsize=fs)



ax2[0].set_xlabel(r'Host growth rate (day $^{-1}$)',fontsize=fs)

ax2[1].set_xlabel(r'Clearance rate (pl day$^{-1}$)',fontsize=fs)

ax2[2].set_xlabel(r'lysis rate (day $^{-1}$)',fontsize=fs)

ax2[3].set_xlabel(r'Burst size',fontsize=fs)



ax2[0].set_ylabel(r'Frequency',fontsize=fs)

ax2[1].set_ylabel(r'Frequency',fontsize=fs)

ax2[2].set_ylabel(r'Frequency',fontsize=fs)

ax2[3].set_ylabel(r'Frequency',fontsize=fs)



ax1[1].set_ylim([0,2])



l1 = ax1[0].legend(loc='upper right',prop={'size':14})

l1.draw_frame(False)

l3 = ax2[0].legend(loc='upper center',prop={'size':14})

l3.draw_frame(False)



ax1[0].text(0.07,0.9,'a',ha='center',va='center',color='k',transform=ax1[0].transAxes)

ax1[1].text(0.07,0.9,'b',ha='center',va='center',color='k',transform=ax1[1].transAxes)



f1.savefig('temp_dynamics')

f2.savefig('temp_params')



show()

