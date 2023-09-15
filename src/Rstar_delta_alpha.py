from scipy import *
from pylab import *
from matplotlib import *
import sys

# create figures
f,ax = subplots()

# model time array
delt,ndays = 1.0/24.0,100.0
times = linspace(0,ndays,int(ndays/delt)) # model time-series

# model parameters
alpha = .000002  #Producer's affinity for the resource
delta = 0.00002 #loss rete (1/day) of producer 

# steady-states 
Rstar = delta/alpha


# initial conditions
R = 1e+5
P = 1e+1

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

#axhline(Rstar,c='r',ls='--',label='Rstar')
plot(times,Rs,c='r',label='Resource')
plot(times,Ps,c='g',label='Producers')
ax.set_ylabel('Concentration')
ax.set_xlabel('Time')

l = ax.legend(loc='best')
#l = ax.legend(loc='center left')
#plt.ylim(100,60000000000000)
#plt.xlim(0,300)
plt.semilogy()
show()



