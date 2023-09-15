from scipy import *
from pylab import *
from numpy import *
from matplotlib import *
import sys

rcParams["legend.loc"]


#create figures
f,ax = subplots()

#model time array (sampling for n days with delta being the # of days between each sample?)
delt,ndays = 1.0/24.0,100.0
times = linspace(0,ndays,int(ndays/delt)) #model time series (start at 0, end with ndays, number of steps =  ndays/delta)

#model parameters
alpha = 2e-07  #Producer's (P) affinity for Resource (R)
delta = 0.00002   #loss rate (1/day) of Producer (P)
Vmax = 0.5   #max velocity of enzyme being used for consumption of R

#steady states
Rstar = (delta*Vmax)/(Vmax-(alpha*delta))

#initial conditions
R = 1e+5
P = 1e+1

#biomass / water arrays
Rs,Ps, = r_[[]],r_[[]]
 

for t in times:
        Rs = append(R,Rs) #eight spaces over used for indention 
        Ps = append(P,Ps) 
        dRst = -(R*Vmax/(R+(Vmax/alpha)))*P
        dPst = ((R*Vmax)/(R+(Vmax/alpha)))*P - delta*P
        R = max(R + dRst*delt,1e-4) #numerical integration (euler's method)
        P = max(P + dPst*delt,1e-4) 

#graphing
#axhline(Rstar,c='r',ls='--',label='Rstar')
  
plot(times,Rs,c='g',label='Producers') #are wrong somehow...--_('_')_-- ...had to switch labels 
plot(times,Ps,c='r',label='Resources')
#figure.legend('lower left')........hoq to change location of fig legend? 
ax.set_ylabel('Concentration')
ax.set_xlabel('Time')

l = ax.legend()
l.loc='lower left'
l.draw_frame(False)
plt.semilogy()     #logs y axis 
show()
