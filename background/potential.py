import matplotlib.pyplot as plt
import numpy as np

# Parameters
#----------------------------------------------------
coeff = 0.966
couf = 1./3.3

# Inflaton field
#----------------------------------------------------
phi = np.linspace( 0.01, 10.0, 1000)

# Inflaton potential
#----------------------------------------------------
#V = 0.5*pow(phi,2.)+coeff*couf*phi*np.sin(phi/couf)
V = np.exp(-phi*phi)*phi*phi*phi
plt.plot(phi, V)
plt.yscale('log')
plt.xscale('log')
plt.show()

# Others
#-------------------------------------------------------
'''
phi=3.3
hubble=(phi**2.)/((phi**2.)+0.25)-A*(phi-phi0)/( 1. + (((phi-phi0)**2.)/(2.*(sigma**2.))) )
print(hubble)

dV=0.5*phi/(((phi**2.)+0.25)**2.)-A*(1.-0.5*((phi-phi0)**2.)/(sigma**2.))/((1.+0.5*((phi-phi0)**2.)/(sigma**2.))**2.)
print(dV/(3.*hubble))
print(dV/(3.*(hubble**2.)))'''
