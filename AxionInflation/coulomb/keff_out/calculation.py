import math

pi = math.pi
Mpl_g = 2.17651
Mpl_e = 1.220932
normM = 1./math.sqrt(8.*pi)
c     = 2.99792458
dirac = 6.58211928
Mpc   = 3.085677581
Smass = 1.98855

mOmeg = 0.3062
dOmeg = mOmeg - 0.048
romeg = 4.1577
momeg = 0.1408



# fPBH-beta
'''
value = math.sqrt(3.*Mpl_g/(8.*pi*Smass))*pow(10./(8.*pi*10.75),1./4.)
value = value*3.909/3.363
value = value*Mpl_e*romeg/(0.06*2.3486367)
print(value)
value = value/10.
value = 1./value
print(value)'''
'''
value = math.sqrt(0.6*Mpl_g/(8.*pi*Smass))*pow(10./(8.*pi*106.75),1./4.)
value = value*3.909/3.363
value = value*2*Mpl_e*romeg/(0.1187*2.3486367)
print(value)
value = value/10.
value = 1./value
print(value)'''


# MPBH-T
'''
value = 2.4*math.sqrt(10./106.75)*Mpl_g*Mpl_e*Mpl_e*pow(normM,3.)
value = value/Smass
print(value)'''

# MPBH-k
'''
value = 0.8*pi*pi/(3.*math.sqrt(10)*pow(106.75,1./6.))
value = value*normM*Mpl_g*pow(2.3486367*Mpc/c/dirac,2.)
value = value*pow(3.909,2./3.)
value = value/Smass
print(value)'''


# horizon-problem
'''
value = pow(106.75,1./6.)*pi/math.sqrt(90.)
value = value*pow(3.909,1./3.)
value = value*2.3486367*math.sqrt(8.*pi)/Mpl_e
value = value/(0.678*dirac/Mpc)
print(value)
print(math.log(value*pow(10,25)))'''