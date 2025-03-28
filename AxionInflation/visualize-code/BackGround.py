import matplotlib.pyplot as plt

t_list   = []
f_list   = []
df_list  = []
var_list = []
ske_list = []
kur_list = []

a_list   = []
N_list   = []
H_list   = []

E_list    = []
eps0_list = []
eps_list  = []
eta0_list = []
eta_list  = []
G_list    = []

xi_list   = []
co1_list  = []
co2_list  = []
co1L_list = []
co2L_list = []


# File Loading (float in python = double in C++)
#---------------------------------------------------------------
path  = "/Users/tsuchitashunsuke/Dropbox/lattice/AxionInflation"
path += "/axion-U(1)"
path += "_adjust"
path += "_coulomb"
path0 = path + "/python.dat"
path1 = path + "/background_0.dat"
path2 = path + "/sf_0.dat"
path3 = path + "/energy_0.dat"
path4 = path + "/coupling_0.dat"

# python.dat
num = 0
f0 = open(path0,"r")
for line in f0:
    data = line[:-1].split(' ')
    if(num==0):
        sbackground = int(data[6])
        sexpansion  = int(data[7])
        senergy     = int(data[8])
        sgauge      = int(data[9])
        num += 1
    

# background_0.dat
if(sbackground):
    f1 = open(path1, "r")
    for line in f1:
        data = line[:-1].split(' ')
        t_list.append(float(data[0]))
        f_list.append(float(data[1]))
        df_list.append(float(data[2]))
        var_list.append(float(data[3]))
        ske_list.append(abs(float(data[4])))
        kur_list.append(abs(float(data[5])))

# sf_0.dat
if(sexpansion):
    f2 = open(path2, "r")
    for line in f2:
        data = line[:-1].split(' ')
        a_list.append(float(data[1]))
        H_list.append(float(data[2]))
        N_list.append(float(data[3]))

# energy_0.dat
if(senergy):
    f3 = open(path3, "r")
    for line in f3:
        data = line[:-1].split(' ')
        E_list.append(abs(float(data[6])))
        eps0_list.append(float(data[7]))
        eps_list.append(float(data[8]))
        eta0_list.append(float(data[9]))
        eta_list.append(float(data[10]))
        G_list.append(float(data[11]))

# coupling.dat
if(senergy):
    f4 = open(path4,"r")
    for line in f4:
        data = line[:-1].split(' ')
        xi_list.append(float(data[1]))
        co1_list.append(float(data[2]))
        co2_list.append(float(data[3]))
        co1L_list.append(abs(float(data[4])))
        co2L_list.append(float(data[5]))


# Graph Drawing
#------------------------------------------------------------------
fig = plt.figure(figsize=(10,7))
fig.suptitle('Background')

if(sbackground):
    # Upper Left (Means of Inflaton Field)
    plt.subplot(2,2,1)
    plt.plot(N_list,f_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\phi [M_{pl}]$')

    # Upper Right (Means of Inflaton Velocity)
    plt.subplot(2,2,2)
    plt.plot(N_list,df_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\dot{\phi}/m [M_{pl}]$')

# Lower Left (Hubble Parameter)
plt.subplot(2,2,3)
plt.plot(N_list,H_list, marker='.', markersize='1', linestyle='None')
plt.xlabel(r'$N_e$')
plt.ylabel(r'$H/m$')

if(senergy):
    # Lower Rithg (Energy Consistency Condition)
    plt.subplot(2,2,4)
    plt.plot(N_list,xi_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\xi$')

plt.grid(True)
plt.show()


# Slow-Roll Parameters
if(senergy):
    fig = plt.figure(figsize=(10,4))
    fig.suptitle('Slow-Roll Parameters')

    plt.subplot(1,2,1)
    plt.plot(N_list,eps0_list, marker='.', markersize='2', linestyle='None')
    plt.plot(N_list,eps_list, marker='.', markersize='0.5', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\epsilon$')
    plt.legend(['Homogeneous','+ Space depend, Gauge'])

    plt.subplot(1,2,2)
    plt.plot(N_list,eta0_list, marker='.', markersize='2', linestyle='None')
    plt.plot(N_list,eta_list, marker='.', markersize='0.5', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\eta$')
    plt.legend(['Homogeneous','+ Space depend, Gauge'])

    plt.grid(True)
    plt.show()


# Consistency Condition
if(senergy):
    fig = plt.figure(figsize=(10,4))
    fig.suptitle('Consistency Condition')

    plt.subplot(1,2,1)
    plt.plot(N_list,E_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\|E-1\|$')
    plt.yscale('log',base=10)

    if(sgauge):
        plt.subplot(1,2,2)
        plt.plot(N_list,G_list, marker='.', markersize='1', linestyle='None')
        plt.xlabel(r'$N_e$')
        plt.ylabel(r'$G$')
        plt.yscale('log',base=10)

    plt.grid(True)
    plt.show()


# Backreaction
if(senergy):
    fig = plt.figure(figsize=(10,7))
    fig.suptitle('Backreaction Condition')

    # Upper Left (back reaction条件1---厳しい方)
    plt.subplot(2,2,1)
    plt.plot(N_list,co1_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$H^2\xi^{-3/2}e^{\pi\xi}/26\pi\dot{\phi}$')

    # Upper Right (back reaction条件2---緩い方)
    plt.subplot(2,2,2)
    plt.plot(N_list,co2_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$H\xi^{-3/2}e^{\pi\xi}/146$')

    # Lower Left 
    plt.subplot(2,2,3)
    plt.plot(N_list,co1L_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$|BR/\frac{dV}{d\phi}|$')

    # Lower Right (back reaction条件2---緩い方)
    plt.subplot(2,2,4)
    plt.plot(N_list,co2L_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\rho_{\mathrm{GF}}/3H^2M_{\mathrm{pl}}^2$')

    plt.grid(True)
    plt.show()


# Cumulant
if(sbackground):
    fig = plt.figure(figsize=(6,6))
    fig.suptitle('Skewness(3rd) and Kurtosis(4th)')

    plt.plot(N_list,ske_list, marker='.', markersize='1', linestyle='None')
    plt.plot(N_list,kur_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.yscale('log')
    plt.legend([r'$\frac{|<\delta\phi^3>|}{\sigma_{\delta\phi}^3}$'
                ,r'$\frac{<\delta\phi^4>}{\sigma_{\delta\phi}^4}-3$'])

    plt.grid(True)
    plt.show()