import matplotlib.pyplot as plt

t_list    = []
f_list    = []
df_list   = []
var_list  = []
ske_list  = []
kur_list  = []

a_list    = []
N_list    = []
H_list    = []

E_list    = []
eps0_list = []
eps_list  = []
eta0_list = []
eta_list  = []


# File Loading (float in python = double in C++)
#---------------------------------------------------------------
path  = "/Users/tsuchitashunsuke/Dropbox/lattice/SingleInflation"
path += "/SFInflation"
path0 = path + "/python.dat"
path1 = path + "/means0_0.dat"
path2 = path + "/sf_0.dat"
path3 = path + "/energy_0.dat"

# python.dat
num = 0
f0 = open(path0,"r")
for line in f0:
    data = line[:-1].split(' ')
    if(num==0):
        sbackground = int(data[4])
        sexpansion  = int(data[5])
        senergy     = int(data[6])
        num += 1

# means_0.dat
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
        E_list.append(abs(float(data[5])))
        eps0_list.append(float(data[6]))
        eps_list.append(float(data[7]))
        eta0_list.append(float(data[8]))
        eta_list.append(float(data[9]))


# Graph Drawing
#------------------------------------------------------------------
fig = plt.figure(figsize=(10,7))
fig.suptitle('Background')

# Upper Left (Means of Inflaton Field)
if(sbackground):
    plt.subplot(2,2,1)
    plt.plot(N_list,f_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\phi [M_{pl}]$')

    # Upper Right (Means of Inflaton Volocity)
    plt.subplot(2,2,2)
    plt.plot(N_list,df_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\dot{\phi}/m [M_{pl}]$')

# Lower Left (Hubble Parameter)
plt.subplot(2,2,3)
plt.plot(N_list,H_list, marker='.', markersize='1', linestyle='None')
plt.xlabel(r'$N_e$')
plt.ylabel(r'$H/m$')

# Lower Right (Energy Consistency Condition)
if(senergy):
    plt.subplot(2,2,4)
    plt.plot(N_list, E_list, marker='.', markersize='1', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\|E-1\|$')
    plt.ylim([pow(10,-16),pow(10,-4)])
    plt.yscale('log', base = 10)

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
    plt.ylabel(r'$\epsilon$',rotation=0,fontsize=15)
    #plt.ylim([0.009,0.011])
    plt.legend(['Homogeneous','+ Space Dependence'])

    plt.subplot(1,2,2)
    plt.plot(N_list,eta0_list, marker='.', markersize='2', linestyle='None')
    plt.plot(N_list,eta_list, marker='.', markersize='0.5', linestyle='None')
    plt.xlabel(r'$N_e$')
    plt.ylabel(r'$\eta$',rotation=0,fontsize=15)
    plt.legend(['Homogeneous','+ Space Dependence'])

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