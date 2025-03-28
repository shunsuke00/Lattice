import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

pi = math.pi

# File Loading (float in python = double in C++)
#---------------------------------------------------------------
path = "/Users/tsuchitashunsuke/Dropbox/lattice/AxionInflation"
#path += "/data_new/weak/256"
path += "/axion-U(1)"
path += "_adjust"
path += "_coulomb"
path0 = path + "/python.dat"
path1 = path + "/spectratimes_0.dat"
path2 = path + "/spectra_ax_0.dat"
path3 = path + "/spectra_eff_ax_0.dat"


# python.dat
num = 0
f0 = open(path0, "r")
for line in f0:
    data      = line[:-1].split(' ')
    if(num==0):
        num_point = float(data[3])
        L_pr      = float(data[4])
        m         = float(data[5])
        num += 1
    else:
        fi        = float(data[1])
        Hi        = float(data[2])
        af        = float(data[3])
        adf       = float(data[4])
        fdf       = float(data[5])

dx_pr = L_pr/num_point
# Coefficient of Dimensionless Power Spectra
coeff = pow(pow(dx_pr,2)/L_pr,3)/(pow(pi,2)*2.)


# spectratimes_0.dat
zeta_factor = []
xi_list     = []
N_list      = []
f1 = open(path1, "r")
for line in f1:
    data = line[:-1].split(' ')
    a    = float(data[1])
    ad   = float(data[2])
    fd   = float(data[3])
    zeta_factor.append(pow(ad/(a*fd),2))
    xi_list.append(float(data[4]))
    N_list.append(round(float(data[5]),1))


# Some Setting
taketime = 13
k_lat = [[] for _ in range(taketime)]
f_lat = [[] for _ in range(taketime)]
k_the = [[] for _ in range(taketime)]
f_the = [[] for _ in range(taketime)]
k_ori = [[] for _ in range(taketime)]
f_ori = [[] for _ in range(taketime)]
index_list = [[] for _ in range(taketime)]

length = math.floor(math.sqrt(3)*num_point/2.)+2
for i in range(taketime):
    index_list[i] = i


# spectra_ax_0.dat
f2 = open(path2, "r")
num = 0
for line in f2:
    num += 1
    for index in range(taketime):
        if num>index_list[index]*length and num<index_list[index]*length+length:
            data = line[:-1].split(' ')
            k_pr = float(data[2])
            keff_pr = float(data[3])
            fkeff2_pr = float(data[4])
            if(fkeff2_pr!=0):
                k_lat[index].append(k_pr)
                f_lat[index].append(fkeff2_pr*pow(k_pr,3)*coeff*zeta_factor[index_list[index]])
                k_the[index].append(keff_pr)
                f_the[index].append(fkeff2_pr*pow(keff_pr,3)*coeff*zeta_factor[index_list[index]])


# spectra_eff_ax_0.dat
f3 = open(path3, "r")
num = 0
for line in f3:
    num += 1
    for index in range(taketime):
        if num>index_list[index]*length and num<index_list[index]*length+length:
            data = line[:-1].split(' ')
            keff_pr = float(data[2])
            fkeff2_pr = float(data[3])
            if(fkeff2_pr!=0):
                k_ori[index].append(keff_pr)
                f_ori[index].append(fkeff2_pr*pow(keff_pr,3)*coeff*zeta_factor[index_list[index]])


# Analytical Solution
t      = np.linspace(0.1, 400, 1000)
p_vac  = pow(2.55*pow(adf/af,2)/(pi*fdf),2)*pow(10,-12)
xi_ini = xi_list[0]
#p_cuv1 = p_vac + pow(p_vac,2)*7.5*pow(10,-5)*pow(xi,-6)*math.exp(4*pi*xi)
p_cuv1 = p_vac + pow(p_vac,2)*3.*pow(10,-5)*pow(xi_ini,-5.4)*math.exp(4*pi*xi_ini)
xi_end = xi_list[taketime-1]
#p_cuv2 = p_vac + pow(p_vac,2)*7.5*pow(10,-5)*pow(xi,-6)*math.exp(4*pi*xi)
p_cuv2 = p_vac + pow(p_vac,2)*3.*pow(10,-5)*pow(xi_end,-5.4)*math.exp(4*pi*xi_end)

y = math.pow(10,-6)*np.ones_like(t)


# Graph Drawing
#------------------------------------------------------------------
legend_labels = [r"$N_e=$" + str(N_list[index_list[i]]) for i in range(taketime)]
legend_labels.append(r"$P_{vac}$")
legend_labels.append(r"$P_{\mathcal{R}}$")

# Latticeeasy Method
fig = plt.figure(figsize=(8,6))
plt.suptitle('Curvature Power Spectram (latticeeasy)')

for i in range(taketime):
    plt.plot(k_lat[i],f_lat[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.5',linestyle='solid',linewidth='1.5')
plt.axhline(y=p_vac, color='black', linestyle='-', linewidth='1')
plt.fill_between(t,p_cuv1, p_cuv2, color='gray', alpha=0.4)

plt.legend(legend_labels)

plt.xlabel(r'$k_{eff}/m$')
plt.ylabel(r'$P_{\zeta}$')
plt.xlim([5,200])
#plt.ylim([5*pow(10,-10),pow(10,-6)])
plt.xscale('log', base = 10)
plt.yscale('log', base = 10)

plt.grid(True)
plt.show()


# Angelo Method
fig = plt.figure(figsize=(8,6))
plt.suptitle('Curvature Power Spectram (thesis)')

for i in range(taketime):
    plt.plot(k_the[i],f_the[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.5',linestyle='solid',linewidth='1.5')
plt.axhline(y=p_vac, color='black', linestyle='-', linewidth='1')
plt.fill_between(t,p_cuv1, p_cuv2, color='gray', alpha=0.4)

plt.legend(legend_labels)

plt.xlabel(r'$k_{eff}/m$')
plt.ylabel(r'$P_{\mathcal{R}}$')
plt.xlim([3,200])
#plt.ylim([5*pow(10,-10),pow(10,-6)])
plt.xscale('log', base = 10)
plt.yscale('log', base = 10)

plt.grid(True)
plt.show()


# Original Method
fig = plt.figure(figsize=(8,6))
plt.suptitle('Curvature Power Spectram (original)')

for i in range(taketime):
    plt.plot(k_ori[i],f_ori[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.5',linestyle='solid',linewidth='1.5')
plt.axhline(y=p_vac, color='black', linestyle='-', linewidth='1')
plt.fill_between(t,p_cuv1, p_cuv2, color='gray', alpha=0.4)

plt.legend(legend_labels)

plt.xlabel(r'$k_{eff}/m$')
plt.ylabel(r'$P_{\zeta}$')
plt.xlim([5,200])
#plt.ylim([5*pow(10,-10),pow(10,-6)])
plt.xscale('log', base = 10)
plt.yscale('log', base = 10)

plt.grid(True)
plt.show()
