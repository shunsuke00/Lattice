import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

pi = math.pi

# File Loading (float in python = double in C++)
#---------------------------------------------------------------
path  = "/Users/tsuchitashunsuke/Dropbox/lattice/SingleInflation"
path += "/SFInflation"
path0 = path + "/python.dat"
path1 = path + "/spectratimes_0.dat"
path2 = path + "/spectra0_0.dat"
path3 = path + "/spectra_eff0_0.dat"

# python.dat
num = 0
f0 = open(path0, "r")
for line in f0:
    data      = line[:-1].split(' ')
    if(num==0):
        model     = data[0]
        num_point = float(data[1])
        L_pr      = float(data[2])
        m         = float(data[3])
        num += 1
    else:
        f0        = float(data[0])
        Hi        = float(data[1])
        af        = float(data[2])
        adf       = float(data[3])

dx_pr = L_pr/num_point
# Coefficient of Dimensionless Power Spectra
coeff = pow(pow(dx_pr,2)/L_pr,3)/(pow(pi,2)*2.)


# spectratimes_0.dat
N_list = []
f1 = open(path1, "r")
for line in f1:
    data      = line[:-1].split(' ')
    N_list.append(round(float(data[1]),1))


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

# spectra0_0.dat
f2 = open(path2, "r")
num = 0
for line in f2:
    num += 1
    for index in range(taketime):
        if num>index_list[index]*length and num<(index_list[index]+1)*length:
            data = line[:-1].split(' ')
            k_pr    = float(data[2])
            keff_pr = float(data[3])
            fkeff2_pr = float(data[4])
            if(k_pr!=0):
                k_lat[index].append(k_pr/Hi)
                f_lat[index].append(fkeff2_pr*pow(k_pr,3)*coeff)
            if(fkeff2_pr!=0):
                k_the[index].append(keff_pr/Hi)
                f_the[index].append(fkeff2_pr*pow(keff_pr,3)*coeff)

# spectra0_eff0.dat
f3 = open(path3, "r")
num = 0
for line in f3:
    num += 1
    for index in range(taketime):
        if num>index_list[index]*length and num<(index_list[index]+1)*length:
            data = line[:-1].split(' ')
            keff_pr = float(data[2])
            fkeff2_pr = float(data[3])
            if(fkeff2_pr!=0):
                k_ori[index].append(keff_pr/Hi)
                f_ori[index].append(fkeff2_pr*pow(keff_pr,3)*coeff)

# Analytical Solution
if model == '0':
    t = np.linspace(1, 200, 1000)
    y = pow(m*adf/(2.*pi*af),2)*pow(t,0) #s=1を仮定している


# Graph Drawing
#------------------------------------------------------------------
legend_labels = [r"$N_e=$" + str(N_list[index_list[i]]) for i in range(taketime)]
legend_labels.append(r"$P_{analytic}=\left(\frac{H}{2\pi}\right)^2$")

# Latticeeasy Method
fig = plt.figure(figsize=(8,6))
plt.suptitle('Latticeeasy method')

for i in range(taketime):
    plt.plot(k_lat[i],f_lat[i],color=cm.cool(i/float(taketime)),marker='.',markersize='1',linestyle='solid',linewidth='2')

if model == '0':
    plt.plot(t,y,color='gray',marker='.',markersize='0.1',linewidth='1')

plt.xlabel(r'$|k_{\rm{eff}}|/H_{\rm{initial}}$')
plt.ylabel(r'$P_{\delta\phi}$',rotation=0,fontsize=15)
plt.legend(legend_labels)

if model == '0':
    plt.xlim([1,200])
    plt.ylim([pow(10,-11),pow(10,-7)])
elif model == '1':
    plt.xlim([1.0,200])
    plt.ylim([pow(10,-12),pow(10,-6)])

plt.xscale('log', base = 10)
plt.yscale('log', base = 10)

plt.grid(True)
plt.show()

# Angelo Method
fig = plt.figure(figsize=(8,6))
plt.suptitle('Angelo method')

for i in range(taketime):
    plt.plot(k_the[i],f_the[i],color=cm.cool(i/float(taketime)),marker='.',markersize='1',linestyle='solid',linewidth='2')

if model == '0':
    plt.plot(t,y,color='gray',marker='.',markersize='0.1',linewidth='1')

plt.xlabel(r'$|k_{\rm{eff}}|/H_{\rm{initial}}$')
plt.ylabel(r'$P_{\delta\phi}$',rotation=0,fontsize=15)
plt.legend(legend_labels)

if model == '0':
    plt.xlim([1.0,200])
    plt.ylim([pow(10,-11),pow(10,-7)])
elif model == '1':
    plt.xlim([1.0,200])
    plt.ylim([pow(10,-12),pow(10,-6)])

plt.xscale('log', base = 10)
plt.yscale('log', base = 10)

plt.grid(True)
plt.show()

# Original Method
fig = plt.figure(figsize=(8,6))
plt.suptitle('Original method')

for i in range(taketime):
    plt.plot(k_ori[i],f_ori[i],color=cm.cool(i/float(taketime)),marker='.',markersize='1',linestyle='solid',linewidth='2')

if model == '0':
    plt.plot(t,y,color='gray',marker='.',markersize='0.1',linewidth='1')

plt.xlabel(r'$|k_{\rm{eff}}|/H_{\rm{initial}}$')
plt.ylabel(r'$P_{\delta\phi}$',rotation=0,fontsize=15)
plt.legend(legend_labels)

if model == '0':
    plt.xlim([1.0,200])
    plt.ylim([pow(10,-11),pow(10,-7)])
elif model == '1':
    plt.xlim([1.0,200])
    plt.ylim([pow(10,-12),pow(10,-6)])

plt.xscale('log', base = 10)
plt.yscale('log', base = 10)

plt.grid(True)
plt.show()