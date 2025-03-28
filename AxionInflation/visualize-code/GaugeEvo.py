import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

analytical_plot = False

# File Loading (float in python = double in C++)
#---------------------------------------------------------------
path = "/Users/tsuchitashunsuke/Dropbox/lattice/AxionInflation"
path += "/axion-U(1)"
path += "_adjust"
#path += "_coulomb"
path0 = path + "/python.dat"
path1 = path + "/spectratimes_0.dat"
path2 = path + "/spectra_gf_0.dat"
path3 = path + "/spectra_eff_gf_0.dat"

if analytical_plot:
    path_analy = "/Users/tsuchitashunsuke/Dropbox/lattice/coulomb/mathematica"
    path4 = path_analy + "/gaugespectra_weak_negative.dat"
    path5 = path_analy + "/gaugespectra_weak_positive.dat"

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

dx_pr = L_pr/num_point
coeff = pow(pow(dx_pr,2)/(m*L_pr),3)    # Rescaling

# spectratimes_0.dat
N_list = []
f1 = open(path1,"r")
for line in f1:
    data = line[:-1].split(' ')
    N_list.append(round(float(data[5]),1))


# Some Setting
taketime = 9
k_lat     = [[] for _ in range(taketime)]
fposi_lat = [[] for _ in range(taketime)]
fnega_lat = [[] for _ in range(taketime)]
k_the     = [[] for _ in range(taketime)]
fposi_the = [[] for _ in range(taketime)]
fnega_the = [[] for _ in range(taketime)]
k_ori     = [[] for _ in range(taketime)]
fposi_ori = [[] for _ in range(taketime)]
fnega_ori = [[] for _ in range(taketime)]
index_list = [[] for _ in range(taketime)]

length = math.floor(math.sqrt(3)*num_point/2.)+2
for i in range(taketime):
    index_list[i] = i

# spectra_gf_0.dat
f2 = open(path2, "r")
num = 0
for line in f2:
    num += 1
    for index in range(taketime):
        if num>index_list[index]*length and num<(index_list[index]+1)*length:
            data = line[:-1].split(' ')
            k_pr = float(data[2])
            keff_pr = float(data[3])
            if(k_pr!=0):
                k_lat[index].append(k_pr)
                fposi_lat[index].append(coeff*float(data[4]))
                fnega_lat[index].append(coeff*float(data[5]))
            if(keff_pr!=0):
                k_the[index].append(keff_pr)
                fposi_the[index].append(coeff*float(data[4]))
                fnega_the[index].append(coeff*float(data[5]))

# spectra_eff_gf_0.dat
f3 = open(path3, "r")
num   = 0
for line in f3:
    num += 1
    for index in range(taketime):
        if num>index_list[index]*length and num<(index_list[index]+1)*length:
            data = line[:-1].split(' ')
            keff_pr = float(data[2])
            if(keff_pr!=0):
                k_ori[index].append(keff_pr)
                fposi_ori[index].append(coeff*float(data[3]))
                fnega_ori[index].append(coeff*float(data[4]))

# Analytical Solution
if analytical_plot:
    k = np.arange(10, 251, 2)
    asposi_list = [[] for _ in range(taketime)]
    asnega_list = [[] for _ in range(taketime)]

    f4 = open(path4, "r")
    num = 0
    for line in f4:
        index = num // 121
        asposi_list[index].append(pow(10,float(line)))
        num += 1

    f5 = open(path5, "r")
    num = 0
    for line in f5:
        index = num // 121
        asnega_list[index].append(pow(10,float(line)))
        num += 1

    asdami_list = []
    for i in range(len(k)):
        asdami_list.append(0.)



# Graph Drawing
#------------------------------------------------------------------
legend_labels = [r"$N_e=$" + str(N_list[index_list[i]]) for i in range(taketime)]
if analytical_plot:
    legend_labels.append(r"$A=\frac{1}{\sqrt{2k}}[G_0+iF_0]$")

    def analyticPlot_posi():
        for i in range(taketime):
            plt.plot(k,asposi_list[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.1',linestyle='dashed',linewidth='1',alpha=0.3)
        
    def analyticPlot_nega():
        for i in range(taketime):
            plt.plot(k,asnega_list[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.1',linestyle='dashed',linewidth='1',alpha=0.3)
    
# Latticeeasy Method
fig = plt.figure(figsize=(9,7))
fig.suptitle('Gauge Polarization Spectram (Latticeeasy Method)')

plt.subplot(2,1,1)
for i in range(taketime):
    plt.plot(k_lat[i],fposi_lat[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1.5')

if analytical_plot:
    plt.plot(k,asdami_list,color='black',linestyle='dashed')
    analyticPlot_posi()

plt.legend(legend_labels, bbox_to_anchor=(1, 1))
plt.ylabel(r'$|A_+|^2$', rotation=0,fontsize=12,labelpad=15,y=0.5)
plt.xlim([0,250])
plt.ylim([10,pow(10,11)])
plt.yscale('log')
plt.grid(which='major', color='gray', linestyle='--', linewidth=0.5)

plt.subplot(2,1,2)
for i in range(taketime):
    plt.plot(k_lat[i],fnega_lat[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1.5')

if analytical_plot:
    analyticPlot_nega()

plt.xlabel(r'$k_{eff}/m$')
plt.ylabel(r'$|A_-|^2$', rotation=0,fontsize=12,labelpad=15,y=0.5)
plt.xlim([0,250])
plt.ylim([10,pow(10,11)])
plt.yscale('log')
plt.grid(which='major', color='gray', linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.grid(True)
plt.show()

# Angelo Method
fig = plt.figure(figsize=(9,7))
fig.suptitle('Gauge Polarization Spectram (Angelo Method)')

plt.subplot(2,1,1)
for i in range(taketime):
    plt.plot(k_the[i],fposi_the[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1.5')

if analytical_plot:
    plt.plot(k,asdami_list,color='black',linestyle='dashed')
    analyticPlot_posi()

plt.legend(legend_labels, bbox_to_anchor=(1, 1))
plt.ylabel(r'$|A_+|^2$', rotation=0,fontsize=12,labelpad=15,y=0.5)
plt.xlim([0,250])
plt.ylim([10,pow(10,11)])
plt.yscale('log')
plt.grid(which='major', color='gray', linestyle='--', linewidth=0.5)

plt.subplot(2,1,2)
for i in range(taketime):
    plt.plot(k_the[i],fnega_the[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1.5')

if analytical_plot:
    analyticPlot_nega()

plt.xlabel(r'$k_{eff}/m$')
plt.ylabel(r'$|A_-|^2$', rotation=0,fontsize=12,labelpad=15,y=0.5)
plt.xlim([0,250])
plt.ylim([10,pow(10,11)])
plt.yscale('log')
plt.grid(which='major', color='gray', linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.grid(True)
plt.show()

# Original Method
fig = plt.figure(figsize=(9,7))
fig.suptitle('Gauge Polarization Spectram (Original Method)')

plt.subplot(2,1,1)
for i in range(taketime):
    plt.plot(k_ori[i],fposi_ori[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1.5')

if analytical_plot:
    plt.plot(k,asdami_list,color='black',linestyle='dashed')
    analyticPlot_posi()

plt.legend(legend_labels, bbox_to_anchor=(1, 1))
plt.ylabel(r'$|A_+|^2$', rotation=0,fontsize=12,labelpad=15,y=0.5)
plt.xlim([0,250])
plt.ylim([10,pow(10,11)])
plt.yscale('log')
plt.grid(which='major', color='gray', linestyle='--', linewidth=0.5)

plt.subplot(2,1,2)
for i in range(taketime):
    plt.plot(k_ori[i],fnega_ori[i],color=cm.cool(i/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1.5')

if analytical_plot:
    analyticPlot_nega()

plt.xlabel(r'$k_{eff}/m$')
plt.ylabel(r'$|A_-|^2$', rotation=0,fontsize=12,labelpad=15,y=0.5)
plt.xlim([0,250])
plt.ylim([10,pow(10,11)])
plt.yscale('log')
plt.grid(which='major', color='gray', linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.grid(True)
plt.show()