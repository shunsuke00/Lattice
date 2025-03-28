import matplotlib.pyplot as plt
import math

# Parameter input
#---------------------------------------------------------------
choice    = 0 # 0->thesis or 1->original

gamma = 0.2
thres = 0.4
gstar = 106.75
Hcmb  = 1.943   # this is refferenced value of Ncmb=55
kcmb  = 0.05
Ncmb  = 56.5     # default is 55, Planck range is [50,60]
Ngap  = 68.5-96.756+Ncmb
# we can't correct Hcmb, assume Hcmb is not affected


H_list    = []
aH_list   = []
N_list    = []

k_list = [[],[]]
P_list = [[],[]]


# File Loading (float in python = double in C++)
#---------------------------------------------------------------
path  = "/Users/tsuchitashunsuke/Dropbox/lattice/AxionInflation"
path += "/data_new/bumpy/256"
path += "/axion-U(1)_adjust"
path0 = path + "/python.dat"
path1 = path + "/sf_0.dat"
path2 = path + "/spectratimes_0.dat"
path3 = path + "/spectra_ax_0.dat"
path4 = path + "/spectra_eff_ax_0.dat"

# python.dat
'''かつての設定
num = 0
f0 = open(path0, "r")
for line in f0:
    data = line[:-1].split(' ')
    if(num==0):
        num_point = float(data[0])
        L_pr      = float(data[1])
        num += 1
dx_pr = L_pr/num_point
coeff = pow(pow(dx_pr,2)/L_pr,3)/(pow(math.pi,2)*2.)'''
num = 0
f0 = open(path0, "r")
for line in f0:
    data = line[:-1].split(' ')
    if(num==0):
        num_point = float(data[3])
        L_pr      = float(data[4])
        num += 1
dx_pr = L_pr/num_point
coeff = pow(pow(dx_pr,2)/L_pr,3)/(pow(math.pi,2)*2.)

# sf_0.dat
f1 = open(path1,"r")
for line in f1:
    data = line[:-1].split(' ')
    H_list.append(float(data[2]))
    aH_list.append(float(data[1])*float(data[2]))
    N_list.append(float(data[3]))

# spectratimes_0.dat
line_num = 0
f2  = open(path2,"r")
for line in f2:
    line_num += 1
    data = line[:-1].split(' ')
    a    = float(data[1])
    ad   = float(data[2])
    fd   = float(data[3])
    zeta_factor=pow(ad/(a*fd),2)
f2.close()

# spectra_ax_0.dat
count = 0
f3 = open(path3,"r")
for line in f3:
    if count == line_num-1:
        data = line[:-1].split(' ')
        if len(data) >= 5:
            k = float(data[3])
            if k!=0:
                k_list[0].append(k)
                P_list[0].append(float(data[4])*pow(k,3)*coeff*zeta_factor)
    elif not line.strip():
        count += 1

# spectra_eff_ax_0.dat
count = 0
f4 = open(path4,"r")
for line in f4:
    if count == line_num-1:
        data = line[:-1].split(' ')
        if len(data) >= 4:
            k = float(data[2])
            if k!=0 and float(data[0])!=0:
                k_list[1].append(k)
                P_list[1].append(float(data[3])*pow(k,3)*coeff*zeta_factor)
    elif not line.strip():
        count += 1


# Calculating M_PBH for each k_pr
#------------------------------------------------------------------
M_list  = []
dN_list = []
for j in range(len(k_list[choice])):
    index=0
    k = k_list[choice][j]
    if k<aH_list[0]:
        M_list.append(0.)
    else:
        while k>aH_list[index]:
            index += 1
            if index > 100000:
                print("ERROR: k is too large.\n")
                break
        H_be  = H_list[index-1]
        H_af  = H_list[index]
        aH_be = aH_list[index-1]
        aH_af = aH_list[index]
        N_be  = N_list[index-1]
        N_af  = N_list[index]

        deltaN = (k-aH_be)*(N_af-N_be)/(aH_af-aH_be)
        deltaH = deltaN*(H_af-H_be)/(N_af-N_be)
        if deltaN<=0. or deltaH>=0.:
            M_list.append(0.)
            dN_list.append(0.)
        elif abs((N_af-N_be)/N_af)>1. or abs((H_af-H_be)/H_af>1.):
            print("ERROR: The linear approximation is not accurate.\n")
        else:
            mass = 2.79*pow(10.,-12.)*(gamma/0.2)*pow(106.75/gstar,1./6.)*pow(pow(10.,12.)*Hcmb/(kcmb*(H_be+deltaH)),2.)*math.exp(-2.*(Ngap+N_be+deltaN))
            M_list.append(mass)
            dN_list.append(N_be+deltaN)


# Calculating f_PBH
#------------------------------------------------------------------
M_plotlist = []
f_plotlist = []
xmin = 3
xmax = len(k_list[choice])-5
for index in range(5,len(k_list[choice])-5,1):
    kref  = k_list[choice][index]
    newk_list = []
    for i in range(len(k_list[choice])):
        newk_list.append(k_list[choice][i]/kref)

    # calculating sigma (by Trapezoidal approximation)
    sigma = 0
    for i in range(xmin,xmax,1):
        LBottom = P_list[choice][i]*pow(newk_list[i],3)*math.exp(-pow(newk_list[i],2))
        RBottom = P_list[choice][i+1]*pow(newk_list[i+1],3)*math.exp(-pow(newk_list[i+1],2))
        sigma  += 0.5*(newk_list[i+1]-newk_list[i])*(LBottom+RBottom)
    sigma = math.sqrt(sigma*16./81.)

    beta = gamma*0.5*math.erfc(thres/(math.sqrt(2.)*sigma))
    abundance = beta*math.sqrt(gamma/0.2)/(math.sqrt(M_list[index])*5.91476*pow(10,-9))
    M_plotlist.append(M_list[index])
    f_plotlist.append(abundance)


# Graph Drawing
#------------------------------------------------------------------
if choice==0:
    methodname = "(Thesis method)"
elif choice==1:
    methodname = "(Original method)"

# M_PBH
fig = plt.figure(figsize=(9,4))
plt.suptitle(r'$k-M_\mathrm{PBH}\ and\ k-\Delta N_e\ $' + methodname)

plt.subplot(1,2,1)
plt.plot(k_list[choice], M_list, marker='.', markersize='2',linestyle='None')
plt.xlabel(r'$k_\mathrm{eff,pr}$')
plt.ylabel(r'$M_\mathrm{PBH}$')

plt.subplot(1,2,2)
plt.plot(k_list[choice], dN_list, marker='.', markersize='2',linestyle='None')
plt.xlabel(r'$k_\mathrm{eff,pr}$')
plt.ylabel(r'$\Delta N_e$')

plt.grid(True)
plt.show()

# Curvature Spectrum
fig = plt.figure(figsize=(9,7))
plt.suptitle("Curvature Power Spectrum " + methodname)
plt.plot(k_list[choice], P_list[choice], marker='.', markersize='1')
plt.xlabel(r'$k_\mathrm{eff,pr}$')
plt.ylabel(r'$P_\mathcal{R}$')
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.show()

# PBH abundance
fig = plt.figure(figsize=(9,7))
plt.suptitle("PBH abundance " + methodname + " Ncmb=" + str(Ncmb))
plt.plot(M_plotlist, f_plotlist, marker='.', markersize='1',linestyle='-')
plt.xlabel(r'$M_{PBH}$')
plt.ylabel(r'$f_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([0.00001,1])
plt.grid(True)
plt.show()

'''
with open('abundance.dat', 'w') as f:
    # list1とlist2を対応させてファイルに書き込む
    for item1, item2 in zip(M_plotlist, f_plotlist):
        f.write(f"{item1} {item2}\n")'''