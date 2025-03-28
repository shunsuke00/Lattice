#This file output the graph from mathematica CFW value file

import matplotlib.pyplot as plt
import math


#各値を格納するリスト
x_list = []
growAmp_list = []
growFac_list = []
stabAmp_list = []
stabFac_list = []

#ファイル読み込み
#pythonのfloatはC++のdoubleに対応する
#---------------------------------------------------------------
path0 = "/Users/tsuchitashunsuke/Dropbox/lattice/coulomb/CWFvalue/bumpy/64/H=1.92/L=3.5"
path1 = path0 + "/x0_list.dat"
path2 = path0 + "/growAmp.dat"
path3 = path0 + "/growFac.dat"
path4 = path0 + "/stabAmp.dat"
path5 = path0 + "/stabFac.dat"

# kのリスト
f1 = open(path1, "r")
for line in f1:
    x_list.append(float(line))

# CWFのリスト
f2 = open(path2, "r")
for line in f2:
    growAmp_list.append(float(line))

f3 = open(path3, "r")
for line in f3:
    growFac_list.append(float(line))

f4 = open(path4, "r")
for line in f4:
    stabAmp_list.append(float(line))

f5 = open(path5, "r")
for line in f5:
    stabFac_list.append(float(line))


#グラフの描画
#------------------------------------------------------------------
fig = plt.figure(figsize=(10,6))
plt.suptitle('growAmp')
plt.plot(x_list,growAmp_list,marker='.',markersize='1',linestyle='None',label='F_0(xi,x)')
plt.xlabel('x')
plt.ylabel('data')
plt.ylim([1.,1.5])
plt.grid(True)
plt.show()

fig = plt.figure(figsize=(10,6))
plt.suptitle('growFac')
plt.plot(x_list,growFac_list,marker='.',markersize='1',linestyle='None',label='F_0(xi,x)')
plt.xlabel('x')
plt.ylabel('data')
plt.ylim([0.9,1.])
plt.grid(True)
plt.show()

fig = plt.figure(figsize=(10,6))
plt.suptitle('stabAmp')
plt.plot(x_list,stabAmp_list,marker='.',markersize='1',linestyle='None',label='F_0(xi,x)')
plt.xlabel('x')
plt.ylabel('data')
plt.ylim([0.85,1.0])
plt.grid(True)
plt.show()

fig = plt.figure(figsize=(10,6))
plt.suptitle('stabFac')
plt.plot(x_list,stabFac_list,marker='.',markersize='1',linestyle='None',label='F_0(xi,x)')
plt.xlabel('x')
plt.ylabel('data')
plt.ylim([-0.9,-0.8])
plt.grid(True)
plt.show()
