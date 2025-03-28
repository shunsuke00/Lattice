import matplotlib.pyplot as plt

#モデルを決めるパラメータ
#"single"(単一場)、"axion-U(1)"(axion-U(1))
#model = "axion-U(1)"
model = "single"

#各値を格納するリスト
t_list    = []
N_list    = []
a_list    = []
H_list    = []
f_list    = []
df_list   = []
ener_list = []
E_list    = []
epsi_list = []
eta_list  = []

xi_list   = []
BR_list   = []

list_list = [t_list, N_list, a_list, H_list, f_list, df_list, 
             ener_list, E_list, epsi_list, eta_list]

#ファイル読み込み
#---------------------------------------------------------------
path = "/Users/tsuchitashunsuke/Dropbox/lattice/background/status.dat"
f = open(path, "r")

for line in f:
    #数値以外は捨てる
    if line[0] == "#":
        continue
    
    #数値をリストに追加
    data = line[:-1].split('  ')

    for i in range(10):
        list_list[i].append(float(data[i]))
    if model=="axion-U(1)":
        xi_list.append(float(data[10]))
        BR_list.append(abs(float(data[11])))


#グラフの描画
#------------------------------------------------------------------
#Figure 3.2再現
fig = plt.figure(figsize=(10,7))
fig.suptitle('background')

plt.subplot(2,2,1)
plt.plot(N_list,f_list, '.',markersize='1', linestyle='None')
plt.xlabel(r'$N_e$')
plt.ylabel(r'$\phi [M_{pl}]$')

plt.subplot(2,2,2)
plt.plot(N_list,df_list, '.',markersize='1', linestyle='None')
plt.xlabel(r'$N_e$')
plt.ylabel(r'$df/m [M_{pl}]$')
#plt.ylim([0.8144,0.8152])

plt.subplot(2,2,3)
plt.plot(N_list,H_list, '.',markersize='1', linestyle='None')
plt.xlabel(r'$N_e$')
plt.ylabel(r'$H/m$')

plt.subplot(2,2,4)
plt.plot(N_list,E_list, '.',markersize='1', linestyle='None')
plt.xlabel(r'$N_e$')
plt.ylabel(r'$E$')
plt.ylim([0.999,1.001])

plt.grid(True)
plt.show()

#スローロールパラメータ
fig = plt.figure(figsize=(10,5))
fig.suptitle('Slow-Roll parameters')

plt.subplot(1,2,1)
plt.plot(N_list,epsi_list, '.',markersize='1', linestyle='None')
plt.xlabel(r'$N_e$')
plt.ylabel(r'$\epsilon$')
plt.yscale('log')

plt.subplot(1,2,2)
plt.plot(N_list,eta_list, '.',markersize='1', linestyle='None')
plt.xlabel(r'$N_e$')
plt.ylabel(r'$\eta$')
plt.ylim([-6,3])

plt.grid(True)
plt.show()

#xi
if model=="axion-U(1)":
    fig = plt.figure(figsize=(10,8))

    plt.subplot(2,1,1)
    plt.plot(N_list,xi_list,'.',markersize='1',linestyle='None')
    plt.xlabel('N')
    plt.ylabel('xi')

    plt.subplot(2,1,2)
    plt.plot(N_list,BR_list,'.',markersize='1',linestyle='None')
    plt.xlabel('N')
    plt.ylabel('condition')
    plt.yscale('log')

    plt.grid(True)
    plt.show()


