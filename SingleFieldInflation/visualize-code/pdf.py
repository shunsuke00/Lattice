import matplotlib.pyplot as plt
import matplotlib.cm as cm

#ファイル読み込み
#-------------------------------------------------------
#-------------------------------------------------------
path0 = "/path"
path0 += "/SFInflation"
path1 = path0 + "/pdf0_0.dat"
path2 = path0 + "/spectratimes_0.dat"

taketime = 11
pdflist  = [[] for _ in range(taketime)]
xlist    = [[] for _ in range(taketime)]

N_list = []


#読み込み
#-------------------------------------------------------
num = 0
f1 = open(path1, "r")
for line in f1:
    if line.strip() == "":
        num += 1
        if num == taketime:
            break
        else:
            continue
    data = line[:-1].split(' ')
    xlist[num].append(float(data[0]))
    pdflist[num].append(float(data[2]))


f2 = open(path2, "r")
for line in f2:
    data      = line[:-1].split(' ')
    N_list.append(round(float(data[1]),1))


legend_labels = [r"$N_e=$" + str(N_list[i]) for i in range(taketime)]

# ヒストグラム
#-------------------------------------------------------
#-------------------------------------------------------
fig = plt.figure(figsize=(8,6))
plt.suptitle('Inflaton fluctuation PDF (bin=0.1)')

for num in range(taketime):
    plt.plot(xlist[num],pdflist[num],color=cm.cool(num/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1')
plt.xlabel(r'$\delta\phi\ [\sigma]$')
plt.legend(legend_labels)

plt.grid(True)
plt.show()


fig = plt.figure(figsize=(8,6))
plt.suptitle('Inflaton fluctuation PDF (bin=0.1)')

for num in range(taketime):
    plt.plot(xlist[num],pdflist[num],color=cm.cool(num/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1')
plt.xlabel(r'$\delta\phi\ [\sigma]$')
plt.legend(legend_labels)
plt.yscale('log')

plt.grid(True)
plt.show()
