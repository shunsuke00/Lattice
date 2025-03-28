import matplotlib.pyplot as plt
import matplotlib.cm as cm

# File Loading (float in python = double in C++)
#-------------------------------------------------------
path0 = "/Users/tsuchitashunsuke/Dropbox/lattice/AxionInflation"
path0 += "/axion-U(1)"
path0 += "_adjust"
path0 += "_coulomb"
path1 = path0 + "/pdf_ax_0.dat"
path2 = path0 + "/spectratimes_0.dat"

taketime = 11
pdflist  = [[] for _ in range(taketime)]
xlist    = [[] for _ in range(taketime)]

N_list = []

# pdf_ax_0.dat
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

# spectratimes_0.dat
f2 = open(path2, "r")
for line in f2:
    data      = line[:-1].split(' ')
    N_list.append(round(float(data[5]),1))


# Graph Drawing
#-------------------------------------------------------
legend_labels = [r"$N_e=$" + str(N_list[i]) for i in range(taketime)]

# Inflaton PDF
fig = plt.figure(figsize=(8,6))
plt.suptitle('Inflaton fluctuation PDF (bin=0.1)')

for num in range(taketime):
    plt.plot(xlist[num],pdflist[num],color=cm.cool(num/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1')
plt.xlabel(r'$\delta\phi\ [\sigma]$')
plt.legend(legend_labels)

plt.grid(True)
plt.show()

# Inflaton PDF (Log scale)
fig = plt.figure(figsize=(8,6))
plt.suptitle('Inflaton fluctuation PDF (bin=0.1)')

for num in range(taketime):
    plt.plot(xlist[num],pdflist[num],color=cm.cool(num/float(taketime)),marker='.',markersize='0.1',linestyle='solid',linewidth='1')
plt.xlabel(r'$\delta\phi\ [\sigma]$')
plt.legend(legend_labels)
plt.yscale('log')

plt.grid(True)
plt.show()
