from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)

width = 0.30

def main():

    with open(argv[1], "r") as f:

        names = []
        time1 = []
        time2 = []
        size = None

        for line in f:
            a = np.array(line.split()).reshape(-1, 2)
            time1.append(float(a[2, 1])/10**6)
            time2.append(float(a[3, 1])/10**6)
            names.append(int(a[1, 1]))
            size = a[0, 1]

        names = names[::-1]
        time1 = time1[::-1]
        time2 = time2[::-1]

    tot_time = [i+j for i, j in zip(time1,time2)]
    x_bars = np.arange(len(names))

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15,5), sharex=True)

    ax[0].bar(x_bars, tot_time, align = 'center', color="#9966ff", width=width, label = "Time")

    header = "Time and number of processors, matrix size: "+size

    ax[0].set_xticks(x_bars)
    ax[0].set_xticklabels(names)
    ax[0].set_xlabel("Number of Processors")
    ax[0].set_yscale("log", nonposy='clip')
    ax[0].set_ylabel("Time, seconds")
    ax[0].set_title(header)
    ax[0].legend()


    speedup = [tot_time[0]/i for i in tot_time]
    eff = [s/p for p, s in zip(2**x_bars, speedup)]

    ax[1].bar(x_bars, speedup, align = 'center', color="#9966ff", width=width, label = "Speedup")
    ax[1].bar(x_bars+width, eff, align = 'center', color="#66ff99", width=width, label = "Efficiency")

    header = "Speedup and efficiency, matrix size: "+size

    ax[1].set_xticks(x_bars, names)
    ax[1].set_xticklabels(names)
    ax[1].set_xlabel("Number of Processors")
    ax[1].set_yscale("log", nonposy='clip')
    #ax[1].set_ylim(0,500)
    ax[1].set_ylabel("Speedup/Efficiency coefficient")
    ax[1].set_title(header)
    ax[1].legend()

    plt.tight_layout()
    plt.savefig("graph_"+size+".png")
    plt.show()

if __name__=="__main__":
    main()