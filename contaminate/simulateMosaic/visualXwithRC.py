# plot simulation results for simXwithRC.py/sh
import numpy as np
import matplotlib.pyplot as plt

def readResult(filePath):
    # return two numpy array conMLE and conCI
    # conMLE: sorted array for MLE of contamination rate
    # conCI: sorted according to the order of conMLE, it's half the width of 95% CI, that is, 1.92*se
    conMLEs = []
    errBars = []
    with open(filePath) as file:
        for line in file:
            if line.startswith("#"):
                continue
            else:
                _, conMLE, lower, _ = line.strip().split("\t")
                conMLEs.append(float(conMLE))
                errBars.append(float(conMLE) - float(lower))
    conMLEs = np.array(conMLEs)
    errBars = np.array(errBars)
    sortedIndex = np.argsort(conMLEs)
    return conMLEs[sortedIndex], errBars[sortedIndex]

if __name__ == '__main__':
    for contamination in [0.0, 0.05, 0.1]:
        prefix = "/mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleX11/"
        confix = ""
        if contamination == 0.0:
            confix += "con0"
        elif contamination == 0.05:
            confix += "con5"
        elif contamination == 0.1:
            confix += "con10"

        prefix += confix
        conMLE_cov1over20, errBars_cov1over20 = readResult(f'{prefix}/chrX_cov1over20/batchresults_bfgs.txt')
        conMLE_cov1over10, errBars_cov1over10 = readResult(f'{prefix}/chrX_cov1over10/batchresults_bfgs.txt')
        conMLE_cov1over2, errBars_cov1over2 = readResult(f'{prefix}/chrX_cov1over2/batchresults_bfgs.txt')
        conMLE_cov1, errBars_cov1 = readResult(f'{prefix}/chrX_cov1/batchresults_bfgs.txt')
        conMLE_cov2, errBars_cov2 = readResult(f'{prefix}/chrX_cov2/batchresults_bfgs.txt')
        conMLE_cov5, errBars_cov5 = readResult(f'{prefix}/chrX_cov5/batchresults_bfgs.txt')

        xs_cov1over20 = np.linspace(2.5, 7.5, num=len(conMLE_cov1over20))
        xs_cov1over10 = np.linspace(12.5, 17.5, num=len(conMLE_cov1over10))
        xs_cov1over2 = np.linspace(22.5, 27.5, num=len(conMLE_cov1over2))
        xs_cov1 = np.linspace(32.5, 37.5, num=len(conMLE_cov1))
        xs_cov2 = np.linspace(42.5, 47.5, num=len(conMLE_cov2))
        xs_cov5 = np.linspace(52.5, 57.5, num=len(conMLE_cov5))

        plt.scatter(xs_cov1over20, conMLE_cov1over20, marker='o', c='blue', s=5, zorder=3, label="MLE for contamination rate")
        plt.errorbar(xs_cov1over20, conMLE_cov1over20, yerr=errBars_cov1over20, fmt='none', ecolor='#8c8c8c', zorder=1, label="95% CI")

        plt.scatter(xs_cov1over10, conMLE_cov1over10, marker='o', c='blue', s=5, zorder=3)
        plt.errorbar(xs_cov1over10, conMLE_cov1over10, yerr=errBars_cov1over10, fmt='none', ecolor='#8c8c8c', zorder=1)

        plt.scatter(xs_cov1over2, conMLE_cov1over2, marker='o', c='blue', s=5, zorder=3)
        plt.errorbar(xs_cov1over2, conMLE_cov1over2, yerr=errBars_cov1over2, fmt='none', ecolor='#8c8c8c', zorder=1)

        plt.scatter(xs_cov1, conMLE_cov1, marker='o', c='blue', s=5, zorder=3)
        plt.errorbar(xs_cov1, conMLE_cov1, yerr=errBars_cov1, fmt='none', ecolor='#8c8c8c', zorder=1)

        plt.scatter(xs_cov2, conMLE_cov2, marker='o', c='blue', s=5, zorder=3)
        plt.errorbar(xs_cov2, conMLE_cov2, yerr=errBars_cov2, fmt='none', ecolor='#8c8c8c', zorder=1)

        plt.scatter(xs_cov5, conMLE_cov5, marker='o', c='blue', s=5, zorder=3)
        plt.errorbar(xs_cov5, conMLE_cov5, yerr=errBars_cov5, fmt='none', ecolor='#8c8c8c', zorder=1)

        plt.title(f'contamination rate: {contamination}')
        plt.ylabel('estimated contamination')
        plt.xlabel('coverage')
        plt.axhline(y=contamination, xmin=0, xmax=1, zorder=2, c='red', linestyle='-', label="true contamination rate")
        plt.legend(loc="upper right")
        plt.xticks([5, 15, 25, 35, 45, 55], ["0.05", "0.1", "0.5", "1", "2", "5"])
        plt.savefig(f'{prefix}/{confix}_bfgs.png', dpi=300)
        plt.clf()







        

