# plot simulation results for simXwithRC.py/sh
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.function_base import extract
from visualXwithRC import readResult

if __name__ == '__main__':
    for cov in [0.1, 0.5, 1.0]:
        prefix = "/mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXjump/"
        confix = ""
        if cov == 0.5:
            confix += "chrX_cov1over2"
        elif cov == 0.1:
            confix += "chrX_cov1over10"
        elif cov == 1.0:
            confix += "chrX_cov1"

        prefix += confix

        conMLEs = []
        errBars = []
        xs = []
        left = 2 - 0.25*(1/9)
        right = 2 + 0.25*(1/9)
        for i in range(1, 11):
            conMLE, errBar = readResult(f'{prefix}/jump{i}/batchresults.txt')
            if i == 1:
                plt.scatter(np.logspace(left, right, len(conMLE)), conMLE, marker='o', c='blue', s=5, zorder=3, label="MLE for contamination rate")
                plt.errorbar(np.logspace(left, right, len(conMLE)), conMLE, yerr=errBar, fmt='none', ecolor='#8c8c8c', zorder=1, label="95% CI")
            else:
                plt.scatter(np.logspace(left, right, len(conMLE)), conMLE, marker='o', c='blue', s=5, zorder=3)
                plt.errorbar(np.logspace(left, right, len(conMLE)), conMLE, yerr=errBar, fmt='none', ecolor='#8c8c8c', zorder=1)

            left += 1/9
            right += 1/9

        plt.title(f'1240k SNP coverage: {cov}X')
        plt.ylabel('estimated contamination')
        #plt.ylim((0.0, 0.18)) # use the same y range to demonstrate the effects of coverage on the estimates
        plt.xlabel('simulated copying jump rate')
        plt.axhline(y=0.075, xmin=0, xmax=1, zorder=2, c='red', linestyle='-', label="true contamination rate")
        plt.axvline(x=300, ymin=0, ymax=1, zorder=2, c='black', linestyle='dashdot', label='jump rate used in inference')
        plt.legend(loc="upper left", fontsize='small')
        plt.xscale('log', basex=10)
        plt.savefig(f'{prefix}/{confix}.png', dpi=300)
        plt.clf()

        

