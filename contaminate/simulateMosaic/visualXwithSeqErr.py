# plot simulation results for simXwithRC.py/sh
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.function_base import extract
from visualXwithRC import readResult

if __name__ == '__main__':
    for cov in [0.5, 0.1]:
        prefix = "/mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXseqErr/"
        confix = ""
        if cov == 0.5:
            confix += "chrX_cov1over2"
        elif cov == 0.1:
            confix += "chrX_cov1over10"

        prefix += confix

        conMLEs = []
        errBars = []
        xs = []
        left = -3 - 0.25*(2/9)
        right = -3 + 0.25*(2/9)
        for i in range(1, 11):
            conMLE, errBar = readResult(f'{prefix}/err{i}/batchresults.txt')
            if i == 1:
                plt.scatter(np.logspace(left, right, len(conMLE)), conMLE, marker='o', c='blue', s=5, zorder=3, label="MLE for contamination rate")
                plt.errorbar(np.logspace(left, right, len(conMLE)), conMLE, yerr=errBar, fmt='none', ecolor='#8c8c8c', zorder=1, label="95% CI")
            else:
                plt.scatter(np.logspace(left, right, len(conMLE)), conMLE, marker='o', c='blue', s=5, zorder=3)
                plt.errorbar(np.logspace(left, right, len(conMLE)), conMLE, yerr=errBar, fmt='none', ecolor='#8c8c8c', zorder=1)

            left += 2/9
            right += 2/9

        plt.title(f'1240k SNP coverage on chrX : {cov}X')
        plt.ylabel('estimated contamination')
        plt.xlabel('simulated genotyping error rate')
        plt.axhline(y=0.075, xmin=0, xmax=1, zorder=2, c='red', linestyle='-', label="true contamination rate")
        plt.axvline(x=1e-2, ymin=0, ymax=1, zorder=2, c='black', linestyle='dashdot', label='genotyping error rate used in inference')
        plt.legend(loc="upper left", fontsize='x-small')
        plt.xscale('log', basex=10)
        plt.savefig(f'{prefix}/{confix}.png', dpi=300)
        plt.clf()

        

