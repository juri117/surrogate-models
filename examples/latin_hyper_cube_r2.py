
import numpy as np
from mylibs.latin_hyper_cube import LatinHyperCube

sam = LatinHyperCube()

import matplotlib.pyplot as plt
samples = sam.generate_sample_plan(14, 2, [(5, 20), (0.01, 0.05), (1, 5)])
#samples = sam.enhanced_latin_hypercube_2_pow_x(16)
for i in range(0, len(samples)):
    plt.plot([samples[i][0]], [samples[i][1]], 'bo')
plt.show()

from myutils.plot_helper import PlotHelper
plt_latin = PlotHelper([], fancy=True, pgf=False)
import matplotlib.pyplot as plt

ax1 = plt_latin.fig.add_subplot(121)
ax2 = plt_latin.fig.add_subplot(122)

sample_mat_full = sam.enhanced_latin_hypercube(2, 16)
xy_full = sam.bool_mat_to_list(sample_mat_full)
plotXY_full = np.array(xy_full).T.tolist()
ax1.plot(plotXY_full[0], plotXY_full[1], 'bo', markersize=5)

deleteX = [0, 15, 4, 1]
deleteY = [0, 15, 1, 4]
for i in range(0, len(deleteX)):
    #ax1.plot([deleteX[i]], [deleteY[i]], 'rx', mew=2, ms=10)
    ax1.plot([-1, 16], [deleteY[i], deleteY[i]], 'r-', linewidth=2)
    ax1.plot([deleteX[i], deleteX[i]], [-1, 16], 'r-', linewidth=2)
    ax1.text(deleteX[i] - 0.3, deleteY[i] - 0.3, str(i+1), plt_latin.font, fontweight='bold')

ax1.set_xticks(range(0, 16), minor=False)
ax1.set_yticks(range(0, 16), minor=False)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])
ax1.set_xlim(16, -1)
ax1.set_ylim(16, -1)
ax1.grid(True)

#ax = fig.add_subplot(1, 2, 2)

sample_mat = sam.enhanced_latin_hypercube(2, 12)
xy = sam.bool_mat_to_list(sample_mat)
plotXY = np.array(xy).T.tolist()

ax2.plot(plotXY[0], plotXY[1], 'bo', markersize=5)
ax2.set_xticks(range(0,12), minor=False)
ax2.set_yticks(range(0,12), minor=False)
ax2.xaxis.set_ticklabels([])
ax2.yaxis.set_ticklabels([])
ax2.set_xlim(12, -1)
ax2.set_ylim(12, -1)
ax2.grid(True)
ax2.grid(True)
#plt.show()
plt_latin.fig.set_size_inches(5, 2.5)
plt.tight_layout()
plt_latin.save('../data_out/plot/latinHyper.pdf')
plt_latin.show()