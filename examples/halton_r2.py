
from mylibs.halton import Halton

hal = Halton()
print(hal.halton(10, 2))

#for i in range(0, 10):
#    point = hal.halton(i, 2)

#sys.exit(15)
import matplotlib.pyplot as plt
samples = hal.generate_sample_plan(14, 2, [(5, 20), (0.01, 0.05)], base=[2,3])
for i in range(0, len(samples)):
    plt.plot([samples[i][0]], [samples[i][1]], 'bo')
plt.show()

from myutils.plot_helper import PlotHelper
pltHalton = PlotHelper([], fancy=True, pgf=False)
import matplotlib.pyplot as plt



ax1 = pltHalton.fig.add_subplot(121)
ax2 = pltHalton.fig.add_subplot(122)

#ax = fig.add_subplot(1, 2, 1)
for i in range(0, 100):
    point = [hal.halton(i, 11), hal.halton(i, 29)]
    ax1.plot([point[1]], [point[0]], 'bo', markersize=3)
ax1.xaxis.set_ticklabels([])
ax1.yaxis.set_ticklabels([])

#ax = fig.add_subplot(1, 2, 2)
for i in range(0, 100):
    point = [hal.halton(i, 2), hal.halton(i, 19)]
    ax2.plot([point[1]], [point[0]], 'bo', markersize=3)
ax2.xaxis.set_ticklabels([])
ax2.yaxis.set_ticklabels([])

pltHalton.fig.set_size_inches(5, 2.5)
plt.tight_layout()
pltHalton.save('../data_out/plot/halton.pdf')
pltHalton.show()

#samples = hal.generate_sample_plan(14, 2, [(5, 20), (0.01, 0.05)])
#for i in range(0, 14):
#    plt.plot([samples[i][0]], [samples[i][1]], 'bo')
    #plt.show()