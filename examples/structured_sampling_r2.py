
from mylibs.structured_sample import StructuredSample
from myutils.plot_helper import PlotHelper

sam = StructuredSample()

str_plt = PlotHelper(['', ''], fancy=True, pgf=False)
import matplotlib.pyplot as plt

samples = sam.generate_sample_plan(30, 2, [(0, 30), (0, 30)])
for i in range(0, len(samples)):
    str_plt.ax.plot([samples[i][0]], [samples[i][1]], 'bo')

str_plt.ax.xaxis.set_ticklabels([])
str_plt.ax.yaxis.set_ticklabels([])
#str_plt.ax.set_xticks(range(0, 31), minor=False)
#str_plt.ax.set_yticks(range(0, 31), minor=False)
str_plt.ax.locator_params(nbins=4, axis='x')
str_plt.ax.locator_params(nbins=4, axis='y')
str_plt.ax.grid(True)

str_plt.finalize(width=1.5, height=1.5, show_legend=False)
str_plt.save('../data_out/plot/structuredSamp.pdf')
str_plt.show()