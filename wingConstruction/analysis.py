__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :analysis of surrogate with different methods
# author          :Juri Bieler
# date            :2018-09-10
# notes           :
# python_version  :3.6
# ==============================================================================

from datetime import datetime
import numpy as np

from wingConstruction.utils.defines import *
from wingConstruction.surrogateV2 import surrogate_analysis
from wingConstruction.surrogateV2 import SurroResults
from wingConstruction.utils.Constants import Constants
from utils.PlotHelper import PlotHelper

def run_analysis():
    surroMethods = [SURRO_KRIGING]  # SURRO_KRIGING, SURRO_RBF
    sampleMethods = [SAMPLE_LATIN, SAMPLE_HALTON, SAMPLE_STRUCTURE]  # SAMPLE_LATIN, SAMPLE_HALTON
    samplePointCount = list(range(2, 41))
    useAbaqus = False
    usePGF = False
    jobCount = len(surroMethods) * len(sampleMethods) * sum(samplePointCount)
    jobsDone = 0
    print('required runs: {:d}'.format(jobCount))

    output_file_name = 'surro_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
    output_f = open(Constants().WORKING_DIR + '/'
                    + output_file_name,
                    'w')
    output_f.write('SampleMethod,SampleMethodID,SamplePointCound,SurroMehtod,SurroMehtodID,deviation,optRib,optShell,optWight,runtime,errorStr\n')

    for surroM in surroMethods:
        for sampleM in sampleMethods:
            for samplePoints in samplePointCount:
                try:
                    res = surrogate_analysis(sampleM,
                                             samplePoints,
                                             surroM,
                                             use_abaqus=useAbaqus,
                                             pgf=usePGF,
                                             show_plots=False,
                                             force_recalc=False)
                except Exception as e:
                    print('ERROR ' + str(e))
                    res = SurroResults()
                    res.errorStr = 'general fail: ' + str(e)
                output_f.write(SAMPLE_NAMES[sampleM] + ','
                               + '{:d}'.format(sampleM) + ','
                               + '{:d}'.format(samplePoints) + ','
                               + SURRO_NAMES[surroM] + ','
                               + '{:d}'.format(surroM) + ','
                               + '{:f}'.format(res.deviation) + ','
                               + '{:f}'.format(res.optimumRib) + ','
                               + '{:f}'.format(res.optimumShell) + ','
                               + '{:f}'.format(res.optimumWeights) + ','
                               + '{:f}'.format(res.runtime) + ','
                               + res.errorStr + '\n')
                output_f.flush()
                jobsDone += samplePoints
                print('#########################################')
                print('### jobs done: {:d}/{:d} -> {:f}%'.format(jobsDone, jobCount, 100. * jobsDone / jobCount))

    output_f.close()
    return output_file_name

def plot_sample_point_analysis(file_name):
    file_path = Constants().WORKING_DIR + '/' + file_name
    data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    sampling_plan_id = data[:, 1]
    sampling_point_count = data[:, 2]
    deviation = data[:, 5]
    sampling_data = {}
    for samp in SAMPLE_NAMES:
        sampling_data[samp] = []
    for i in range(0, len(sampling_plan_id)):
        sampling_data[SAMPLE_NAMES[int(sampling_plan_id[i])]].append((sampling_point_count[i], deviation[i]))
    sampPlot = PlotHelper(['Anzahl der Sampling Punkte', 'Abweichung in %'], fancy=False)
    for key in sampling_data:
        x = [x for x,y in sampling_data[key]]
        y = [y for x,y in sampling_data[key]]
        y = np.array(y) * 100. # make it percent
        sampPlot.ax.plot(x, y, 'x-', label=key)
    sampPlot.finalize()
    sampPlot.show()


if __name__ == '__main__':
    file = run_analysis()
    plot_sample_point_analysis(file)
    #plot_sample_point_analysis('surro_2018-09-13_13_45_10.csv')
