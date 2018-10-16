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

from wingConstruction.wingUtils.defines import *
from wingConstruction.surrogate_v2 import surrogate_analysis
from wingConstruction.surrogate_v2 import SurroResults
from wingConstruction.wingUtils.constants import Constants
from myutils.plot_helper import PlotHelper

def run_analysis():
    surro_methods = [SURRO_POLYNOM, SURRO_KRIGING]  # SURRO_KRIGING, SURRO_RBF, SURRO_POLYNOM
    sample_methods = [SAMPLE_LATIN, SAMPLE_HALTON, SAMPLE_STRUCTURE]  # SAMPLE_LATIN, SAMPLE_HALTON
    sample_point_count = list(range(2, 26))
    use_abaqus = False
    use_pgf = False
    job_count = len(surro_methods) * len(sample_methods) * sum(sample_point_count)
    jobs_done = 0
    print('required runs: {:d}'.format(job_count))

    output_file_name = 'surro_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
    output_f = open(Constants().WORKING_DIR + '/'
                    + output_file_name,
                    'w')
    output_f.write('SampleMethod,SampleMethodID,SamplePointCound,SurroMehtod,SurroMehtodID,deviation,rmse,mae,press,optRib,optShell,optWight,runtime,errorStr\n')

    for surro_m in surro_methods:
        for sample_m in sample_methods:
            for sample_points in sample_point_count:
                print('##################################################################################')
                print('next run: surro: {:s}, sample: {:s}, points: {:d}'.format(SURRO_NAMES[surro_m], SAMPLE_NAMES[sample_m], sample_points))
                try:
                    res, _ = surrogate_analysis(sample_m,
                                             sample_points,
                                             surro_m,
                                             use_abaqus=use_abaqus,
                                             pgf=use_pgf,
                                             show_plots=False,
                                             force_recalc=False)
                except Exception as e:
                    print('ERROR ' + str(e))
                    res = SurroResults()
                    res.errorStr = 'general fail: ' + str(e)
                output_f.write(SAMPLE_NAMES[sample_m] + ','
                               + '{:d}'.format(sample_m) + ','
                               + '{:d}'.format(sample_points) + ','
                               + SURRO_NAMES[surro_m] + ','
                               + '{:d}'.format(surro_m) + ','
                               + '{:f}'.format(res.valiResults.deviation) + ','
                               + '{:f}'.format(res.valiResults.rmse) + ','
                               + '{:f}'.format(res.valiResults.mae) + ','
                               + '{:f}'.format(res.valiResults.press) + ','
                               + '{:f}'.format(res.optimumRib) + ','
                               + '{:f}'.format(res.optimumShell) + ','
                               + '{:f}'.format(res.optimumWeights) + ','
                               + '{:f}'.format(res.runtime) + ','
                               + res.errorStr.replace(',', ';') + '\n')
                output_f.flush()
                jobs_done += sample_points
                print('jobs done: {:d}/{:d} -> {:f}%'.format(jobs_done, job_count, 100. * jobs_done / job_count))

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
    samp_plot = PlotHelper(['Anzahl der Sampling Punkte', 'Abweichung in %'], fancy=False, pgf=False)
    # plot one % line
    samp_plot.ax.plot([0, max(sampling_point_count)], [1., 1.], 'k--', label='1%-Linie')
    for key in sampling_data:
        x = [x for x,y in sampling_data[key]]
        y = [y for x,y in sampling_data[key]]
        y = np.array(y) * 100. # make it percent
        samp_plot.ax.plot(x, y, 'x-', label=key)
    samp_plot.ax.set_ylim([0, 8.])
    samp_plot.ax.set_xlim([0, max(sampling_point_count)])
    samp_plot.finalize()
    samp_plot.save(Constants().PLOT_PATH + 'samplePlanCompare.pdf')
    samp_plot.show()


if __name__ == '__main__':
    file = run_analysis()
    #plot_sample_point_analysis(file)
    plot_sample_point_analysis('surro_2018-09-15_12_59_09.csv')
