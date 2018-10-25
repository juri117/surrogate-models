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

from wingconstruction.wingutils.defines import *
from wingconstruction.surrogate_v3 import Surrogate
from wingconstruction.surrogate_v3 import SurroResults
from wingconstruction.wingutils.constants import Constants
from myutils.plot_helper import PlotHelper

def run_analysis():
    surro_methods = [SURRO_POLYNOM]  # SURRO_KRIGING, SURRO_RBF, SURRO_POLYNOM, SURRO_PYKRIGING, SURRO_RBF_SCIPY
    sample_methods = [SAMPLE_STRUCTURE, SAMPLE_LATIN, SAMPLE_HALTON]  # SAMPLE_LATIN, SAMPLE_HALTON
    sample_point_count = list(range(3, 30+1))
    use_abaqus = True
    use_pgf = False
    job_count = len(surro_methods) * len(sample_methods) * sum(sample_point_count)
    jobs_done = 0
    print('required runs: {:d}'.format(job_count))

    output_file_name = 'analysis_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
    output_f = open(Constants().WORKING_DIR + '/'
                    + output_file_name,
                    'w')
    output_f.write('SampleMethod,SampleMethodID,SamplePointCound,SurroMehtod,SurroMehtodID,deviation,rmse,mae,press,optRib,optShell,optWight,optStress,optiParam1,optiParam2,runtime,errorStr\n')

    for surro_m in surro_methods:
        for sample_m in sample_methods:
            for sample_points in sample_point_count:
                print('##################################################################################')
                print('next run: surro: {:s}, sample: {:s}, points: {:d}'.format(SURRO_NAMES[surro_m], SAMPLE_NAMES[sample_m], sample_points))
                try:
                    sur = Surrogate(use_abaqus=True, pgf=use_pgf, show_plots=False, scale_it=True)
                    res, _ = sur.auto_run(sample_m, sample_points, surro_m, run_validation=True)
                except Exception as e:
                    print('ERROR ' + str(e))
                    res = SurroResults()
                    res.errorStr = 'general fail: ' + str(e)

                opti_param1 = -1
                opti_param2 = -1
                if len(res.opti_params) > 0:
                    opti_param1 = res.opti_params[0]
                if len(res.opti_params) > 1:
                    opti_param2 = res.opti_params[1]
                output_f.write(SAMPLE_NAMES[sample_m] + ','
                               + '{:d}'.format(sample_m) + ','
                               + '{:d}'.format(sample_points) + ','
                               + SURRO_NAMES[surro_m] + ','
                               + '{:d}'.format(surro_m) + ','
                               + '{:f}'.format(res.vali_results.deviation) + ','
                               + '{:f}'.format(res.vali_results.rmse) + ','
                               + '{:f}'.format(res.vali_results.mae) + ','
                               + '{:f}'.format(res.vali_results.press) + ','
                               + '{:f}'.format(res.optimum_rib) + ','
                               + '{:f}'.format(res.optimum_shell) + ','
                               + '{:f}'.format(res.optimum_weight) + ','
                               + '{:f}'.format(res.optimum_stress) + ','
                               + str(opti_param1) + ','
                               + str(opti_param2) + ','
                               + '{:f}'.format(res.runtime) + ','
                               + res.errorStr.replace(',', ';') + '\n')
                output_f.flush()
                jobs_done += sample_points
                print('jobs done: {:d}/{:d} -> {:f}%'.format(jobs_done, job_count, 100. * jobs_done / job_count))

    output_f.close()
    return output_file_name


def plot_sample_point_analysis(file_name, ax=None, title='', data_i=8):
    DEVIATION = 5 # $\O$ -Abweichung in $\%$
    RMSE = 6 # RMSE in $\%$
    MAE = 7
    PRESS = 8
    RIBS = 9
    SHELL = 10
    WEIGHT = 11
    STRESS = 12
    OPTI_PARAM = 13
    #data_i = RMSE
    file_path = Constants().WORKING_DIR + '/' + file_name
    data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    sampling_plan_id = data[:, 1]
    sampling_point_count = data[:, 2]
    deviation = data[:, data_i]
    sampling_data = {}
    for samp in SAMPLE_NAMES[:-1]:
        sampling_data[samp] = []
    for i in range(0, len(sampling_plan_id)):
        sampling_data[SAMPLE_NAMES[int(sampling_plan_id[i])]].append((sampling_point_count[i], deviation[i]))
    samp_plot = PlotHelper(['Anzahl der Stützstellen', 'RMSE in $\%$ von der max. Mises Spannung'], fancy=True, pgf=False, ax=ax)
    # plot one % line
    if data_i == DEVIATION or data_i == RMSE or data_i == PRESS:
        samp_plot.ax.plot([0, max(sampling_point_count)], [1., 1.], '--', color='lightgray', label='1\%-Linie')
    for key in sampling_data:
        x = [x for x,y in sampling_data[key]]
        y = [y for x,y in sampling_data[key]]
        if data_i == DEVIATION:
            y = np.array(y) * 100. # make it percent
        if data_i == RMSE or data_i == PRESS:
            y = np.array(y) * (100./max_shear_strength)
        samp_plot.ax.plot(x, y, 'x-', label=key)
    legend_loc = 'upper right'
    if data_i == DEVIATION:
        samp_plot.ax.set_ylim([0, 3.])
    elif data_i == RMSE:
        samp_plot.ax.set_ylim([0, 3.])
    elif data_i == MAE:
        samp_plot.ax.set_ylim([0, 5e+7])
    elif data_i == PRESS:
        samp_plot.ax.set_ylabel('PRESS in $\%$ von der max. Mises Spannung')
        samp_plot.ax.set_ylim([0, 5])
    elif data_i == OPTI_PARAM:
        legend_loc = 'center left'
        samp_plot.ax.set_ylabel('Polynom-Grad')
    samp_plot.ax.set_xlim([5, 30])#max(sampling_point_count)])

    if title != '':
        height = (samp_plot.ax.get_ylim()[1] - samp_plot.ax.get_ylim()[0]) * 0.75
        samp_plot.ax.text(22, height, title, fontdict=samp_plot.font)
    show_legend = True
    if ax != None:
        show_legend = False
    samp_plot.ax.locator_params(nbins=4, axis='y')
    samp_plot.ax.locator_params(nbins=7, axis='x')
    #samp_plot.finalize(height=4, width=3, legendLoc=legend_loc, show_legend=show_legend, bbox_to_anchor=(1., 0.5), tighten_layout=False)
    #samp_plot.save(Constants().PLOT_PATH + file_name.replace('csv', 'pdf'))
    if ax == None:
        samp_plot.show()
    return samp_plot


if __name__ == '__main__':
    #file = run_analysis()
    #plot_sample_point_analysis(file)
    order_plot = plot_sample_point_analysis('analysis_KrigV001.csv', data_i=9)
    order_plot.ax.locator_params(nbins=15, axis='y')
    #order_plot.save(Constants().PLOT_PATH + 'polymOrderComp.pdf')
    order_plot.show()

    plt_comp = PlotHelper([], fancy=True, pgf=True)
    import matplotlib.pyplot as plt
    ax1 = plt_comp.fig.add_subplot(322)
    ax2 = plt_comp.fig.add_subplot(324)
    ax3 = plt_comp.fig.add_subplot(326)

    poly_file = 'analysis_PolyV001.csv'
    rbf_file = 'analysis_RbfV001.csv'
    krig_file = 'analysis_KrigV001.csv'

    plot_sample_point_analysis(poly_file, ax=ax1, title='Polynom', data_i=8)
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    empty_string_labels = [''] * len(labels)
    ax1.set_xticklabels(empty_string_labels)
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    plot_sample_point_analysis(rbf_file, ax=ax2, title='RBF (mq)', data_i=8)
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    empty_string_labels = [''] * len(labels)
    ax2.set_xticklabels(empty_string_labels)
    ax2.set_xlabel('')
    plot_sample_point_analysis(krig_file, ax=ax3, title='Kriging', data_i=8)
    ax3.set_ylabel('')

    ax4 = plt_comp.fig.add_subplot(321)
    ax5 = plt_comp.fig.add_subplot(323)
    ax6 = plt_comp.fig.add_subplot(325)

    plot_sample_point_analysis(poly_file, ax=ax4, title='Polynom', data_i=6)
    labels = [item.get_text() for item in ax4.get_xticklabels()]
    empty_string_labels = [''] * len(labels)
    ax4.set_xticklabels(empty_string_labels)
    ax4.set_xlabel('')
    ax4.set_ylabel('')
    plot_sample_point_analysis(rbf_file, ax=ax5, title='RBF (mq)', data_i=6)
    labels = [item.get_text() for item in ax5.get_xticklabels()]
    empty_string_labels = [''] * len(labels)
    ax5.set_xticklabels(empty_string_labels)
    ax5.set_xlabel('')
    plot_sample_point_analysis(krig_file, ax=ax6, title='Kriging', data_i=6)
    ax6.set_ylabel('')

    handles, labels = ax1.get_legend_handles_labels()
    legend = plt_comp.fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.01), ncol=4, fancybox=True)
    import matplotlib.pyplot as plt
    plt_comp.finalize(show_legend=False, height=4.5, width=6, tighten_layout=False)
    plt.subplots_adjust(hspace=0.18, bottom=0.2, top=.98, wspace=0.35, right=.98, left=0.1)
    plt_comp.save(Constants().PLOT_PATH + 'sampleComparePressRmse.pdf')
    plt_comp.show()