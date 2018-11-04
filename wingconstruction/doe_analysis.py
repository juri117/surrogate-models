__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :runs simple DoE analysis
# author          :Juri Bieler
# date            :2018-09-10
# notes           :
# python_version  :3.6
# ==============================================================================


from mylibs.doe import DoE
from wingconstruction.wingutils.defines import *
from myutils.plot_helper import PlotHelper


input_names = ['ribs', 'shell']

ranges = [range_rib, range_shell]
from wingconstruction.multi_run import MultiRun
runner = MultiRun(use_calcu=True, use_aba=False, non_liner=False, force_recalc=False, project_name_prefix='DoE')
d2 = DoE(input_names, ranges, runner.calc_stress, level_count=2)
d2.corellation()
d2.print_res_table(ref=max_shear_strength)

#d3 = DoE(input_names, ranges, runner.calc_stress, level_count=3)
#d3.corellation()
#d3.print_res_table()

if True:
    FANCY = True
    PGF = False

    avg = 0.5 * (d2._results[0].res + d2._results[3].res) - 0.5 * (d2._results[1].res + d2._results[2].res)
    print('interaction: {:f}'.format(avg))


    pl = PlotHelper([], fancy=FANCY, pgf=PGF)
    import matplotlib.pyplot as plt

    ax1 = pl.fig.add_subplot(221)
    ax2 = pl.fig.add_subplot(223)
    ax3 = pl.fig.add_subplot(224)

    pl1 = PlotHelper(['Level', 'Ausgang'], fancy=FANCY, pgf=PGF, ax=ax1)
    line_b = pl1.ax.plot([-1, 1], [d2._results[0].res, d2._results[1].res], '-', label='Blech-Einfl.(Rippen-L.: -1)')
    line_r = pl1.ax.plot([-1, 1], [d2._results[0].res, d2._results[2].res], '-', label='Rippen-Einfl.(Blech-L.: -1)')
    pl1.finalize(show_legend=False, legendLoc='upper right', bbox_to_anchor=(1.2, 1.))

    pl2 = PlotHelper(['Level', 'Ausgang'], fancy=FANCY, pgf=PGF, ax=ax2)
    pl2.ax.plot([-1, 1], [d2._results[0].res, d2._results[1].res], '-', color=line_b[0].get_color(), label='Blech-Einfl.(Rippen-L.: -1)')
    pl2.ax.plot([-1, 1], [d2._results[2].res, d2._results[3].res], '--', color=line_b[0].get_color(), label='Blech-Einfl.(Rippen-L.: +1)')

    pl3 = PlotHelper(['Level', 'Ausgang'], fancy=FANCY, pgf=PGF, ax=ax3)
    pl3.ax.plot([-1, 1], [d2._results[0].res, d2._results[2].res], '-', color=line_r[0].get_color(), label='Rippen-Einfl.(Blech-L.: -1)')
    pl3.ax.plot([-1, 1], [d2._results[1].res, d2._results[3].res], '--', color=line_r[0].get_color(), label='Rippen-Einfl.(Blech-L.: +1)')

    handles2, labels2 = ax2.get_legend_handles_labels()
    handles3, labels3 = ax3.get_legend_handles_labels()
    legend = pl.fig.legend(handles2 + handles3, labels2 + labels3, loc='upper right', ncol=1, fancybox=True, bbox_to_anchor=(.97, .94))
    pl.finalize(height=3.5, show_legend=False)

    #pl.ax.plot([-1, 0, 1], [d3._results[0].res, d3._results[1].res, d3._results[2].res], label='Blechdickeneinfluss')
    #pl.ax.plot([-1, 0, 1], [d3._results[0].res, d3._results[3].res, d3._results[6].res], label='Rippenanzahleinfluss')
    #pl.ax.plot([-1, 0, 1], [d3._results[3].res, d3._results[4].res, d3._results[5].res],
    #           label='Blechdickeneinfluss')
    #pl.ax.plot([-1, 0, 1], [d3._results[1].res, d3._results[4].res, d3._results[7].res],
    #           label='Rippenanzahleinfluss')
    #pl.ax.plot([-1, 0, 1], [d3._results[6].res, d3._results[7].res, d3._results[8].res],
    #           label='Blechdickeneinfluss')
    #pl.ax.plot([-1, 0, 1], [d3._results[2].res, d3._results[5].res, d3._results[8].res],
    #           label='Rippenanzahleinfluss')

    from wingconstruction.wingutils.constants import Constants
    pl.save(Constants().PLOT_PATH + 'wingDoE.pdf')
    pl.finalize()
    pl.show()
