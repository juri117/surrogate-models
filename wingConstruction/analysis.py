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

from wingConstruction.utils.defines import *
from wingConstruction.surrogateV2 import surrogate_analysis
from wingConstruction.surrogateV2 import SurroResults
from wingConstruction.utils.Constants import Constants

if __name__ == '__main__':
    surroMethods = [KRIGING] # KRIGING, RBF
    sampleMethods = [HALTON] # LATIN, HAMMERS, HALTON
    samplePointCount = [14]#[10, 11, 12, 13, 14, 15, 16]
    useAbaqus = False
    usePGF = False

    output_file_name = 'surro_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
    output_f = open(Constants().WORKING_DIR + '/'
                    + output_file_name,
                    'w')
    output_f.write('SampleMethod,SamplePointCound,SurroMehtod,deviation,optRib,optShell,optWight,runtime\n')

    for surroM in surroMethods:
        for sampleM in sampleMethods:
            for samplePoints in samplePointCount:
                res = surrogate_analysis(sampleM,
                                         samplePoints,
                                         surroM,
                                         use_abaqus=useAbaqus,
                                         pgf=usePGF,
                                         show_plots=False)
                output_f.write(SAMPLE_NAMES[sampleM] + ','
                               + str(samplePoints) + ','
                               + SURRO_NAMES[surroM] + ','
                               + '{:f}'.format(res.deviation) + ','
                               + '{:f}'.format(res.optimumRib) + ','
                               + '{:f}'.format(res.optimumShell) + ','
                               + '{:f}'.format(res.optimumWeights) + ','
                               + '{:f}'.format(res.runtime) + '\n')
    output_f.close()
