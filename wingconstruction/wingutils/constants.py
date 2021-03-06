__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :singleton holding global parameters like paths and setting from data_in/setup.ini
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import os.path
import sys
import configparser

from wingconstruction.wingutils.singleton import Singleton


# this is needed for python 2, 3 compatibility
def with_metaclass(mcls):
    def decorator(cls):
        body = vars(cls).copy()
        # clean out class body
        body.pop('__dict__', None)
        body.pop('__weakref__', None)
        return mcls(cls.__name__, cls.__bases__, body)
    return decorator


@with_metaclass(Singleton)
class Constants(object):
    __metaclass__ = Singleton

    # tolerance for float comparison
    FLOAT_TOLERANCE = 0.00000001
    PROJECT_ROOT_DIR = os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + '/../../')
    INPUT_DIR = PROJECT_ROOT_DIR + '/data_in'

    def __init__(self):
        print('init Constants...')
        self.config = configparser.ConfigParser()
        self.config.read(self.INPUT_DIR+'/setup.ini')

        ### PATHS
        self.WORKING_DIR = self.PROJECT_ROOT_DIR + '/' + self.config['meta']['working_dir']
        if not os.path.isdir(self.WORKING_DIR):
            os.mkdir(self.WORKING_DIR)

        self.PLOT_PATH = self.WORKING_DIR + '\\plot'
        if not os.path.isdir(self.PLOT_PATH):
            os.mkdir(self.PLOT_PATH)
        self.PLOT_PATH = self.PLOT_PATH + '\\'

        self.CALCULIX_BIN_PATH = self.config['fem']['calculix_path']
        if not os.path.isdir(self.CALCULIX_BIN_PATH):
            print('ERROR: calculix could not be found at given location: ' + self.CALCULIX_BIN_PATH)
            print('open the file: ' + self.SETUP_INI_FILE_NAME + ' and adapt the entry for calculix_path to the real location')
            sys.exit(1)

        self.ABAQUS_EXE_PATH = self.config['fem']['abaqus_exe_path']

        self.CALCULIX_CCX_EXE_PATH = self.CALCULIX_BIN_PATH + '/' + self.config['fem']['calculix_ccx_executable']
        self.CALCULIX_CGX_EXE_PATH = self.CALCULIX_BIN_PATH + '/' + self.config['fem']['calculix_cgx_executable']
