__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :helper class for tracking runntime
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================


from timeit import default_timer as timer


class TimeTrack:

    def __init__(self, name='defaultTimer'):
        self._name = name
        self._stop = None
        self._start = timer()

    def tic(self):
        self._start = timer()

    def toc(self, print_it=False):
        self._stop = timer()
        if print_it:
            self.print_time()
        return self._stop - self._start

    def print_time(self):
        print(self._name + ': ' + str(self._stop - self._start) + ' s')

    def get_time(self):
        return timer() - self._start
