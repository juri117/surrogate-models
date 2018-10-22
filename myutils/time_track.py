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