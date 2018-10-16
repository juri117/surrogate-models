from timeit import default_timer as timer


class TimeTrack:

    def __init__(self, name='defaultTimer'):
        self._name = name
        self._stop = None
        self._start = timer()

    def tic(self):
        self._start = timer()

    def toc(self):
        self._stop = timer()
        print(self._name + ': ' + str(self._stop - self._start) + ' s')
        return self._stop - self._start

    def get_time(self):
        return timer() - self._start