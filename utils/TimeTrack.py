from timeit import default_timer as timer


class TimeTrack:

    def __init__(self, name='defaultTimer'):
        self.name = name
        self.stop = None
        self.start = timer()

    def tic(self):
        self.start = timer()

    def toc(self):
        self.stop = timer()
        print(self.name + ': ' + str(self.stop - self.start) + ' s')
        return self.stop - self.start