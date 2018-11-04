__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :abstract of subprocess to support real time output (work in progress...)
# date            :2018-01-11
# notes           :
# python_version  :3.6
# ==============================================================================

import subprocess
import os
import signal
import time
import sys
from threading import Thread
from timeit import default_timer as timer

try:
    from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty  # python 3.x


class RealTimeSubprocess:

    def __init__(self):
        self.useTimeout = False
        print('done')

    def wait_for_keyword(self, que, word):
        while (True):
            try:
                line = que.get_nowait()  # or q.get(timeout=.1)
                if sys.version_info[0] >= 3:
                    line = line.decode('UTF-8')
            except Empty:
                #print('no output yet')
                time.sleep(0.01)
            else:  # got line
                print(line)
                if word in line:
                    return True

    def wait(self):
        while not self.terminated:
            if not self.timeout <= 0.:
                if timer() - self.startTime > self.timeout:
                    os.kill(os.getpgid(self.process.pid), signal.SIGTERM)
                    print('kill')
            print(str(timer() - self.startTime))
            time.sleep(0.001)

    def enqueue_output(self, out, queue):
        for line in iter(out.readline, b''):
            try:
                queue.put(line.decode('UTF-8'))
                print(line.decode('UTF-8'))
            except:
                queue.put(line)
                print(line)
        self.terminated = True
        print('closing')
        out.close()

    def write_to_console_and_enter(self, p, str):
        outStr = str + '\n'
        if sys.version_info[0] >= 3:
            outStr = outStr.encode('UTF-8')
        p.stdin.write(outStr)
        p.stdin.flush()

    def execute(self, commandList, working_dir='.', timeout=0.):
        self.errorFlag = False
        self.terminated = False
        self.startTime = timer()
        self.timeout = timeout
        ON_POSIX = 'posix' in sys.builtin_module_names
        self.process = subprocess.Popen(commandList,
                             cwd=working_dir,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             bufsize=1,
                             close_fds=ON_POSIX)
        self.q = Queue()
        t = Thread(target=self.enqueue_output, args=(self.process.stdout, self.q))
        #t.daemon = True  # thread dies with the program
        t.start()


if __name__ == '__main__':
    prc = RealTimeSubprocess()
    prc.execute(['D:/prog/portable/Luftfahrt/bConverged/CalculiX/bin/cgx.bat',
                 '-b',
                 'wing_post.fbd'],
                working_dir='../../data_out/firstTry',
                timeout=.1)
    #while not prc.terminated:
    #    time.sleep(0.01)
    prc.wait()
    print('done')