# ==============================================================================
# description     :creates plot of openMDAO model in open_mdao.py
# author          :Juri Bieler
# date            :2018-04-19
# version         :0.01
# notes           :
# python_version  :2.7
# ==============================================================================


import subprocess
import sys


pythonPath = sys.executable
openMdaoPath = pythonPath.replace('python.exe', 'Scripts/openmdao.exe')

p = subprocess.Popen([openMdaoPath, 'view_model', 'open_mdao.py'], cwd='.')
p.wait()