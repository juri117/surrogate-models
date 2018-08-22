__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :interface to Abaqus and solver-input-file generation
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import subprocess
import os
from shutil import copyfile

from wingConstruction.utils.Constants import Constants

##############################################
# set up variables

#SHELL_THICKNESS = 0.008 #in mm
#SHELL_MATERIAL_YOUNG = 60e9 #N/m^2
#SHELL_MATERIAL_POISSON = 0.34

class Abaqus():

    ##############################################
    # system variables [DO NOT TOUCH]
    errorFlag = False

    # init function prepares working environment
    def __init__(self, workingDir='data'):
        self._workingDir = workingDir
        if not os.path.isdir(self._workingDir):
            os.mkdir(self._workingDir)
        self.dispD3Min = 0
        self.dispD3Max = 0
        self.stressMisesMin = 0
        self.stressMisesMax = 0
        self.stressMisesMaxFixed = 0

    ##############################################
    # fem pre-processing

    def calculix_to_abaqus(self, shell_thickness):
        # open a new file
        aba_f = open(self._workingDir + '/' + 'abaqusJob.inp', 'w')
        aba_f.write('*Heading\n')
        aba_f.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
        aba_f.write('*Part, name=PART-1\n')
        # copy the mesh
        mesh_f = open(self._workingDir + '/' + 'all.msh', 'r')
        node_ids = []
        node_section = False
        for l in mesh_f.readlines():
            if '*NODE' in l:
                aba_f.write('*Node\n')
                node_section = True
            elif '*ELEMENT' in l:
                if 'S4' in l:
                    aba_f.write('*Element, type=S4\n')
                else:
                    print('ERROR, unknown element in mesh file, adapt parser...')
                node_section = False
            elif l != '':
                aba_f.write(l)
                if node_section:
                    vals = l.split(',')
                    node_ids.append(int(vals[0]))
        mesh_f.close()
        max_elem_id = int(l.split(',')[0])
        # create elset and define shell
        aba_f.write('*Elset, elset=EALL, generate\n')
        aba_f.write('\t{:d},\t{:d},\t1\n'.format(1, max_elem_id))
        aba_f.write('*Shell Section, elset=EALL, material=ALU\n')
        aba_f.write('{:f},{:d}\n'.format(shell_thickness, 5))
        aba_f.write('*End Part\n')
        # start assembly
        aba_f.write('*Assembly, name=Assembly\n')
        aba_f.write('*Instance, name=PART-1-1, part=PART-1\n')
        aba_f.write('*End Instance\n')
        # generate NALL set
        aba_f.write('*Nset, nset=NALL, instance=PART-1-1\n')
        for n_id in node_ids:
            aba_f.write('{:d},\n'.format(n_id))
        # load all node sets
        for file in os.listdir(self._workingDir):
            if file.endswith('.nam'):
                nam_f = open(self._workingDir + '/' + file, 'r')
                for l in nam_f.readlines():
                    if '**' in l:
                        pass
                    elif '*NSET' in l:
                        n_set_name = (l.split('=')[1]).strip().upper()
                        aba_f.write('*Nset, nset='+n_set_name+', instance=PART-1-1\n')
                    elif l != '':
                        aba_f.write(l)
        #load all frc files containing loads and nodes
        loads = []
        for file in os.listdir(self._workingDir):
            if file.endswith('.frc'):
                nam_f = open(self._workingDir + '/' + file, 'r')
                file_name = (file.split('/')[-1]).split('.')[0]
                aba_f.write('*Nset, nset = '+file_name+', internal, instance = PART-1-1\n')
                for l in nam_f.readlines():
                    if '*' in l:
                        pass
                    elif l != '':
                        vals = l.split(',')
                        aba_f.write(vals[0] + ',\n')
                loads.append(Load(file_name, int(vals[1]), float(vals[2])))
        # load all surface sets
        surf_names = []
        for file in os.listdir(self._workingDir):
            if file.endswith('.sur'):
                nam_f = open(self._workingDir + '/' + file, 'r')
                file_name = (file.split('/')[-1]).split('.')[0]
                aba_f.write('*Nset, nset = '+file_name+', internal, instance = PART-1-1\n')
                min_id = 9999999999
                max_id = 0
                for l in nam_f.readlines():
                    if '**' in l:
                        pass
                    elif '*SURFACE' in l:
                        surf_name = (l.split('=')[-1]).strip().upper()
                        aba_f.write('*Elset, elset=_'+surf_name+'_SPOS, internal, instance=PART-1-1, generate\n')
                    elif l != '':
                        vals = l.split(',')
                        min_id = min(int(vals[0]), min_id)
                        max_id = max(int(vals[0]), max_id)
                aba_f.write('\t{:d},\t{:d},\t{:d}\n'.format(min_id, max_id, 1))
                aba_f.write('*Surface, type=ELEMENT, name='+surf_name+'\n')
                aba_f.write('_'+surf_name+'_SPOS, SPOS\n')
                if not 'SII' in surf_name:
                    surf_names.append(surf_name)
        # write constraints
        for n in surf_names:
            aba_f.write('*Tie, name='+n+'-1, adjust=yes, position tolerance=0.01\n')
            aba_f.write(n+', SII\n')
        aba_f.write('*End Assembly\n')
        # write materials
        material_young = Constants().config.getfloat('defaults', 'material_young')
        material_poisson = Constants().config.getfloat('defaults', 'material_poisson')
        aba_f.write('*Material, name=ALU\n')
        aba_f.write('*Elastic\n')
        aba_f.write(' {:.3e}, {:.3f}\n'.format(material_young, material_poisson))
        # boundary conditions
        aba_f.write('*Boundary\n')
        aba_f.write('NBC, 1, 1\n')
        aba_f.write('*Boundary\n')
        aba_f.write('NBC2, 1, 1\n')
        aba_f.write('NBC2, 2, 2\n')
        aba_f.write('NBC2, 3, 3\n')
        # step
        aba_f.write('*Step, name=Step-1, nlgeom=NO\n')
        aba_f.write('*Static\n')
        aba_f.write('1., 1., 1e-05, 1.\n')
        # loads
        for load in loads:
            aba_f.write('*Cload\n')
            aba_f.write(load.name+', {:d}, {:f}\n'.format(load.direction, load.val))
        # get done already
        aba_f.write('*Restart, write, frequency=0\n')
        aba_f.write('*Output, field, variable=PRESELECT\n')
        aba_f.write('*Output, history, variable=PRESELECT\n')
        aba_f.write('*End Step\n')
        aba_f.close()


    ##############################################
    # fem post-processing

    def post_processing(self, save_to_file=False):
        print('run fem solver ccx(' + self._workingDir + ')')
        p = subprocess.Popen([Constants().ABAQUS_EXE_PATH, 'odbreport', 'odb=abaqusJob', 'mode=CSV', 'results', 'invariants'],
                             cwd=self._workingDir,
                             stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        out = out.decode('UTF-8')
        if save_to_file:
            out_f = open(self._workingDir + '/results.txt', 'w')
            out_f.write(out)
            out_f.close()
        print('done')
        lines = out.split('\n')
        for i in range(0, len(lines)):
            if 'Invariants of field \'S\'' in lines[i]:
                i += 2
                break

        misis = []
        for i2 in range(i, len(lines)):
            if lines[i2] == '':
                break
            vals = lines[i2].split(',')
            if len(vals) != 15:
                break
            misis.append(float(vals[5].strip()))

        self.stressMisesMax = max(misis)
        self.stressMisesMin = min(misis)




    ##############################################
    # fem calls

    # calls the fem solver (all input files must be present in the working directory)
    def solve_model(self):
        # print('--- start ccx output ---------------------------------------')
        print('run fem solver ccx('+self._workingDir+')')
        p = subprocess.Popen([Constants().ABAQUS_EXE_PATH, 'job=abaqusJob', 'cpus={:d}'.format(Constants().config.getint('meta', 'used_cores')), 'int', 'ask=off'], cwd=self._workingDir,
                             stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        out = out.decode('UTF-8')
        print(out)
        if not 'COMPLETED' in out:
            self.errorFlag = True


    def process_cgx_output(self, strOut):
        lines = strOut.split(b'\n')
        for i in range(0, len(lines) - 3):
            if b'DISP' in lines[i] and b'D3' in lines[i]:
                if b'max:' in lines[i + 2]:
                    str = lines[i + 2].replace(b' max:', b'')
                    str = str.split(b' ')[0]
                    self.dispD3Max = float(str)
                if b'min:' in lines[i + 3]:
                    str = lines[i + 3].replace(b' min:', b'')
                    str = str.split(b' ')[0]
                    self.dispD3Min = float(str)
            if b'STRESS' in lines[i] and b'Mises' in lines[i]:
                if b'max:' in lines[i + 2]:
                    str = lines[i + 2].replace(b' max:', b'')
                    str = str.split(b' ')[0]
                    self.stressMisesMax = float(str)
                if b'min:' in lines[i + 3]:
                    str = lines[i + 3].replace(b' min:', b'')
                    str = str.split(b' ')[0]
                    self.stressMisesMin = float(str)
            if b'mode:h' in lines[i]:
                if b'node' in lines[i+1] and b'value' in lines[i+1]:
                    str = lines[i + 1].replace(b' node:', b'')
                    str = str.replace(b'value:', b'')
                    str = str.replace(b'dist:', b'')
                    vals = str.split(b' ')
                    self.stressMisesMaxFixed = float(vals[1])

    ##############################################
    # helper functions

    def run_ccx(self, fileName, pipe_response=False):
        if pipe_response:
            p = subprocess.Popen([Constants().CALCULIX_CCX_EXE_PATH, fileName], cwd=self._workingDir,
                                 stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            p = subprocess.Popen([Constants().CALCULIX_CCX_EXE_PATH, fileName], cwd=self._workingDir)
            p.wait()
        return p

    def run_cgx(self, file_name, pipe_response=False):
        if pipe_response:
            p = subprocess.Popen([Constants().CALCULIX_CGX_EXE_PATH, '-b', file_name], cwd=self._workingDir,
                                 stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        else:
            p = subprocess.Popen([Constants().CALCULIX_CGX_EXE_PATH, '-b', file_name], cwd=self._workingDir)
            p.wait()
        return p

    '''
    def run_cgx_postprocessing(self, partName):
        p = subprocess.Popen([self.pathCGX, '-b', partName + '.fbd'], cwd=self.workingDir)
        p.wait()
        print("done")
        return p
    '''

    def get_file_path(self, fileName):
        return self._workingDir + '/' + fileName


class Load:
    def __init__(self, name, direction, val):
        self.name = name
        self.direction = direction
        self.val = val


if __name__ == '__main__':
    clx = Calculix(workingDir='../dataOut/test01')
    clx.generate_mesh('test')
    #clx.generate_geometry()
    #clx.generate_mesh()
    #clx.generate_bc()
    #clx.generate_load()
    #clx.generate_inp()
    #clx.solve_model()
