__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :interface to Calculix and solver-input-file generation
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

from wingConstruction.utils.Constants import Constants


class WingConstruction:

    def __init__(self, project_path, half_span, box_depth, box_height, ribs, shell_thickness, engine_pos, box_overhang=0.):
        self.projectPath = project_path
        self.halfSpan = half_span
        self.boxDepth = box_depth
        self.boxHeight = box_height
        self.ribs = ribs
        self.shellThickness = shell_thickness
        self.enginePos = engine_pos
        self.pylonHeight = 0.3
        self.boxOverhang = box_overhang
        self.elementSize = 0.1

    def calc_weight(self, density):
        v_box = self.halfSpan * 2. * (self.boxHeight + self.boxDepth) * self.shellThickness
        v_ribs = self.ribs * self.boxHeight * self.boxDepth * self.shellThickness
        w = (v_box + v_ribs) * density
        return w

    def calc_span_division(self, length):
        return self.calc_division(length)
        #div = int(length / self.elementSize)
        #if self.ribs <= 1:
        #    return self.calc_division(length)
        #while div % (self.ribs - 1) > 0 or div % 2 > 0:
        #    div += 1
        #return div

    # the division shouldn't be odd or 0
    def calc_division(self, length):
        div = int(length / self.elementSize)
        if div == 0:
            div = 2
        if div % 2 > 0:
            div += 1
        return div

    def generate_wing(self, force_top, force_bot, engine_weight, element_size, element_type='qu4'):
        self.elementSize = element_size
        out_lines = []
        out_lines.append('# draw flat T as lines')
        out_lines.append('# top right corner')
        out_lines.append('pnt pc 0 0 0')
        out_lines.append('seta pc all')
        out_lines.append('swep pc new tra 0 0 {:f} {:d}'.format(self.boxHeight, self.calc_division(self.boxHeight)))
        out_lines.append('swep pc new tra 0 {:f} 0 {:d}'.format(self.boxDepth, self.calc_division((self.boxDepth))))
        if self.boxOverhang > 0.:
            out_lines.append('swep pc new tra 0 {:f} 0 {:d}'.format(-1*self.boxOverhang, self.calc_division(self.boxOverhang)))
        out_lines.append('')
        out_lines.append('# top left corner')
        out_lines.append('pnt pc2 0 {:f} 0'.format(self.boxDepth))
        out_lines.append('seta pc2 pc2')
        out_lines.append('swep pc2 new tra 0 0 {:f} {:d}'.format(self.boxHeight, self.calc_division(self.boxHeight)))
        if self.boxOverhang > 0.:
            out_lines.append('swep pc2 new tra 0 {:f} 0 {:d}'.format(self.boxOverhang, self.calc_division(self.boxOverhang)))
        out_lines.append('')
        out_lines.append('# lower right corner')
        out_lines.append('pnt lowRightCorn 0 0 {:f}'.format(self.boxHeight))
        out_lines.append('seta lowRightCorn lowRightCorn')
        if self.boxOverhang > 0.:
            out_lines.append('swep lowRightCorn lowLeftCorn tra 0 {:f} 0 {:d}'.format(-1*self.boxOverhang, self.calc_division(self.boxOverhang)))
        out_lines.append('swep lowRightCorn lowLeftCorn tra 0 {:f} 0 {:d}'.format(self.boxDepth, self.calc_division(self.boxDepth)))
        out_lines.append('')
        if self.boxOverhang > 0.:
            out_lines.append('# lower left corner')
            out_lines.append('pnt pc4 0 {:f} {:f}'.format(self.boxDepth, self.boxHeight))
            out_lines.append('seta pc4 pc4')
            out_lines.append('swep pc4 new tra 0 {:f} 0 {:d}'.format(self.boxOverhang, self.calc_division(self.boxOverhang)))
            out_lines.append('')

        '''
        out_lines.append('# stringer')
        out_lines.append('pnt str0 0 0.3 0')
        out_lines.append('seta str0 str0')
        out_lines.append('swep str0 str0 tra 0 0 0.005 2')
        out_lines.append('')
        out_lines.append('pnt str1 0 0.9 0')
        out_lines.append('seta str1 str1')
        out_lines.append('swep str1 str1 tra 0 0 0.005 2')
        out_lines.append('')
        
        out_lines.append('pnt str2 0 0.3 0.55')
        out_lines.append('seta str2 str2')
        out_lines.append('swep str2 str2 tra 0 0 -0.005 2')
        out_lines.append('')
        out_lines.append('pnt str3 0 0.9 0.55')
        out_lines.append('seta str3 str3')
        out_lines.append('swep str3 str3 tra 0 0 -0.005 2')
        out_lines.append('')
        '''

        out_lines.append('seta II2d all')
        out_lines.append('')

        out_lines.append('# extrude the II')
        spanDiv = self.calc_span_division(self.halfSpan)
        out_lines.append('swep II2d II2dOppo tra {:f} 0 0 {:d}'.format(self.halfSpan, spanDiv))
        out_lines.append('# make surfaces face outside')
        out_lines.append('seta toflip s A002 A005 A006')
        out_lines.append('flip toflip')
        out_lines.append('# new set for II beam')
        out_lines.append('seta II all')
        out_lines.append('')
        out_lines.append('# define top and bottom shell for load')
        out_lines.append('seta loadTop s A002')
        out_lines.append('comp loadTop d')
        if self.boxOverhang > 0.:
            out_lines.append('seta loadBot s A007')
        else:
            out_lines.append('seta loadBot s A004')
        out_lines.append('comp loadBot d')

        out_lines.append('#generate engine pylon')
        out_lines.append('seto pyl')
        out_lines.append('pnt pylP {:f} 0 {:f}'.format(self.enginePos, self.boxHeight))
        out_lines.append('swep pyl pyl tra 0 0 {:f} 6 a'.format(self.pylonHeight))
        out_lines.append('swep pyl pyl tra 0 {:f} 0 {:d} a'.format(self.boxDepth, self.calc_division(self.boxDepth)))
        out_lines.append('setc pyl')
        out_lines.append('')

        for i in range(0, self.ribs):
            if self.ribs <= 1:
                span_pos = 0
            else:
                span_pos = i * (self.halfSpan / (self.ribs - 1))
            pt_name = 'rp{:d}'.format(i)
            rib_name = 'rib{:d}'.format(i)
            out_lines.append('')
            out_lines.append('# generate a rib{:d}'.format(i))
            out_lines.append('seto ' + rib_name)
            out_lines.append('pnt '+pt_name+' {:f} 0 0'.format(span_pos))
            out_lines.append('swep '+rib_name+' '+rib_name+' tra 0 0 {:f} {:d} a'.format(self.boxHeight, self.calc_division((self.boxHeight))))
            out_lines.append('swep '+rib_name+' '+rib_name+' tra 0 {:f} 0 {:d} a'.format(self.boxDepth, self.calc_division(self.boxDepth)))
            out_lines.append('setc ' + rib_name)

        out_lines.append('')
        out_lines.append('# mesh it')
        out_lines.append('elty all '+element_type)
        out_lines.append('mesh all')
        out_lines.append('')

        out_lines.append('# merge beam nodes to get one big body')
        out_lines.append('seta nodes n II rib0')
        out_lines.append('merg n nodes')
        out_lines.append('')

        out_lines.append('# write surface files for TIEs')
        out_lines.append('send II abq sur')
        out_lines.append('')
        out_lines.append('seta pylL l pyl')
        out_lines.append('comp pylL do')
        out_lines.append('comp pylL do')
        out_lines.append('send pylL abq sur')
        out_lines.append('')

        for i in range(1, self.ribs):
            out_lines.append('seta ribL{:d} l rib{:d}'.format(i, i))

            out_lines.append('comp ribL{:d} do'.format(i))
            out_lines.append('comp ribL{:d} do'.format(i))
            out_lines.append('send ribL{:d} abq sur'.format(i))
        out_lines.append('')

        out_lines.append('# write bc')
        out_lines.append('seta bc n rib0')
        out_lines.append('send bc abq nam')
        out_lines.append('enq bc bc2 rec 0 _ 0.275 0.1')
        out_lines.append('send bc2 abq nam')
        out_lines.append('')

        out_lines.append('# write msh file')
        out_lines.append('send all abq')
        out_lines.append('')

        node_count = (self.calc_span_division(self.halfSpan)+1) * (self.calc_division(self.boxDepth)+1)
        node_count_engine = self.calc_division(self.boxDepth)+1
        if element_type == 'qu8':
            node_count -= 0.5*self.calc_span_division(self.halfSpan) * 0.5*self.calc_division(self.boxDepth)
        noad_load_top = force_top/node_count
        noad_load_bot = force_bot / node_count
        node_load_engine = engine_weight / node_count_engine
        out_lines.append('# load application')

        out_lines.append('# top')
        out_lines.append('send loadTop abq force 0 0 {:f}'.format(noad_load_top))
        out_lines.append('# bottom')
        out_lines.append('send loadBot abq force 0 0 {:f}'.format(noad_load_bot))
        out_lines.append('#engine weight')
        out_lines.append('seta engNodes n pyl')
        out_lines.append('enq engNodes engLoad rec {:f} _ {:f} 0.01'.format(self.enginePos, self.boxHeight+self.pylonHeight))
        out_lines.append('send engLoad abq force 0 0 {:f}'.format(node_load_engine))
        out_lines.append('')
        out_lines.append('')

        out_lines.append('# plot it')
        out_lines.append('rot -y')
        out_lines.append('rot r 110')
        out_lines.append('rot u 20')
        out_lines.append('seta ! all')
        out_lines.append('frame')
        out_lines.append('zoom 2')
        out_lines.append('view elem')
        out_lines.append('plus n loadTop g')
        out_lines.append('plus n loadBot y')
        out_lines.append('plus n bc r')
        out_lines.append('plus n II2dOppo b')
        out_lines.append('hcpy png')
        out_lines.append('sys mv hcpy_1.png mesh.png')
        out_lines.append('')
        out_lines.append('quit')
        f = open(self.projectPath + '/wingGeo.fbl', 'w')
        f.writelines(line + '\n' for line in out_lines)
        f.close()

    def generate_inp(self, nonlinear=False):
        material_young = Constants().config.getfloat('defaults', 'material_young')
        material_poisson = Constants().config.getfloat('defaults', 'material_poisson')
        out_lines = []
        out_lines.append('** load mesh- and bc-file')
        out_lines.append('*include, input=all.msh')
        out_lines.append('*include, input=bc.nam')
        out_lines.append('*include, input=bc2.nam')
        out_lines.append('*include, input=II.sur')
        out_lines.append('*include, input=pylL.sur')
        for i in range(1, self.ribs):
            out_lines.append('*include, input=ribL{:d}.sur'.format(i))
        out_lines.append('')
        out_lines.append('** constraints')
        out_lines.append('*boundary')
        out_lines.append('Nbc,1')
        out_lines.append('Nbc2,1,6')
        out_lines.append('')
        out_lines.append('** material definition')
        out_lines.append('*MATERIAL,NAME=ALU')
        out_lines.append('*ELASTIC')
        out_lines.append('{:.3e}, {:.3f}'.format(material_young, material_poisson))
        out_lines.append('')
        out_lines.append('** define surfaces')
        out_lines.append('*shell section, elset=Eall, material=ALU')
        out_lines.append('{:f}'.format(self.shellThickness))
        out_lines.append('')
        out_lines.append('*tie,name=t1,position tolerance=0.1')
        out_lines.append('SpylL,SII')
        out_lines.append('')
        for i in range(1, self.ribs):
            out_lines.append('*tie,name=t{:d},position tolerance=0.1'.format(i))
            out_lines.append('SribL{:d},SII'.format(i))
            out_lines.append('')
        out_lines.append('** step')
        if nonlinear:
            out_lines.append('*step, nlgeom')
        else:
            out_lines.append('*step')
        out_lines.append('*static')
        out_lines.append('')
        out_lines.append('** load')
        out_lines.append('*cload')
        out_lines.append('*include, input=loadTop.frc')
        out_lines.append('*include, input=loadBot.frc')
        out_lines.append('*include, input=engLoad.frc')
        out_lines.append('')
        out_lines.append('*node file')
        out_lines.append('U')
        out_lines.append('*el file')
        out_lines.append('S')
        out_lines.append('*end step')
        out_lines.append('')
        f = open(self.projectPath + '/wingRun.inp', 'w')
        f.writelines(line + '\n' for line in out_lines)
        f.close()


if __name__ == '__main__':
    geo = WingConstruction()
    geo.generate_wing('../dataOut/test01/test', 5, 0.1, 0.1, 0.1)
    geo.generate_inp('../dataOut/test01/test', 2.)