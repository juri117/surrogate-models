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

    def __init__(self, project_path, half_span, box_depth, box_height, ribs, shell_thickness, box_overhang=0.):
        self.projectPath = project_path
        self.halfSpan = half_span
        self.boxDepth = box_depth
        self.boxHeight = box_height
        self.ribs = ribs
        self.shellThickness = shell_thickness
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

    def generate_wing(self, force_top, force_bot, element_size, element_type='qu4'):
        #elemType = 'qu8'
        #elemType = 'qu4'
        #force in N at wingtip
        #force = -0.5*(77000. * 9.81)
        # outer geometry in m
        #overhang = 0.5
        #height = 1.
        #halfSpan = 17.
        #boxDepth = 2.
        #elementSize = 0.25
        self.elementSize = element_size
        outLines = []
        outLines.append('# draw flat T as lines')
        outLines.append('# top right corner')
        outLines.append('pnt pc 0 0 0')
        outLines.append('seta pc all')
        outLines.append('swep pc new tra 0 0 {:f} {:d}'.format(self.boxHeight, self.calc_division(self.boxHeight)))
        outLines.append('swep pc new tra 0 {:f} 0 {:d}'.format(self.boxDepth, self.calc_division((self.boxDepth))))
        if self.boxOverhang > 0.:
            outLines.append('swep pc new tra 0 {:f} 0 {:d}'.format(-1*self.boxOverhang, self.calc_division(self.boxOverhang)))
        outLines.append('')
        outLines.append('# top left corner')
        outLines.append('pnt pc2 0 {:f} 0'.format(self.boxDepth))
        outLines.append('seta pc2 pc2')
        outLines.append('swep pc2 new tra 0 0 {:f} {:d}'.format(self.boxHeight, self.calc_division(self.boxHeight)))
        if self.boxOverhang > 0.:
            outLines.append('swep pc2 new tra 0 {:f} 0 {:d}'.format(self.boxOverhang, self.calc_division(self.boxOverhang)))
        outLines.append('')
        outLines.append('# lower right corner')
        outLines.append('pnt lowRightCorn 0 0 {:f}'.format(self.boxHeight))
        outLines.append('seta lowRightCorn lowRightCorn')
        if self.boxOverhang > 0.:
            outLines.append('swep lowRightCorn lowLeftCorn tra 0 {:f} 0 {:d}'.format(-1*self.boxOverhang, self.calc_division(self.boxOverhang)))
        outLines.append('swep lowRightCorn lowLeftCorn tra 0 {:f} 0 {:d}'.format(self.boxDepth, self.calc_division((self.boxDepth))))
        outLines.append('')
        if self.boxOverhang > 0.:
            outLines.append('# lower left corner')
            outLines.append('pnt pc4 0 {:f} {:f}'.format(self.boxDepth, self.boxHeight))
            outLines.append('seta pc4 pc4')
            outLines.append('swep pc4 new tra 0 {:f} 0 {:d}'.format(self.boxOverhang, self.calc_division(self.boxOverhang)))
            outLines.append('')
        outLines.append('seta II2d all')
        outLines.append('')

        outLines.append('# extrude the II')
        spanDiv = self.calc_span_division(self.halfSpan)
        outLines.append('swep II2d II2dOppo tra {:f} 0 0 {:d}'.format(self.halfSpan, spanDiv))
        outLines.append('# make surfaces face outside')
        outLines.append('seta toflip s A002 A005 A006')
        outLines.append('flip toflip')
        outLines.append('# new set for II beam')
        outLines.append('seta II all')
        outLines.append('')
        outLines.append('# define top and bottom shell for load')
        outLines.append('seta loadTop s A002')
        outLines.append('comp loadTop d')
        if self.boxOverhang > 0.:
            outLines.append('seta loadBot s A007')
        else:
            outLines.append('seta loadBot s A004')
        outLines.append('comp loadBot d')

        for i in range(0, self.ribs):
            if self.ribs <= 1:
                spanPos = 0
            else:
                spanPos = i * (self.halfSpan / (self.ribs - 1))
            ptName = 'rp{:d}'.format(i)
            ribName = 'rib{:d}'.format(i)
            outLines.append('')
            outLines.append('# generate a rib{:d}'.format(i))
            outLines.append('seto ' + ribName)
            #outLines.append('seta prevAll all')
            outLines.append('pnt '+ptName+' {:f} 0 0'.format(spanPos))
            #outLines.append('seta '+ribName+' p '+ribName+'')
            outLines.append('swep '+ribName+' '+ribName+' tra 0 0 {:f} {:d} a'.format(self.boxHeight, self.calc_division((self.boxHeight))))
            #outLines.append('seta '+ribName+' l all')
            #outLines.append('setr '+ribName+' l prevAll')
            outLines.append('swep '+ribName+' '+ribName+' tra 0 {:f} 0 {:d} a'.format(self.boxDepth, self.calc_division(self.boxDepth)))
            #outLines.append('comp '+ribName+' u')
            outLines.append('setc ' + ribName)

        outLines.append('')
        outLines.append('# mesh it')
        outLines.append('elty all '+element_type)
        outLines.append('mesh all')
        outLines.append('')

        outLines.append('# merge beam nodes to get one big body')
        outLines.append('seta nodes n II rib0')
        outLines.append('merg n nodes')
        outLines.append('')

        outLines.append('# write surface files for TIEs')
        outLines.append('send II abq sur')
        for i in range(1, self.ribs):
            outLines.append('seta ribL{:d} l rib{:d}'.format(i, i))

            outLines.append('comp ribL{:d} do'.format(i))
            outLines.append('comp ribL{:d} do'.format(i))
            outLines.append('send ribL{:d} abq sur'.format(i))
        outLines.append('')

        outLines.append('# write bc')
        #outLines.append('seta nodes n all')
        #outLines.append('enq nodes bc rec 0 _ _')
        outLines.append('seta bc n rib0')
        outLines.append('send bc abq nam')
        #outLines.append('seta bc2 n rib0')
        #outLines.append('setr bc2 n II2d')
        outLines.append('enq bc bc2 rec 0 _ 0.275 0.1')
        outLines.append('send bc2 abq nam')
        outLines.append('')
        #outLines.append('seta lowCorns n lowRightCorn lowLeftCorn')
        #outLines.append('send lowCorns abq nam')
        #outLines.append('send lowRightCorn abq nam')

        outLines.append('# write msh file')
        outLines.append('send all abq')
        outLines.append('')

        nodeCount = (self.calc_span_division(self.halfSpan)+1) * (self.calc_division(self.boxDepth)+1)
        if element_type == 'qu8':
            nodeCount -= 0.5*self.calc_span_division(self.halfSpan) * 0.5*self.calc_division(self.boxDepth)
        noadLoadTop = force_top/nodeCount
        noadLoadBot = force_bot / nodeCount
        outLines.append('# load application')

        #temp test with point load
        #outLines.append('seta nodes n all')
        #outLines.append('enq nodes bc1 rec {:f} 0 0 0.01'.format(self.halfSpan))
        #outLines.append('enq nodes bc2 rec {:f} {:f} 0 0.01'.format(self.halfSpan, self.boxDepth))
        #outLines.append('seta load bc1 bc2')
        #outLines.append('send load abq force 0 0 {:f}'.format((force_top+force_bot)/2.))

        outLines.append('# top')
        outLines.append('send loadTop abq force 0 0 {:f}'.format(noadLoadTop))
        outLines.append('# bottom')
        outLines.append('send loadBot abq force 0 0 {:f}'.format(noadLoadBot))
        outLines.append('')
        outLines.append('# plot it')
        outLines.append('rot -y')
        outLines.append('rot r 110')
        outLines.append('rot u 20')
        outLines.append('seta ! all')
        outLines.append('frame')
        outLines.append('zoom 2')
        outLines.append('view elem')
        outLines.append('plus n loadTop g')
        outLines.append('plus n loadBot y')
        outLines.append('plus n bc r')
        outLines.append('plus n II2dOppo b')
        outLines.append('hcpy png')
        outLines.append('sys mv hcpy_1.png mesh.png')
        outLines.append('')
        outLines.append('quit')
        f = open(self.projectPath + '/wingGeo.fbl', 'w')
        f.writelines(line + '\n' for line in outLines)
        f.close()

    def generate_inp(self, nonlinear=False):
        material_young = Constants().config.getfloat('defaults', 'material_young')
        material_poisson = Constants().config.getfloat('defaults', 'material_poisson')
        outLines = []
        outLines.append('** load mesh- and bc-file')
        outLines.append('*include, input=all.msh')
        outLines.append('*include, input=bc.nam')
        outLines.append('*include, input=bc2.nam')
        outLines.append('*include, input=II.sur')
        for i in range(1, self.ribs):
            outLines.append('*include, input=ribL{:d}.sur'.format(i))
        #outLines.append('*include, input=lowCorns.nam')
        #outLines.append('*include, input=lowRightCorn.nam')
        outLines.append('')
        outLines.append('** constraints')
        outLines.append('*boundary')
        outLines.append('Nbc,1')
        outLines.append('Nbc2,1,6')
        #outLines.append('NlowCorns,3')
        #outLines.append('NlowRightCorn,2')
        outLines.append('')
        outLines.append('** material definition')
        outLines.append('*MATERIAL,NAME=alu')
        outLines.append('*ELASTIC')
        outLines.append('{:.3e}, {:.3f}'.format(material_young, material_poisson))
        outLines.append('')
        outLines.append('** define surfaces')
        outLines.append('*shell section, elset=Eall, material=alu')
        outLines.append('{:f}'.format(self.shellThickness))
        outLines.append('')
        for i in range(1, self.ribs):
            outLines.append('*tie,name=t{:d},position tolerance=0.1'.format(i))
            outLines.append('SribL{:d},SII'.format(i))
            outLines.append('')
        #outLines.append('')
        outLines.append('** step')
        if nonlinear:
            outLines.append('*step, nlgeom')
        else:
            outLines.append('*step')
        outLines.append('*static')
        outLines.append('')
        outLines.append('** load')
        outLines.append('*cload')
        #outLines.append('*include, input=load.frc')
        outLines.append('*include, input=loadTop.frc')
        outLines.append('*include, input=loadBot.frc')
        outLines.append('')
        outLines.append('*node file')
        outLines.append('U')
        outLines.append('*el file')
        outLines.append('S')
        outLines.append('*end step')
        outLines.append('')
        f = open(self.projectPath + '/wingRun.inp', 'w')
        f.writelines(line + '\n' for line in outLines)
        f.close()


if __name__ == '__main__':
    geo = WingConstruction()
    geo.generate_wing('../dataOut/test01/test', 5, 0.1, 0.1, 0.1)
    geo.generate_inp('../dataOut/test01/test', 2.)