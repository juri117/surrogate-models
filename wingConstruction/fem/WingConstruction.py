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


class WingConstruction:

    def __init__(self, project_path, half_span, box_depth, box_height, n_ribs, box_overhang=0.):
        self.projectPath = project_path
        self.halfSpan = half_span
        self.boxDepth = box_depth
        self.boxHeight = box_height
        self.nRibs = n_ribs
        self.boxOverhang = box_overhang
        self.elementSize = 0.1
        print('done')

    def calc_span_division(self, length):
        div = int(length / self.elementSize)
        if self.nRibs <= 1:
            return self.calc_division(length)
        while div % (self.nRibs-1) > 0 or div % 2 > 0:
            div += 1
        return div

    # the division shouldn't be odd or 0
    def calc_division(self, length):
        div = int(length / self.elementSize)
        if div == 0:
            div = 2
        if div % 2 > 0:
            div += 1
        return div

    def generate_wing(self, force_top, force_but, element_size, element_type='qu4'):
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
        outLines.append('pnt pc 0 0 0')
        outLines.append('seta pc all')
        outLines.append('swep pc new tra 0 0 {:f} {:d}'.format(self.boxHeight/2., self.calc_division(self.boxHeight/2.)))
        outLines.append('swep pc new tra 0 {:f} 0 {:d}'.format(-1*self.boxOverhang, self.calc_division(self.boxOverhang)))
        outLines.append('swep pc new tra 0 {:f} 0 {:d}'.format(self.boxDepth/2., self.calc_division((self.boxDepth/2.))))
        outLines.append('seta T1o all')
        outLines.append('')
        outLines.append('# mirro to get the TT, still flat lines')
        outLines.append('copy T1o T1u mir z a')
        outLines.append('move T1u tra 0 0 {:f}'.format(self.boxHeight))
        outLines.append('')
        outLines.append('#new set for the TT (flat lines)')
        outLines.append('seta TT T1o T1u')
        outLines.append('')
        outLines.append('# extrude the TT')
        # here a more suffisticated calculation of the division is needed, since it has to fit the n_rib
        spanDiv = self.calc_span_division(self.halfSpan)
        outLines.append('swep TT xL tra {:f} 0 0 {:d}'.format(self.halfSpan, spanDiv))
        outLines.append('seta toflip s A001 A004 A006')
        outLines.append('flip toflip')
        outLines.append('# new set for TT body')
        outLines.append('seta T1 all')

        for i in range(0, self.nRibs):
            if self.nRibs <= 1:
                spanPos = 0
            else:
                spanPos = i * (self.halfSpan / (self.nRibs-1))
            prName = 'rp{:d}'.format(i)
            ribName = 'rib{:d}'.format(i)
            outLines.append('')
            outLines.append('# generate a rib{:d}'.format(i))
            outLines.append('seta prevAll all')
            outLines.append('pnt '+prName+' {:f} {:f} 0'.format(spanPos, -1*self.boxOverhang))
            outLines.append('seta '+prName+' p '+prName+'')
            outLines.append('swep '+prName+' new tra 0 0 {:f} {:d}'.format(self.boxHeight, self.calc_division((self.boxHeight))))
            outLines.append('seta '+ribName+' l all')
            outLines.append('setr '+ribName+' l prevAll')
            outLines.append('swep '+ribName+' new tra 0 {:f} 0 {:d}'.format((self.boxDepth/2.)+self.boxOverhang, self.calc_division(((self.boxDepth/2.)+self.boxOverhang))))

        outLines.append('')
        outLines.append('# mesh it')
        outLines.append('elty all '+element_type)
        outLines.append('mesh all')
        outLines.append('')
        outLines.append('# merge TT, like welding')
        outLines.append('seta nodes n all')
        outLines.append('enq nodes tomerg rec _ _ {:f} 0.01'.format(self.boxHeight/2.))
        outLines.append('merg n tomerg')

        for i in range(0, self.nRibs):
            if self.nRibs <= 1:
                spanPos = 0
            else:
                spanPos = i * (self.halfSpan / (self.nRibs-1))
            prName = 'rp{:d}'.format(i)
            ribName = 'rib{:d}'.format(i)
            outLines.append('')
            outLines.append('# merge '+ribName)
            outLines.append('seta nodes n all')
            outLines.append('enq nodes tomerg rec {:f} _ _ 0.01'.format(spanPos))
            outLines.append('merg n tomerg')

        outLines.append('')
        outLines.append('# new set for the already connected TT with rib')
        outLines.append('seta I1 all')
        outLines.append('')
        outLines.append('# mirror TT so we get the beam')
        outLines.append('copy all I2 mir y a')
        outLines.append('move I2 tra 0 {:f} 0'.format(self.boxDepth))
        outLines.append('')
        outLines.append('# merge nodes in x and z direction')
        outLines.append('seta nodes n all')
        outLines.append('enq nodes tomerg rec _ {:f} _ 0.01'.format(self.boxDepth/2.))
        outLines.append('merg n tomerg')
        outLines.append('')
        outLines.append('# write msh file')
        outLines.append('seta nodes n all')
        outLines.append('enq nodes x0 rec 0 _ _')
        outLines.append('send x0 abq nam')
        outLines.append('send all abq')
        outLines.append('')
        #outLines.append('enq nodes load1 rec {:f} 0 40 0.1 a'.format(halfSpan))
        nodeCount = ((self.halfSpan/self.elementSize)+1) * ((((2.*self.boxOverhang)+self.boxDepth)/self.elementSize)+1)
        if element_type == 'qu8':
            nodeCount -= 0.5*(self.halfSpan/self.elementSize) * 0.5*(((2.*self.boxOverhang)+self.boxDepth)/self.elementSize)
        noadLoadTop = force_top/nodeCount
        noadLoadBut = force_but / nodeCount
        outLines.append('# load application')
        outLines.append('# top')
        outLines.append('seta nodes n all')
        outLines.append('merg n nodes')  # merge duplicate nodes
        outLines.append('enq nodes loadTop rec _ _ 0')  # .format(halfSpan))
        outLines.append('send loadTop abq force 0 0 {:f}'.format(noadLoadTop))
        outLines.append('# top')
        outLines.append('seta nodes n all')
        #outLines.append('merg n nodes')  # merge duplicate nodes
        outLines.append('enq nodes loadBut rec _ _ {:f}'.format(self.boxHeight))
        outLines.append('send loadBut abq force 0 0 {:f}'.format(noadLoadBut))
        outLines.append('')
        outLines.append('# plot it')
        outLines.append('rot -y')
        outLines.append('rot r 110')
        outLines.append('rot u 20')
        outLines.append('seta ! all')
        outLines.append('frame')
        outLines.append('view elem')
        outLines.append('plus n x0 r')
        outLines.append('plus n xL b')
        outLines.append('plus n load1 g')
        outLines.append('hcpy png')
        outLines.append('sys mv hcpy_1.png mesh.png')
        outLines.append('')
        outLines.append('quit')
        f = open(self.projectPath + '/wingGeo.fbl', 'w')
        f.writelines(line + '\n' for line in outLines)
        f.close()

    def generate_inp(self, shell_thickness):
        material_young = 69000000000.
        material_poisson = 0.32
        outLines = []
        outLines.append('** load mesh- and bc-file')
        outLines.append('*include, input=all.msh')
        outLines.append('*include, input=x0.nam')
        outLines.append('')
        outLines.append('** constraints')
        outLines.append('*boundary')
        outLines.append('Nx0,1,3')
        outLines.append('')
        outLines.append('** material definition')
        outLines.append('*MATERIAL,NAME=alu')
        outLines.append('*ELASTIC')
        outLines.append('{:f}, {:f}'.format(material_young, material_poisson))
        outLines.append('')
        outLines.append('** define surfaces')
        outLines.append('*shell section, elset=Eall, material=alu')
        outLines.append('{:f}'.format(shell_thickness))
        outLines.append('')
        outLines.append('** step')
        outLines.append('*step')
        outLines.append('*static')
        outLines.append('')
        outLines.append('** load')
        outLines.append('*cload')
        outLines.append('*include, input=loadTop.frc')
        outLines.append('*include, input=loadBut.frc')
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