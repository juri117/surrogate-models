__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :interface to Calculix and solver-input-file generation
# author          :Juri Bieler
# date            :2017-12-13
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================


class GeometryGenerator:

    def __init__(self):
        print('done')

    def generate_wing(self, file_path, n_ribs, ribs_thickness, beam_thickness, skin_thickness):
        elemType = 'qu8'
        #elemType = 'qu4'
        #force in N at wingtip
        force = -100
        # outer geometry in m
        overhang = 10
        height = 60
        halfSpan = 400
        boxDepth = 40
        elementSize = 5
        outLines = []
        outLines.append('# draw flat T as lines')
        outLines.append('pnt pc 0 0 0')
        outLines.append('seta pc all')
        outLines.append('swep pc new tra 0 0 {:f} {:d}'.format(height/2., int((height/2.)/elementSize)))
        outLines.append('swep pc new tra 0 {:f} 0 {:d}'.format(-1*overhang, int(overhang/elementSize)))
        outLines.append('swep pc new tra 0 {:f} 0 {:d}'.format(boxDepth/2., int((boxDepth/2.)/elementSize)))
        outLines.append('seta T1o all')
        outLines.append('')
        outLines.append('# mirro to get the TT, still flat lines')
        outLines.append('copy T1o T1u mir z a')
        outLines.append('move T1u tra 0 0 {:f}'.format(height/2.))
        outLines.append('')
        outLines.append('#new set for the TT (flat lines)')
        outLines.append('seta TT T1o T1u')
        outLines.append('')
        outLines.append('# extrude the TT')
        outLines.append('swep TT xL tra {:f} 0 0 {:d}'.format(halfSpan, int(halfSpan/elementSize)))
        outLines.append('seta toflip s A001 A004 A006')
        outLines.append('flip toflip')
        outLines.append('# new set for TT body')
        outLines.append('seta T1 all')

        for i in range(0, n_ribs):
            spanPos = i * (halfSpan / (n_ribs-1))
            prName = 'rp{:d}'.format(i)
            ribName = 'rib{:d}'.format(i)
            outLines.append('')
            outLines.append('# generate a rib{:d}'.format(i))
            outLines.append('seta prevAll all')
            outLines.append('pnt '+prName+' {:f} {:f} 0'.format(spanPos, -1*overhang))
            outLines.append('seta '+prName+' p '+prName+'')
            outLines.append('swep '+prName+' new tra 0 0 {:f} {:d}'.format(height/2., int((height/2.)/elementSize)))
            outLines.append('seta '+ribName+' l all')
            outLines.append('setr '+ribName+' l prevAll')
            outLines.append('swep '+ribName+' new tra 0 {:f} 0 {:d}'.format((boxDepth/2.)+overhang, int(((boxDepth/2.)+overhang)/elementSize)))

        outLines.append('')
        outLines.append('# mesh it')
        outLines.append('elty all '+elemType)
        outLines.append('mesh all')
        outLines.append('')
        outLines.append('# merge TT, like welding')
        outLines.append('seta nodes n all')
        outLines.append('enq nodes tomerg rec 180 _ _ 0.1 #(last param is the tolerance)')
        outLines.append('merg n tomerg')

        for i in range(0, n_ribs):
            spanPos = i * (halfSpan / (n_ribs - 1))
            prName = 'rp{:d}'.format(i)
            ribName = 'rib{:d}'.format(i)
            outLines.append('')
            outLines.append('# merge '+ribName)
            outLines.append('seta nodes n all')
            outLines.append('enq nodes tomerg rec {:f} _ _ 0.1'.format(spanPos))
            outLines.append('merg n tomerg')

        outLines.append('')
        outLines.append('# new set for the already connected TT with rib')
        outLines.append('seta I1 all')
        outLines.append('')
        outLines.append('# mirror TT so we get the beam')
        outLines.append('copy all I2 mir y a')
        outLines.append('move I2 tra 0 {:f} 0'.format(boxDepth))
        outLines.append('')
        outLines.append('# merge nodes in x and z direction')
        outLines.append('seta nodes n all')
        outLines.append('enq nodes tomerg rec _ {:f} _ 0.1'.format(boxDepth/2.))
        outLines.append('merg n tomerg')
        outLines.append('')
        outLines.append('# write msh file')
        outLines.append('seta nodes n all')
        outLines.append('enq nodes x0 rec 0 _ _')
        outLines.append('send x0 abq nam')
        outLines.append('send all abq')
        outLines.append('')
        outLines.append('# load application')
        outLines.append('seta nodes n all')
        outLines.append('merg n nodes') #merge duplicate nodes
        outLines.append('enq nodes load1 rec _ _ 0')#.format(halfSpan))
        #outLines.append('enq nodes load1 rec {:f} 0 40 0.1 a'.format(halfSpan))
        nodeCount = ((halfSpan/elementSize)+1) * ((((2.*overhang)+boxDepth)/elementSize)+1)
        if elemType == 'qu8':
            nodeCount -= 0.5*(halfSpan/elementSize) * 0.5*(((2.*overhang)+boxDepth)/elementSize)
        noadLoad = force/nodeCount
        outLines.append('send load1 abq force 0 0 {:f}'.format(noadLoad))
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
        f = open(file_path + '.fbl', 'w')
        f.writelines(line + '\n' for line in outLines)
        f.close()

if __name__ == '__main__':
    geo = GeometryGenerator()
    geo.generate_wing('../dataOut/test01/test', 5, 0.1, 0.1, 0.1)