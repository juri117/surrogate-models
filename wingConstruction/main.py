

from wingConstruction.WingConstruction import WingConstruction
from wingConstruction.Calculix import Calculix

geo = WingConstruction()
geo.generate_wing('../dataOut/test01/test', 18, 0.1, 0.1, 0.1)
geo.generate_inp('../dataOut/test01/test', .002)

clx = Calculix(workingDir='../dataOut/test01')
clx.generate_mesh('test')

clx.solve_model('test')