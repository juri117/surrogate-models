#import importlib.util
#spec = importlib.util.spec_from_file_location("module.name", "/path/to/file.py")
#mdao = importlib.util.module_from_spec(spec)
#spec.loader.exec_module(mdao)

import sys
import os

sys.path.insert(0, './OpenMDAO')
sys.path.insert(0, './pyDOE2-master')

import numpy as np

from openmdao.api import Problem
from openmdao.api import MetaModelUnStructuredComp
from openmdao.api import FloatKrigingSurrogate

trig = MetaModelUnStructuredComp()
trig.add_input('x', 0.)
trig.add_output('sin_x', 0., surrogate=FloatKrigingSurrogate())
trig.add_output('cos_x', 0.)

trig.options['default_surrogate'] = FloatKrigingSurrogate()

# provide training data
trig.options['train:x'] = np.linspace(0,10,20)
trig.options['train:sin_x'] = .5*np.sin(trig.options['train:x'])
trig.options['train:cos_x'] = .5*np.cos(trig.options['train:x'])

# add it to a Problem, run and check the predicted values
prob = Problem()
prob.model.add_subsystem('trig', trig)
prob.setup(check=False)

prob['trig.x'] = 2.1
prob.run_model()
print(prob['trig.sin_x'])