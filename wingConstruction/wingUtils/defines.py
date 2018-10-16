__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :definition of constants
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

#######################################################
# PARAMETERS

max_g = 1. #2.5
safety_fac = 1.5 # 1.5
mtow = 27987.
fuel_mass_in_wings = 2 * 2659.
first_wing_struct_mass = 2 * 1000.
wing_load = ((mtow - fuel_mass_in_wings - first_wing_struct_mass) * 9.81) * max_g * 0.5
engine_weight = (873.1 + 251.9 + 112.7 + 62.8) * 9.81  # engine, prop, wheel, brake
engine_pos_y = 3.
wing_length = 12.87
chord_length = 3.
chord_height = 0.55
shell_thickness = 0.005

max_shear_strength = 5.72e8 / safety_fac  # 3.31e8 #Pa

density = 2810 #kg/m^3

#shear_strength = 3.31e8 / safety_fac

#range_rib = (5, 18) #calcu
range_rib = (5, 25) #range_rib = (8, 24) #aba
range_shell = (0.002, 0.0033)

#######################################################
# CONSTANTS FOR INDEXING AND CHOICE

SAMPLE_NAMES = ['LatinHyperCube', 'Halton', 'StructuredSampling', 'OptimizedLatinHyperCubePyKriging']
SAMPLE_LATIN = 0
SAMPLE_HALTON = 1
SAMPLE_STRUCTURE = 2
SAMPLE_OPTI_LATIN_HYPER = 3

SURRO_NAMES = ['Kriging', 'RBF', 'Polynom', 'PyKriging']
SURRO_KRIGING = 0
SURRO_RBF = 1
SURRO_POLYNOM = 2
SURRO_PYKRIGING = 3