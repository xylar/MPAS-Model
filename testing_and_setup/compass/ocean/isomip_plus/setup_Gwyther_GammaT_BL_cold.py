#!/usr/bin/env python
"""
Sets up a paramter study that varies GammaT and GammaS (the
heat and salt transfer coefficients) according to the
specification of the ISOMIP+ Ocean0 experiment
"""

import numpy
import subprocess

gammaTs = numpy.array(
    [0.01, 0.02, 0.03, 0.04, 0.05, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11,
     0.12, 0.15, 0.2])
gammaSs = gammaTs/35.

blFlux = numpy.array([2., 20., 2., 20.])
blTracer = numpy.array([2., 2., 20., 20.])

GammaTs, BLFluxes = numpy.meshgrid(gammaTs, blFlux)
GammaSs, BLTracers = numpy.meshgrid(gammaSs, blTracer)

GammaTStrings = ['{:g}'.format(GammaT) for GammaT in GammaTs.ravel()]
GammaSStrings = ['{:g}'.format(GammaS) for GammaS in GammaSs.ravel()]
BLFluxStrings = ['{:g}'.format(depth) for depth in BLFluxes.ravel()]
BLTracerStrings = ['{:g}'.format(depth) for depth in BLTracers.ravel()]

args = ['../../utility_scripts/make_parameter_study_configs.py', '-t', 
        'template_Gwyther_GammaT_BL.xml', '-o', '2km/Ocean0_cold/config_GammaT_BL',
        '-p', 'GammaT={}'.format(','.join(GammaTStrings)),
        'GammaS={}'.format(','.join(GammaSStrings)), 
         'blFlux={}'.format(','.join(BLFluxStrings)),
         'blTracer={}'.format(','.join(BLTracerStrings))]
        
print('running {}'.format(' '.join(args)))
subprocess.check_call(args)
