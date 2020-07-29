"""
CBMPy: CBModelTools module
==========================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2018 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last edit: $Author: bgoli $ ($Id: CBModelTools.py 710 2020-04-27 14:22:34Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals


from .CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

import re, numpy
from . import CBModel

def addReversibilityBounds(reactions, bounds, infinity=numpy.inf):
    bounds = bounds.copy()
    for R in reactions:
        if not reactions[R]['id'] in bounds:
            if not reactions[R]['reversible']:
                bounds.update({reactions[R]['id'] : {'lower' : 0, 'upper' : infinity}})
            else:
                bounds.update({reactions[R]['id'] : {'lower' : -infinity, 'upper' : infinity}})
        else:
            print('Reaction {} already has bounds: {}'.format(reactions[R]['id'], str(reactions[R])))
    return bounds

def addReversibilityBoundsIgnoreReversible(reactions, bounds):
    bounds = bounds.copy()
    for R in reactions:
        if not reactions[R]['id'] in bounds:
            if not reactions[R]['reversible']:
                bounds.update({reactions[R]['id'] : {'lower' : 0}})
        else:
            print('Reaction {} already has bounds: {}'.format(reactions[R]['id'], str(reactions[R])))
    return bounds

def addSpecies(model, species):
    skeys = list(species)
    skeys.sort()
    for S in skeys:
        comp = 'Cell'
        ids = name = species[S].pop('id')
        if 'compartment' in species[S]:
            comp = species[S].pop('compartment')
        if 'name' in species[S]:
            name = species[S].pop('name')
        sObj = CBModel.Species(ids, boundary=species[S].pop('boundary'),\
                            name=name, value=0, compartment=comp)
        if len(species[S]) > 0:
            for a_ in species[S]:
                sObj.annotation.update({a_ : species[S][a_]})
        model.addSpecies(sObj)

def addBounds(model, bounds):
    cntr = 0
    for B in bounds:
        if 'lower' in bounds[B]:
            model.addFluxBound(CBModel.FluxBound('%s_lb' % B, B, 'greaterEqual', bounds[B]['lower']))
            cntr += 1
        if 'upper' in bounds[B]:
            model.addFluxBound(CBModel.FluxBound('%s_ub' % B, B, 'lessEqual', bounds[B]['upper']))
            cntr += 1

def addReactions(model, reactions):
    specId = model.getSpeciesIds()
    rkeys = list(reactions)
    rkeys.sort()
    for R in rkeys:
        Bounds = {}
        exchange = False
        rid = name = reactions[R]['id']
        if 'lower' in reactions[R]:
            Bounds.update({R : {'lower' : float(reactions[R]['lower'])}})
        if 'upper' in reactions[R]:
            Bounds.update({R : {'upper' : float(reactions[R]['upper'])}})
        if 'exchange' in reactions[R]:
            if reactions[R]['exchange'] in ['True','true','TRUE',True]:
                exchange = True
        if 'name' in reactions[R]:
            name = reactions[R]['name']
            ##  else:
                ##  print 'Ignoring exchange status for reaction %s, %s is not one of [True, TRUE, true]' % (R, reactions[R]['exchange'])
        if len(Bounds) > 0:
            addBounds(model, Bounds)

        revers = reactions[R]['reversible']
        ##  reagents = []
        react = CBModel.Reaction(rid, name=name, reversible=revers)
        model.addReaction(react, create_default_bounds=False)
        for RG in reactions[R]['reagents']:
            if RG[1] in specId:
                react.addReagent(CBModel.Reagent(rid+RG[1], RG[1], RG[0]))
                ##  reagents.append(CBModel.Reagent(RG[1], RG[0]))
            else:
                print(reactions[R]['reagents'])
                raise RuntimeWarning('\n%s is not a species!' % RG[1])

        # more generic way of doing this
        for k in reactions[R]:
            if k not in ['lower', 'upper', 'exchange', 'name', 'id', 'reversible', 'reagents']:
                react.annotation.update({k : reactions[R][k]})
        ##  if 'SUBSYSTEM' in reactions[R]:
            ##  react.annotation.update({'SUBSYSTEM' : reactions[R]['SUBSYSTEM']})
        ##  else:
            ##  react.annotation.update({'SUBSYSTEM' : 'metabolism'})
        react.is_exchange = exchange


def addObjectiveFunction(model, objective_function):
    for o in objective_function:
        objf = objective_function[o]
        id = objf['id']
        flux = objf['flux']
        coef = objf['coefficient']
        sense = objf['sense']
        active = objf['active']
        objF = CBModel.Objective(id, sense)
        model.addObjective(objF, active=active)
        fObj = CBModel.FluxObjective(flux+'_obj', reaction=flux, coefficient=coef)
        objF.addFluxObjective(fObj)


def quickDefaultBuild(model_name, Reactions, Species, Bounds, Objective_function, infinity=numpy.inf):
    fba = CBModel.Model(model_name)
    addSpecies(fba, Species)
    addReactions(fba, Reactions)
    addObjectiveFunction(fba, Objective_function)
    Bounds = addReversibilityBounds(Reactions, Bounds, infinity)
    addBounds(fba, Bounds)
    return fba

def printSolution(fba, wait=False):
    for R in fba.reactions:
        print(R.getId(), R.value, end=" ")
        if 'SUBSYSTEM' in  R.annotation:
            print(R.annotation['SUBSYSTEM'])
        else:
            print(' ')
    if wait: raw_input('Press <enter> to continue ...')
