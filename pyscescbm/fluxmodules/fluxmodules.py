"""
CBMPy: fluxmodules main module

This module contains convenience methods for module related computations.
For most users this is the only relevant module.
=====================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)

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

Author: Arne C. Reimers
Contact email: arne.c.reimers@gmail.com
"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

from math import isnan

import fractions

import scipy
from scipy.sparse.linalg import spsolve

from . import matroid
from . import sparserationals


def computeModules(cmod, variable=None, tol=1e-5):
    """Computes flux modules for a given model and set of variable reactions.
    
    This method only computes the minimal 0-modules with constant interface 
    flux. 
    
    This method is a convenience method for computeModulesMatroid().
    This method requires that a flux variability analysis has been performed
    beforehand. However, the set of variable reactions must not be passed
    explicitly, but it is sufficient that the variability data is encoded in
    the model.
    
    Parameters:
    - cmod a pyscescbm metabolic network model
    - variable list of reaction ids of reactions with variable flux rate. If
        this parameter is not given, the flux variability data is directly
        fetched from the model.
    - tol gives tives the tolerance for checking variability. If span <= tol, 
        then the reaction is assumed to have fixed flux. This parameter
        is only used if not a list of variable reactions is given.
    
    Returns:
    A list of lists, where each list contains the ids of
    the reactions in the module
    """ 
    
    # identify variable reactions if not explicitly given
    if variable == None:
        variable = []
        for i in range(len(cmod.reactions)):
            r = cmod.reactions[i]
            assert r.getFVAdata() != None, 'You must run FVA first'
            span = r.getFVAdata()[3]
            if (span > tol) or isnan(span):
                variable.append(r.getId())
    
    # compute modules
    modmatroids = computeModulesMatroid(cmod, variable)
    # turn into output format
    mods = []
    for m in modmatroids:
        mods.append(m.elems)
    return mods


def computeModulesMatroid(cmod, variable):
    """Computes flux modules for a given model and set of variable reactions.
    
    This method only computes the minimal 0-modules with constant interface 
    flux. Only the stoichiometry data of the model is used. 
    
    Parameters:
    - cmod a pyscescbm metabolic network model
    - variable list of reaction ids of reactions with variable flux rate
    
    Returns:
    A list of matroids, where each matroid represents a module. 
    """ 
 #   cmod = cmod.clone()
 #   for r in cmod.getReactionIds():
 #       if r not in variable:
 #           cmod.deleteReactionAndBounds(r)
 #           
 #   cmod.buildStoichMatrix(matrix_type='sympy')
#    print("finished building matrix")   
#    matrix = cmod.N.array
    matrix = sparserationals.Matrix()
    matrix.addMetabolicNetwork(cmod)
    var = []
    varNames = []
    for ri, r in enumerate(cmod.reactions):
        if r.getId() in variable:
            var.append(ri)
            varNames.append(r.getId())
    
    matrix = matrix[:,var]
        
    m = matroid.fromMatrix(matrix, varNames)
    mods = m.findModules()
    for mod in mods:
        print(mod.elems)
    
    return mods
    
def getInterface(cmod, matroid, separation):
    """ compute the interface of the separation in the given metabolic network
    
    We assume that the matroid is describing the variable reactions of cmod.
    TODO: Allow separation to also contain reactions with fixed flux rate
    and do not ignore the affine part of the interface
    """
    circuitInterface = matroid.getInterface(separation)
    
    matrix = sparserationals.Matrix()
    reactions, metabolites = matrix.addMetabolicNetwork(cmod)
    interface = []
    for c in circuitInterface:
        interfaceVector = sparserationals.Matrix()
        test = sparserationals.Matrix()
        for r,v in c.viewitems():
            rxnIdx = reactions.index(r)
            if r in separation:
                interfaceVector += matrix[:, rxnIdx] * v
            test += matrix[:, rxnIdx] * v
        #verify that it is a circuit
        for i in range(len(metabolites)):
            assert fractions.Fraction(test[i,0]) == 0
        
        # first try to simplify (compute smallest common multiple)
        scm = 1
        for i in range(len(metabolites)):
            v = fractions.Fraction(interfaceVector[i,0])
            if v != 0:
                scm *= v.denominator / fractions.gcd(scm, v.denominator)
        interfaceVector *= scm
        # now get the numerators small
        gcd = 0
        for i in range(len(metabolites)):
            v = fractions.Fraction(interfaceVector[i,0])
            if v != 0:
                gcd = fractions.gcd(gcd, v.numerator)
        interfaceVector /= gcd
        
        interfaceMap = {}
        for i in range(len(metabolites)):
            if interfaceVector[i,0] != 0:
                interfaceMap[metabolites[i]] = interfaceVector[i,0]
        interface.append(interfaceMap)
    return interface
    
    # The following is commented out, because
    # we can do this more elegantly in exact arithmetic (see code above)!
#     cmod.buildStoichMatrix('scipy_csr')
#     for c in circuitInterface:
#         # first we compute a flux vector for the circuit
#         # select one reaction and fix its flux to 1
#         # then compute flux values by solving S_c v_c = 0
#         assert len(c) > 0
#         firstRxn = c[0]
#         b = -cmod.N.array[:,cmod.N.col.index(firstRxn)]
#         rxnIdx = [cmod.N.col.index(r) for r in c[1:]]
#         A = cmod.N.array[:, rxnIdx]
#         v = spsolve(A, b)
#         sepIdx = [i for i in range(len(c)-1) if c[i+1] in separation]
#         interface = A[:,sepIdx] * v
#         if c[0] in separation:
#             interface -= b
