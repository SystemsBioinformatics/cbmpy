"""
CBMPy: CBSandBox module
=======================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2018 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last edit: $Author: bgoli $ ($Id: CBSandBox.py 660 2018-09-24 14:57:04Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os, time, re
try:
    import pickle
except ImportError:
    import cPickle as pickle
import numpy


from .CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

# Set Suffix new function extended by Ruchir

def setSuffix_extended(model, suffix='', target='all', ignore_list=[]):
    print('\nThis function has been replaced by the new model method: cmod.renameObjectIds()\nPlease update your code!\n')
    t = 6
    while t > 0:
        print(t)
        time.sleep(1)
        t = t - 1
    if hasattr(model,'N') and hasattr(model.N,'col'):
        assert len(model.N.col)==len(model.reactions), "\nStoichiometry (number of reactions) doesn't match the model"
    if hasattr(model,'N') and hasattr(model.N,'row'):
        assert len(model.N.row)==len(model.species), "\nStoichiometry (number of species) doesn't match the model"

    assert type(suffix) == str, "\nSuffix must be a string not %s" % suffix
    assert target in ['species', 'reactions','all'], "\nTarget can only be one of ['species', 'reactions','all'] not %s" % target

    model.suffix = suffix

    if target == 'species' or target == 'all':
        for s in model.species:
            if s.getId() not in ignore_list:
                s.id = s.id+suffix
    if hasattr(model,'N'):
        if hasattr(model.N,'row'):
            rowidx = model.N.row.index(s)
            model.N.row[rowidx] = model.N.row[rowidx]+suffix

    for r in range(len(s.isReagentOf())):
        if s.reagent_of[r] not in ignore_list:
            s.reagent_of[r] = s.reagent_of[r]+suffix

    if target == 'reactions' or target == 'all':
        for r in model.reactions:
            if r.getId() not in ignore_list:
                r.id = r.id + suffix
            for reag in r.reagents:
                if reag.species_ref not in ignore_list:
                    reag.species_ref = reag.species_ref + suffix
            if hasattr(model,'N'):
                if hasattr(model.N,'col'):
                    colidx = model.N.col.index(r)
                    model.N.col[colidx] = model.N.col[colidx] + suffix
    if target == 'constraints' or target == 'reactions' or target == 'all' :
        for c in model.flux_bounds:
            if c.reaction not in ignore_list:
                c.id = c.id + suffix
                c.reaction = c.reaction + suffix

    if target == 'objectives' or target == 'all':
        for obj in model.objectives:
            obj.id = obj.id + suffix
        for f in obj.fluxObjectives:
            if f.reaction not in ignore_list:
                f.id = f.id + suffix
                f.reaction = f.reaction + suffix


