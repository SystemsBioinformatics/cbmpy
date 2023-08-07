"""
CBMPy: CBMultiModel module
==========================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2022 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: CBTools.py 544 2022-01-12 16:31:50Z bgoli $)

"""
# gets rid of "invalid variable name" info
# pylint: disable=C0103
# gets rid of "line to long" info
# pylint: disable=C0301
# use with caution: gets rid of module xxx has no member errors (run once enabled)
# pylint: disable=E1101

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# from __future__ import unicode_literals

from .CBConfig import __CBCONFIG__ as __CBCONFIG__

__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']


def copySpecies(m_src, m_targ, sid, is_boundary=False):
    """
    Copy a species from a source model to a target model

     - *m_src* the source model
     - *m_targ* the target model
     - *sid* the species to copy
     - *is_boundary* [default] add the species as a variable metabolite or boundary if True

    """
    out = {}
    if m_targ.getSpecies(sid) is not None:
        print('ERROR: species with id \"{}\" exists in target model'.format(sid))
        out = None
    if m_src.getSpecies(sid) is None:
        print(
            'ERROR: species with id \"{}\" does not exist in source model'.format(sid)
        )
        out = None
    if out is None:
        return None

    S = m_src.getSpecies(sid).clone()
    if is_boundary:
        S.setBoundary()
    else:
        S.unsetBoundary()
    m_targ.addSpecies(S)
    del S
    return True


def copyReaction(m_src, m_targ, rid, altrid=None):
    """
    Copy a reaction from a source model to a target model, if the required species exist in the target
    then they are mapped as reagents, otherwise new metabolites are added as boundary species.

     - *m_src* the source model
     - *m_targ* the target model
     - *rid* the reaction id to copy
     - *altrid* if the reaction name exists in the target, try use this one instead

    """
    out = {}
    targ_exists = False
    # print('1', m_src, m_targ, rid, altrid)
    # print('2', m_targ.getReaction(rid), altrid, m_targ.getReaction(altrid))
    # print('3', m_targ.getReaction(rid) is not None and altrid is not None and m_targ.getReaction(altrid) is not None)
    if m_src.getReaction(rid) is None:
        print(
            'ERROR: reaction with id \"{}\" does not exist in source model'.format(rid)
        )
        out = None
    if out is not None and m_targ.getReaction(rid) is not None and altrid is not None:
        print(
            'INFO: reaction with id \"{}\" exists in target model trying alternate id'.format(
                rid
            )
        )
        if m_targ.getReaction(altrid) is not None:
            print(
                'ERROR: alternative reaction with id \"{}\" already exists in target model'.format(
                    rid
                )
            )
            out = None
        else:
            targ_exists = True
    if out is None:
        return None

    R = m_src.getReaction(rid).clone()
    if targ_exists and altrid is not None:
        R.setId(altrid)
        for re in R.reagents:
            re.setId('{}_{}'.format(altrid, re.getSpecies()))
    tSpecies = m_targ.getSpeciesIds()
    out['new_species'] = []
    out['existing_species'] = []
    out['reagents'] = []
    out['unmapped_species_reactions'] = {}
    for s in R.getSpeciesIds():
        if s not in tSpecies:
            S = m_src.getSpecies(s).clone()
            S.setBoundary()
            m_targ.addSpecies(S)
            out['new_species'].append(s)
            reag = m_src.getSpecies(s).isReagentOf()
            reag.remove(rid)
            reag.sort()
            out['unmapped_species_reactions'][s] = reag
        else:
            out['existing_species'].append(s)
        out['reagents'].append(s)
    m_targ.addReaction(R, create_default_bounds=True, silent=True)
    del R

    return out
