"""
CBMPy: MultiCoreFVA module
==========================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2024 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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

Author: Brett G. Olivier PhD
Contact developers: https://github.com/SystemsBioinformatics/cbmpy/issues
Last edit: $Author: bgoli $ ($Id: _multicorefva.py 710 2020-04-27 14:22:34Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# from __future__ import unicode_literals

import os
import time
import math
import copy
import multiprocessing

try:
    import pickle
except ImportError:
    import cPickle as pickle

import cbmpy


def funcFVA(x, r):
    try:
        Args = x.FVAARGS
        del x.FVAARGS
        Args.insert(0, x)
        Args[1] = r
        fva, fvan = cbmpy.CBSolver.FluxVariabilityAnalysis(*Args)
    except Exception as ex:
        print(ex)
        fva = fvan = None
        raise RuntimeError(ex)
    return fva, fvan


def multiCoreFVA(fvamod, procs=2, timeout=None):
    # assert procs >= 2, '\nMinimum 2 processes required'
    rid = fvamod.getReactionIds()
    PRs = [
        d_[0]
        for d_ in list(
            cbmpy.CBMultiCore.grouper(
                int(math.ceil(len(rid) / float(procs))), range(len(rid))
            )
        )
    ]
    print(procs, len(PRs), PRs)

    s2time = time.time()
    TP = multiprocessing.Pool(processes=procs)
    results = []
    if procs == 1:
        results.append(TP.apply_async(funcFVA, (fvamod.clone(), rid)))
    else:
        for p_ in range(len(PRs) - 1):
            print('{} -> {}'.format(PRs[p_], PRs[p_ + 1]))
            print('{} -> {}'.format(rid[PRs[p_]], rid[PRs[p_ + 1]]))
            # investigate using without clone
            # results.append(TP.apply_async(funcFVA, (fvamod, rid[PRs[p_]:PRs[p_+1]])))
            results.append(
                TP.apply_async(funcFVA, (fvamod.clone(), rid[PRs[p_] : PRs[p_ + 1]]))
            )
        results.append(TP.apply_async(funcFVA, (fvamod.clone(), rid[PRs[-1] :])))
    fva = []
    fvan = []
    for r_ in results:
        out = r_.get(timeout)
        fva.append(out[0])
        fvan.append(out[1])
    e2time = time.time()
    del results

    print(
        '\nMulticore FVA took: {} min ({} processes)'.format(
            (e2time - s2time) / 60.0, procs
        )
    )
    return fva, fvan


if __name__ == '__main__':

    print(os.sys.argv)

    cores = int(os.sys.argv[1])
    F = open(os.sys.argv[2], 'rb')
    cmod = pickle.load(F)
    F.close()

    cbmpy.CBSolver.analyzeModel(cmod, oldlpgen=False)
    res = multiCoreFVA(cmod, procs=cores)

    F = open(os.sys.argv[2], 'wb')
    pickle.dump(res, F, protocol=-1)
    F.flush()
    F.close()
