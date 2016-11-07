"""
CBMPy: MultiCoreEnvFVA module
=============================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2016 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: _multicoreenvfva.py 515 2016-11-07 14:20:11Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals


import os, time, math
try:
    import pickle
except ImportError:
    import cPickle as pickle
import multiprocessing
import cbmpy

def funcFVA(x, r, tol, rhs_sense, optPercentage, work_dir, debug):
    x = cbmpy.CBCPLEX.cplx_getModelFromLP(x)
    x.parameters.threads.set(4)
    try:

        fva, fvan = cbmpy.CBCPLEX.cplx_MultiFluxVariabilityAnalysis(x, r, tol, rhs_sense, optPercentage, work_dir, debug)
    except Exception as ex:
        print(ex)
        fva = fvan = None
        raise RuntimeError(ex)
    return fva, fvan

import copy
def multiCoreEnvFVA(MEargs, procs=2, timeout=None):
    #assert procs >= 2, '\nMinimum 2 processes required'
    rid = MEargs[1]
    work_dir = MEargs[5]
    PRs = [d_[0] for d_ in list(cbmpy.CBMultiCore.grouper(int(math.ceil(len(rid)/float(procs))), range(len(rid))))]
    print(procs, len(PRs), PRs)

    s2time = time.time()
    TP = multiprocessing.Pool(processes=procs)
    results = []
    if procs == 1:
        MEargs2 = copy.deepcopy(MEargs)
        results.append(TP.apply_async(funcFVA, MEargs2))
    else:
        for p_ in range(len(PRs)-1):
            print('{} -> {}'.format(PRs[p_], PRs[p_+1]))
            print('{} -> {}'.format(rid[PRs[p_]], rid[PRs[p_+1]]))
            print(rid[PRs[p_]:PRs[p_+1]])
            MEargs2 = copy.deepcopy(MEargs)
            MEargs2[1] = rid[PRs[p_]:PRs[p_+1]]
            results.append(TP.apply_async(funcFVA, MEargs2))
        MEargs2 = copy.deepcopy(MEargs)
        MEargs2[1] = rid[PRs[-1]:]
        print(rid[PRs[-1]:])
        results.append(TP.apply_async(funcFVA, MEargs2))
    fva = []
    fvan = []
    for r_ in results:
        out = r_.get(timeout)
        fva.append(out[0])
        fvan.append(out[1])
    e2time = time.time()
    del results

    print('\nMulticore FVA took: {}.2f min ({} processes)'.format((e2time-s2time)/60., procs))
    return fva, fvan


if __name__ == '__main__':
    print(os.sys.argv)

    cores = int(os.sys.argv[1])
    F = file(os.sys.argv[2],'rb')
    MEargs = pickle.load(F)
    print('\t', MEargs)
    F.close()

    res = multiCoreEnvFVA(MEargs, procs=cores)

    F = file(os.sys.argv[2], 'wb')
    pickle.dump(res, F, protocol=-1)
    F.flush()
    F.close()
    os.remove(MEargs[0])
