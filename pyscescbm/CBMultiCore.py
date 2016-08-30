"""
CBMPy: CBMultiCore module
=========================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2015 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: CBMultiCore.py 305 2015-04-23 15:18:31Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os, time, subprocess, itertools, shutil, numpy

try:
    import pickle
except ImportError:
    import cPickle as pickle

from . import CBSolver

from . import _multicorefva
MULTIFVAFILE = __file__.replace('CBMultiCore', '_multicorefva')
del _multicorefva

try:
    from . import _multicoreenvfva
    MULTIENVFVAFILE = __file__.replace('CBMultiCore','_multicoreenvfva')
    del _multicoreenvfva
    HAVE_MULTIENV = True
except ImportError as ex:
    print(ex)
    HAVE_MULTIENV = False

from .CBConfig import __CBCONFIG__ as __CBCONFIG__

def grouper(n, iterable, padvalue=None):
    "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
    return itertools.izip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

def runMultiCoreFVA(fba, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower',\
                    optPercentage=100.0, work_dir=None, quiet=True, debug=False, oldlpgen=False, markupmodel=True, procs=2):
    """
    Run a multicore FVA where:

     - *fba* is an fba model instance
     - *procs* [default=2] number of processing threads (optimum seems to be about the number of physical cores)

    """

    #CBSolver.analyzeModel(fba, oldlpgen=False)
    #cplx_FluxVariabilityAnalysis(fba, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower', optPercentage=100.0, work_dir=None, quiet=True, debug=False, oldlpgen=False, markupmodel=True)
    fba.FVAARGS = [selected_reactions, pre_opt, tol, objF2constr, rhs_sense, optPercentage, work_dir, True, False, False, False]
    fN = str(time.time()).split('.')[0]
    fba.serializeToDisk(fN, protocol=-1)
    fN = os.path.abspath(fN)
    subprocess.call(['python', MULTIFVAFILE, str(procs), fN])
    F = file(fN, 'rb')
    res = pickle.load(F)
    F.close()
    os.remove(fN)

    fva = res[0]
    fvan = res[1]

    if len(fva) == 1:
        fva = fva[0]
        fvan = fvan[0]
    elif len(fva) > 1:
        fva = numpy.vstack(fva)
        fvan2 = []
        for n_ in fvan:
            fvan2 += n_
        fvan = fvan2

    # print 'Reaction, Reduced Costs, Variability Min, Variability Max, abs(Max-Min), MinStatus, MaxStatus'
    if markupmodel:
        for R in range(len(fvan)):
            REAC = fba.getReaction(fvan[R])
            REAC.setValue(fva[R][0])
            REAC.fva_max = fva[R][3]
            REAC.fva_min = fva[R][2]
            REAC.reduced_costs = fva[R][1]
    return fva, fvan

if HAVE_MULTIENV:
    def runMultiCoreMultiEnvFVA(lp, selected_reactions=None, tol=None, rhs_sense='lower', optPercentage=100.0, work_dir=None, debug=False, procs=2):
        """
        Run a multicore FVA where:

         - *lp* is a multienvironment lp model instance
         - *procs* [default=2] number of processing threads (optimum seems to be about the number of physical cores)

        """
        fN = os.path.join(work_dir, str(time.time()).split('.')[0])
        lp.write(fN+'.lp', filetype='lp')
        MEargs = [fN+'.lp', selected_reactions, tol, rhs_sense, optPercentage, work_dir, debug]
        print(MEargs)
        print(fN)
        F = file(fN, 'wb')
        pickle.dump(MEargs, F, protocol=-1)
        F.close()
        subprocess.call(['python', MULTIENVFVAFILE, str(procs), fN])
        F = file(fN, 'rb')
        res = pickle.load(F)
        F.close()
        os.remove(fN)

        fva = res[0]
        fvan = res[1]

        if len(fva) == 1:
            fva = fva[0]
            fvan = fvan[0]
        elif len(fva) > 1:
            fva = numpy.vstack(fva)
            fvan2 = []
            for n_ in fvan:
                fvan2 += n_
            fvan = fvan2

        return fva, fvan
else:
    def runMultiCoreMultiEnvFVA(lp, selected_reactions=None, tol=None, rhs_sense='lower', optPercentage=100.0, work_dir=None, debug=False, procs=2):
        raise RuntimeError('\nMultiCore module not present')




# look at using if __name__ != '__main__'

"""
import os, time, math, pickle
import multiprocessing
import CBSolver, CBMultiCore

def funcFVA(x,r):
    try:
        fva, fvan = CBSolver.FluxVariabilityAnalysis(x, selected_reactions=r, pre_opt=True, quiet=True, oldlpgen=False)
    except Exception as ex:
        print ex
        fva = fvan = None
        raise RuntimeError(ex)
    return fva, fvan

def multiCoreFVA(fvamod, procs=2, timeout=None):
    #assert procs >= 2, '\nMinimum 2 processes required'
    rid = fvamod.getReactionIds()
    PRs = [d_[0] for d_ in list(CBMultiCore.grouper(int(math.ceil(len(rid)/float(procs))), range(len(rid))))]
    print procs, len(PRs), PRs

    s2time = time.time()
    TP = multiprocessing.Pool(processes=procs)
    results = []
    if procs == 1:
        results.append(TP.apply_async(funcFVA, (fvamod.clone(), rid)))
    else:
        for p_ in range(len(PRs)-1):
            print '%s -> %s' % (PRs[p_], PRs[p_+1])
            results.append(TP.apply_async(funcFVA, (fvamod.clone(), rid[PRs[p_]:PRs[p_+1]])))
        results.append(TP.apply_async(funcFVA, (fvamod.clone(), rid[PRs[-1]:])))
    fva = []
    fvan = []
    for r_ in results:
        out = r_.get(timeout)
        fva.append(out[0])
        fvan.append(out[1])
    e2time = time.time()
    del results

    print '\nMulticore FVA took: %.2f min (%s processes)' % ((e2time-s2time)/60., procs)
    return fva, fvan


if __name__ == '__main__':
    print os.sys.argv

    cores = int(os.sys.argv[1])
    F = file(os.sys.argv[2],'rb')
    cmod = pickle.load(F)
    F.close()


    CBSolver.analyzeModel(cmod, oldlpgen=False)
    res = multiCoreFVA(cmod, procs=cores)

    F = file(os.sys.argv[2], 'wb')
    pickle.dump(res, F, protocol=-1)
    F.flush()
    F.close()

"""