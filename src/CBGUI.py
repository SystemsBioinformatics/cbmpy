"""
CBMPy: CBGUI module
===================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2017 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: CBGUI.py 575 2017-04-13 12:18:44Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os

__HAVE_QT4__ = False
__HAVE_WX__ = False
try:
    from . import CBWx
    __HAVE_WX__ = True
except ImportError as ex:
    print('WARNING: Error importing wxpython module')
    #print(ex)
#try:
    #if os.sys.platform == 'win32':
        #from . import CBQt4
        #__HAVE_QT4__ = True
    #else:
        #print('\nINFO: Qt currently only available on Windows')
#except ImportError as ex:
    #print('WARNING: Error importing CBQt4 module')
    ##print(ex)

if __HAVE_WX__:
    def loadCBGUI(mod, version=2):
        """
        Load an FBA model instance into the quick editor to view or change basic model properties

         - *mod* a PySCeS CBMPy model instance

        """
        if version == 2:
            loadCBGUI = CBWx.runModelEditor(mod)
        elif version == 3:
            CBWx.runMyAUIApp()
else:
    print('\nWX GUI is not available please make sure WxPython is installed (http://www.wxpython.org)')

if __HAVE_QT4__:
    def createReaction(mod):
        """
        Load the QT4 reaction creator widget

         - *mod* a PySCeS CBMPy model instance

        """
        CBQt4.createReaction(mod)

    def openFileName(work_dir=None):
        """
        Load the QT4 file open selection dialogue

         - *work_dir* the optional initial directory

        """
        if work_dir == None:
            work_dir = os.getcwd()
        return CBQt4.fileDialogue(work_dir, 'open', None)

    def saveFileName(work_dir=None):
        """
        Load the QT4 file save selection dialogue

         - *work_dir* the optional initial directory

        """
        if work_dir == None:
            work_dir = os.getcwd()
        return CBQt4.fileDialogue(work_dir, 'save', None)
else:
    print('\nQT4 GUI tools are not available please make sure PyQT4 is installed ()')


'''
def loadBasicApp(mod=None):
    """
    Load the QT4 basic app widget

     - *mod* a PySCeS CBMPy model instance

    """
    if __HAVE_QT4__:
        CBQt4.loadBasicApp(mod)
    else:
        print('\nQT4 GUI tools are not available please make sure PyQT4 is installed ()')
'''
