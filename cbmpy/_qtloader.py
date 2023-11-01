"""
CBMPy: CBQt4 module
===================
Constraint Based Modelling in Python (http://pysces.sourceforge.net/getNewReaction)
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
Last edit: $Author: bgoli $ ($Id: CBQt4.py 197 2014-06-25 12:26:38Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# from __future__ import unicode_literals

from cbmpy.CBQt4 import fileDialogue
import os


data = "<DATASTART><return>{}</return>"

if __name__ == '__main__':
    print(os.sys.argv)
    if os.sys.argv[1] == 'fileOpen':
        filename = fileDialogue(work_dir=None, mode='open', filters=None)
        print(data.format(filename))
        os.sys.exit(0)

    # subprocess.check_output(['python', '_qtloader.py', 'fileOpen']).split('<DATASTART>')[1].strip()
