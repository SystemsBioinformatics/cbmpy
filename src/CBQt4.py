"""
CBMPy: CBQt4 module
===================
Constraint Based Modelling in Python (http://pysces.sourceforge.net/getNewReaction)
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
Last edit: $Author: bgoli $ ($Id: CBQt4.py 515 2016-11-07 14:20:11Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os, time, random, math, re, webbrowser, urllib2
from .CBCommon import pp_chemicalFormula
from .CBModel import Reaction as CBMReaction


HAVE_QT4 = False
try:
    import PyQt4
    from PyQt4 import QtCore, QtGui, QtSvg
    from PyQt4.QtWebKit import QGraphicsWebView
    HAVE_QT4 = True
    print('Qt4 GUI tools available')
except ImportError as ex:
    print('\nQt4 GUI tools not available.')
    print(ex)


class ReactionCreator(QtGui.QWidget):

    #_fba = None
    _mlist = None
    _rlist = None
    _flist = None
    _newSubs = None
    _newProds = None
    _fixColour = None
    _errColour = None
    _goodColour = None
    _cfdict = None
    _cndict = None
    _ccdict = None
    IGNORECHECK = False
    ISBALANCED = True
    NewReaction = None
    _Blower = '-inf'
    _Bupper = 'inf'

    def __init__(self, rlist, mlist, flist, cfdict, cndict, ccdict):
        super(ReactionCreator, self).__init__()
        self.mousePos = self.cursor().pos()
        rlist.sort()
        mlist.sort()
        flist.sort()
        self._rlist = rlist
        self._mlist = mlist
        self._flist = flist
        self._cfdict = cfdict
        self._cndict = cndict
        self._ccdict = ccdict
        self._newSubs = []
        self._newProds = []
        self._fixColour = QtGui.QColor(0,0,153,alpha=255)
        self._errColour = QtGui.QColor(255,0,0,alpha=255)
        self._goodColour = QtGui.QColor(0,100,0,alpha=255)
        self.initUI()

    def addSubstrate(self, coeff, sid):
        self.tblSub.insertRow(self.tblSubRow)
        self.tblSub.setItem(self.tblSubRow, 0, QtGui.QTableWidgetItem('{}'.format(coeff)))
        self.tblSub.setItem(self.tblSubRow, 1, QtGui.QTableWidgetItem('{}'.format(sid)))
        CF = 'None'
        if sid in self._cfdict:
            CF = self._cfdict[sid]
        self.tblSub.item(self.tblSubRow, 1).setToolTip(CF)
        if sid in self._flist:
            self.tblSub.item(self.tblSubRow, 1).setForeground(self._fixColour)
        self.tblSub.item(self.tblSubRow, 0).setTextAlignment(QtCore.Qt.AlignCenter)
        self.tblSubRow += 1

    def addSelectedSubstrates(self):
        self.IGNORECHECK = True
        items = [str(it_.text()) for it_ in self.lstSub.selectedItems()]
        self.setFocus(PyQt4.QtCore.Qt.OtherFocusReason)
        #print(items)
        for i_ in items:
            if i_ not in self._newSubs:
                self.addSubstrate(1, i_)
                self._newSubs.append(i_)
        self.IGNORECHECK = False
        self.statusBar.showMessage('Substrates(s) added')
        self.checkBalance()

    def addProduct(self, coeff, sid):
        self.tblProd.insertRow(self.tblProdRow)
        self.tblProd.setItem(self.tblProdRow, 0, QtGui.QTableWidgetItem('{}'.format(coeff)))
        self.tblProd.setItem(self.tblProdRow, 1, QtGui.QTableWidgetItem('{}'.format(sid)))
        CF = 'None'
        if sid in self._cfdict:
            CF = self._cfdict[sid]
        self.tblProd.item(self.tblProdRow, 1).setToolTip(CF)
        if sid in self._flist:
            self.tblProd.item(self.tblProdRow, 1).setForeground(self._fixColour)
        self.tblProd.item(self.tblProdRow, 0).setTextAlignment(QtCore.Qt.AlignCenter)
        self.tblProdRow += 1

    def addSelectedProducts(self):
        self.IGNORECHECK = True
        items = [str(it_.text()) for it_ in self.lstProd.selectedItems()]
        self.setFocus(PyQt4.QtCore.Qt.OtherFocusReason)
        #print(items)
        for i_ in items:
            if i_ not in self._newProds:
                self.addProduct(1, i_)
                self._newProds.append(i_)
        self.IGNORECHECK = False
        self.statusBar.showMessage('Product(s) added')
        self.checkBalance()

    #def keyPressEvent(self, event):
        #print('KeyPress key: {}'.format(str(event.key())))
        #if event.key() == 16777223:
            #print('You pressed the delete key')

    def deleteSubstrates(self):
        self.deleteReagents('substrate')

    def deleteProducts(self):
        self.deleteReagents('product')

    def deleteAllSubstrates(self):
        self.tblSub.clear()
        for r_ in range(self.tblSubRow-1,-1,-1):
            self.tblSub.removeRow(r_)
        self.tblSubRow = 0
        self.checkBalance()

    def deleteAllProducts(self):
        self.tblProd.clear()
        for r_ in range(self.tblProdRow-1,-1,-1):
            self.tblProd.removeRow(r_)
        self.tblProdRow = 0
        self.checkBalance()

    def deleteReagents(self, reagentType):
        selected = None
        PRODACTIVE = False
        SUBACTIVE = False
        self.IGNORECHECK = True
        if reagentType == 'substrate':
            selected = [(it_.row(), it_.column()) for it_ in self.tblSub.selectedItems()]
            SUBACTIVE = True
        elif reagentType == 'product':
            selected = [(it_.row(), it_.column()) for it_ in self.tblProd.selectedItems()]
            PRODACTIVE = True

        if selected != None:
            deleteRow = []
            if len(selected) == 2:
                if selected[0][0] == selected[1][0]:
                    if selected[0][1]+1 == selected[1][1]:
                        deleteRow.append(selected[0][0])
            elif len(selected) > 2:
                for it_ in range(0, len(selected), 2):
                    if selected[it_][1] == selected[it_+1][1]:
                        if selected[it_][0]+1 == selected[it_+1][0]:
                            if selected[it_][0] not in deleteRow:
                                deleteRow.append(selected[it_][0])
                            if selected[it_][0]+1 not in deleteRow:
                                deleteRow.append(selected[it_][0]+1)
            deleteRow.sort()
            for d_ in range(len(deleteRow)-1,-1,-1):
                if SUBACTIVE:
                    print('Deleting Sub table row: {}'.format(deleteRow[d_]))
                    self.statusBar.showMessage('Substrate(s) deleted')
                    #print(self._newSubs)
                    #print(str(self.tblSub.item(d_, 1).text()))
                    self._newSubs.pop(self._newSubs.index(str(self.tblSub.item(d_, 1).text())))
                    self.tblSub.removeRow(deleteRow[d_])
                    self.tblSubRow -= 1
                elif PRODACTIVE:
                    print('Deleting Prod table row: {}'.format(deleteRow[d_]))
                    self.statusBar.showMessage('Product(s) deleted')
                    #print(self._newProds)
                    #print(str(self.tblProd.item(d_, 1).text()))
                    self._newProds.pop(self._newProds.index(str(self.tblProd.item(d_, 1).text())))
                    self.tblProd.removeRow(deleteRow[d_])
                    self.tblProdRow -= 1
        self.IGNORECHECK = False
        self.checkBalance()

    def checkBalance(self):
        if self.IGNORECHECK:
            return
        output = {}
        left = {}
        right = {}
        for r_ in range(self.tblSubRow):
            sid = str(self.tblSub.item(r_, 1).text())
            scoef = float(str(self.tblSub.item(r_, 0).text()))
            #print scoef
            if sid in self._cfdict:
                cf = self._cfdict[sid]
                if cf not in [None, 'None', '', ' ']:
                    cfl = pp_chemicalFormula.parseString(cf).asList()
                else:
                    cfl = []
                #print sid, cf, cfl
                for e_ in cfl:
                    if e_[0] in output:
                        output[e_[0]] = output[e_[0]] + -scoef*float(e_[1])
                        #print scoef*float(e_[1])
                    else:
                        output[e_[0]] = -scoef*float(e_[1])
                        #print scoef*float(e_[1])
                    if e_[0] in left:
                        left[e_[0]] = left[e_[0]] + scoef*float(e_[1])
                    else:
                        left[e_[0]] = scoef*float(e_[1])
        for r_ in range(self.tblProdRow):
            sid = str(self.tblProd.item(r_, 1).text())
            pcoef = float(str(self.tblProd.item(r_, 0).text()))
            #print pcoef, type(pcoef)
            if sid in self._cfdict:
                cf = self._cfdict[sid]
                if cf not in [None, 'None', '', ' ']:
                    cfl = pp_chemicalFormula.parseString(cf).asList()
                else:
                    cfl = []
                #print sid, cf, cfl
                for e_ in cfl:
                    if e_[0] in output:
                        #print -pcoef*float(e_[1])
                        output[e_[0]] = output[e_[0]] + pcoef*float(e_[1])
                    else:
                        #print -pcoef*float(e_[1])
                        output[e_[0]] = pcoef*float(e_[1])
                    if e_[0] in right:
                        #print -pcoef*float(e_[1])
                        right[e_[0]] = right[e_[0]] + pcoef*float(e_[1])
                    else:
                        #print -pcoef*float(e_[1])
                        right[e_[0]] = pcoef*float(e_[1])
        #print output
        self.updateBalance(output, left, right)
        #self.txtBal.setText(str(output))

    def updateBalance(self, bdict, left, right):
        colHead = []
        keys = list(bdict)
        if self.tblBalCol > 0:
            for c_ in range(self.tblBalCol-1,-1,-1):
                self.tblBal.removeColumn(c_)
        self.tblBalCol = 0
        self.ISBALANCED = True
        for k_ in range(len(keys)):
            if not keys[k_] in left:
                left[keys[k_]] = 'None'
            if not keys[k_] in right:
                right[keys[k_]] = 'None'
            self.tblBalCol += 1
            self.tblBal.insertColumn(k_)
            self.tblBal.setItem(0, k_, QtGui.QTableWidgetItem('{}'.format(abs(bdict[keys[k_]]))))
            self.tblBal.setItem(1, k_, QtGui.QTableWidgetItem('{}'.format(left[keys[k_]])))
            self.tblBal.setItem(2, k_, QtGui.QTableWidgetItem('{}'.format(right[keys[k_]])))
            self.tblBal.item(0, k_).setTextAlignment(QtCore.Qt.AlignCenter)
            if k_ == 0:
                boldFont = self.tblBal.item(0, k_).font()
                boldFont.setBold(True)
            self.tblBal.item(0, k_).setFont(boldFont)
            self.tblBal.item(1, k_).setTextAlignment(QtCore.Qt.AlignCenter)
            self.tblBal.item(2, k_).setTextAlignment(QtCore.Qt.AlignCenter)
            if bdict[keys[k_]] != 0.0:
                self.ISBALANCED = False
                self.tblBal.item(0, k_).setForeground(self._errColour)
            else:
                self.tblBal.item(0, k_).setForeground(self._goodColour)
        self.tblBal.setHorizontalHeaderLabels(QtCore.QStringList(keys))


    def showErrorMessage(self, errorMsg, title="Reaction Creator"):
        QtGui.QMessageBox.critical(None, title,\
                                   errorMsg,\
                                   QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default,\
                                   QtGui.QMessageBox.NoButton)
        self.statusBar.showMessage(errorMsg)

    def getNewReaction(self):
        reversible = self.btReverse.isChecked()
        if not reversible and self._Blower == '-inf':
            self._Blower = 0.0
        Id = str(self.txtId.toPlainText()).strip()
        Name = str(self.txtNm.toPlainText()).strip()

        if Id == 'NewReactionId':
            print('\nWARNING: using default reaction id')
        errorMsg = None
        if Id == "":
            errorMsg = 'Reaction ID must be specified.'
        elif Id in self._rlist:
            errorMsg = 'Reaction ID \"{}\" already exists.'.format(Id)
        if errorMsg != None:
            self.showErrorMessage(errorMsg)
            self.NewReaction = None
            return None
        if self.tblSubRow == 0 and self.tblProdRow == 0:
            self.showErrorMessage('At least one reagent must be defined.')
            self.NewReaction = None
            return None

        #print('\nid=\"{}\"\nname=\"{}\"'.format(Id, Name))
        #print self.tblSub.rowCount(), self.tblSubRow
        #print self.tblProd.rowCount(), self.tblProdRow
        Reag = {}
        exReac = False
        for s_ in range(self.tblSub.rowCount()):
            coeff = -abs(float(self.tblSub.item(s_, 0).text()))
            Sid = str(self.tblSub.item(s_, 1).text()).strip()
            if Sid in self._flist:
                exReac = True
            Reag[Sid] = coeff
        for p_ in range(self.tblProd.rowCount()):
            coeff = abs(float(self.tblProd.item(p_, 0).text()))
            Sid = str(self.tblProd.item(p_, 1).text()).strip()
            if Sid in self._flist:
                exReac = True
            if Sid in Reag:
                Reag[Sid] += coeff
            else:
                Reag[Sid] = coeff

        for r_ in tuple(Reag):
            if Reag[r_] == 0.0:
                Reag.pop(r_)
                print('removing zero coefficient reagent: {}'.format(r_))
        self.NewReaction = {'reversible' : reversible,
                           'id' : Id,
                           'name' : Name,
                           'is_exchange' : exReac,
                           'is_balanced'  : self.ISBALANCED,
                           'reagents' : Reag,
                           'upper_bound' : self._Bupper,
                           'lower_bound' : self._Blower
                           }


        sub = ''
        prod = ''
        for r_ in Reag:
            coeff = abs(Reag[r_])
            if Reag[r_] < 0.0:
                if coeff == 1.0:
                    sub += '%s + ' % (r_)
                else:
                    sub += '{%s} %s + ' % (coeff, r_)
            else:
                if coeff == 1.0:
                    prod += '%s + ' % (r_)
                else:
                    prod += '{%s} %s + ' % (coeff, r_)
        if reversible:
            eq = '%s\n\t%s\n%s' % (sub[:-3], '<==>', prod[:-2])
        else:
            eq = '%s\n\t%s\n%s' % (sub[:-3], '-->', prod[:-2])

        #quit_msg = "Add reaction:\n\n{}\n\t{}\n{}".format(sub,rev,prod)
        quit_msg = eq
        reply = QtGui.QMessageBox.question(self, 'Do you want to add the reaction \"{}\" to the model?'.format(Id),
                                           quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            self.statusBar.showMessage('Reaction {} added to model'.format(Id))
            print(self.NewReaction)
            QtGui.qApp.quit()
        else:
            print('Try again')
            self.NewReaction = None

    def initUI(self):
        # create labels
        lblSub = QtGui.QLabel('Substrates')
        lblSub.setAlignment(QtCore.Qt.AlignCenter)
        lblProd = QtGui.QLabel('Products')
        lblProd.setAlignment(QtCore.Qt.AlignCenter)
        lblId = QtGui.QLabel('Reaction id')
        lblId.setAlignment(QtCore.Qt.AlignCenter)
        lblNm = QtGui.QLabel('Reaction name')
        lblNm.setAlignment(QtCore.Qt.AlignCenter)

        # create text boxes
        self.txtId = QtGui.QTextEdit()
        self.txtId.setMaximumHeight(25)
        self.txtId.setText('NewReactionId')
        self.txtNm = QtGui.QTextEdit()
        self.txtNm.setMaximumHeight(25)
        self.txtNm.setText('NewReactionName')
        self.txtBal = QtGui.QTextEdit()
        self.txtBal.setMaximumHeight(40)

        # create static lists
        self.lstSub = QtGui.QListWidget()
        self.lstSub.setSelectionMode(self.lstSub.ExtendedSelection)
        self.lstProd = QtGui.QListWidget()
        self.lstProd.setSelectionMode(self.lstProd.ExtendedSelection)
        # populate lists
        cntr = 0
        for m_ in self._mlist:
            cntr += 1
            name = 'None'
            comp = 'None'
            if m_ in self._ccdict:
                comp = self._ccdict[m_]
            if m_ in self._cndict:
                name = self._cndict[m_]
            item = QtGui.QListWidgetItem(m_)
            item.setToolTip('{}\t{}'.format(name, comp))
            if m_ in self._flist:
                item.setForeground(self._fixColour)
            self.lstSub.addItem(item.clone())
            self.lstProd.addItem(item)
            #if cntr == 20: break

        # create buttons
        self.btAddSub = QtGui.QPushButton('Add substrate(s)')
        QtCore.QObject.connect(self.btAddSub, QtCore.SIGNAL('clicked()'), self.addSelectedSubstrates)
        self.btAddProd = QtGui.QPushButton('Add product(s)')
        QtCore.QObject.connect(self.btAddProd, QtCore.SIGNAL('clicked()'), self.addSelectedProducts)
        self.btReverse = QtGui.QPushButton('Reversible')
        self.btReverse.setCheckable(True)
        self.btReverse.setChecked(True)

        # create tables
        self.tblSub = QtGui.QTableWidget()
        self.tblSub.setSortingEnabled(True)
        self.tblSub.insertColumn(0)
        self.tblSub.insertColumn(1)
        self.tblSub.setHorizontalHeaderLabels(QtCore.QStringList(('Coefficient','Metabolite')))
        self.tblSub.verticalHeader().setVisible(False)
        QtCore.QObject.connect(self.tblSub, QtCore.SIGNAL('cellChanged(int,int)'), self.checkBalance)
        self.tblSubRow = 0
        self.tblProd = QtGui.QTableWidget()
        self.tblProd.setSortingEnabled(True)
        self.tblProd.insertColumn(0)
        self.tblProd.insertColumn(1)
        self.tblProd.setHorizontalHeaderLabels(QtCore.QStringList(('Coefficient','Metabolite')))
        self.tblProd.verticalHeader().setVisible(False)
        self.tblProdRow = 0
        QtCore.QObject.connect(self.tblProd, QtCore.SIGNAL('cellChanged(int,int)'), self.checkBalance)

        self.tblBal  = QtGui.QTableWidget()
        self.tblBal.setMaximumHeight(150)
        self.tblBal.insertRow(0)
        self.tblBal.insertRow(1)
        self.tblBal.insertRow(2)
        self.tblBal.verticalHeader().setVisible(False)
        self.tblBalCol = 0


        # set up menu and status bar
        menuBar = QtGui.QMenuBar()

        exitAction = QtGui.QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit and loose changes')
        exitAction.triggered.connect(QtGui.qApp.quit)
        addReAction = QtGui.QAction('&Add reaction and exit', self)
        addReAction.setShortcut('Ctrl+A')
        addReAction.setStatusTip('Add reaction to model and exit')
        addReAction.triggered.connect(self.getNewReaction)

        fileMenu = menuBar.addMenu('&Model')
        fileMenu.addAction(exitAction)
        fileMenu.addAction(addReAction)

        subAddAction = QtGui.QAction('&Add', self)
        subAddAction.triggered.connect(self.addSelectedSubstrates)
        subDelAction = QtGui.QAction('&Delete selected', self)
        subDelAction.triggered.connect(self.deleteSubstrates)
        subDelAllAction = QtGui.QAction('&Delete all', self)
        subDelAllAction.triggered.connect(self.deleteAllSubstrates)
        subMenu = menuBar.addMenu('&Substrate')
        subMenu.addAction(subAddAction)
        subMenu.addAction(subDelAction)
        subMenu.addAction(subDelAllAction)

        prodAddAction = QtGui.QAction('&Add', self)
        prodAddAction.triggered.connect(self.addSelectedProducts)
        prodDelAction = QtGui.QAction('&Delete selected', self)
        prodDelAction.triggered.connect(self.deleteProducts)
        prodDelAllAction = QtGui.QAction('&Delete all', self)
        prodDelAllAction.triggered.connect(self.deleteAllProducts)
        prodMenu = menuBar.addMenu('&Product')
        prodMenu.addAction(prodAddAction)
        prodMenu.addAction(prodDelAction)
        prodMenu.addAction(prodDelAllAction)

        self.statusBar = QtGui.QStatusBar()
        self.statusBar.showMessage('{} ready'.format('Model'))

        # do layout
        grid = QtGui.QGridLayout()
        grid.setSpacing(10)

        grid.addWidget(menuBar, 0, 0, 1, 2)
        grid.addWidget(self.statusBar, 0, 2, 1, 2)
        grid.addWidget(lblId, 1, 0)
        grid.addWidget(self.txtId, 1, 1)
        grid.addWidget(lblNm, 1, 2)
        grid.addWidget(self.txtNm, 1, 3)
        grid.addWidget(lblSub, 2, 0, 1, 2)
        grid.addWidget(lblProd, 2, 2, 1, 2)
        grid.addWidget(self.lstSub, 3, 0, 1, 2)
        grid.addWidget(self.lstProd, 3, 2, 1, 2)
        grid.addWidget(self.btAddSub, 4, 0)
        grid.addWidget(self.btReverse, 4, 1, 1, 2)
        grid.addWidget(self.btAddProd, 4, 3)
        grid.addWidget(self.tblSub, 5, 0, 1, 2)
        grid.addWidget(self.tblProd, 5, 2, 1, 2)
        grid.addWidget(self.tblBal, 6, 0, 1, 4)

        self.setLayout(grid)
        self.setGeometry(self.mousePos.x()-75, self.mousePos.y()-75, 500, 640)
        self.setWindowTitle('Reaction Creator')
        self.show()

def createReaction(mod):
    """
    Create a reaction using the graphical Reaction Creator

     - *mod* a CBMPy model object

    """
    cfdict = {}
    cndict = {}
    ccdict = {}
    for s_ in mod.species:
        cfdict[s_.getPid()] = s_.getChemFormula()
        cndict[s_.getPid()] = s_.getName()
        ccdict[s_.getPid()] = s_.compartment
    app = QtGui.QApplication([])
    ex = ReactionCreator(mod.getReactionIds(), mod.getSpeciesIds(), mod.getBoundarySpeciesIds(), cfdict, cndict, ccdict)
    app.exec_()
    newR = ex.NewReaction
    del app, ex
    if newR == None:
        return None
    else:
        R = CBMReaction(newR['id'], name=newR['name'], reversible=newR['reversible'])
        R.is_balanced = newR['is_balanced']
        R.is_exchange = newR['is_exchange']
        for r_ in newR['reagents']:
            R.createReagent(r_, newR['reagents'][r_])
            mod.getSpecies(r_).setReagentOf(newR['id'])
        mod.addReaction(R)
        mod.createReactionLowerBound(newR['id'], newR['lower_bound'])
        mod.createReactionUpperBound(newR['id'], newR['upper_bound'])
        return R

class CBFileDialogue(QtGui.QWidget):
    _appTitle = 'Open file'
    work_dir = None
    model_file = None
    mode = None

    def __init__(self, work_dir, mode='open', filters=None):
        super(CBFileDialogue, self).__init__()
        self.mousePos = self.cursor().pos()
        self.work_dir = work_dir
        self.mode = mode
        if mode == 'save':
            self._appTitle = 'Save file'
        self.initUI()

    def initUI(self):
        #self.setAttribute(QtCore.Qt.WA_DeleteOnClose, True)
        #self.setHidden(True)
        self.__dlg__ = QtGui.QFileDialog(self)
        if self.mode == 'open':
            #self.model_file = str(self.__dlg__.getOpenFileName(self, 'Open file', self.work_dir, options=QtGui.QFileDialog.DontUseNativeDialog))
            self.model_file = str(self.__dlg__.getOpenFileName(self, 'Open file', self.work_dir))
        elif self.mode == 'save':
            #self.model_file = str(self.__dlg__.getSaveFileName(self, 'Save file as', self.work_dir, options=QtGui.QFileDialog.DontUseNativeDialog))
            self.model_file = str(self.__dlg__.getSaveFileName(self, 'Save file as', self.work_dir))
        self.model_file = os.path.normpath(self.model_file)


def fileDialogue(work_dir=None, mode='open', filters=None):
    if work_dir == None:
        work_dir = os.getcwd()
    if mode in ['open', 'save']:
        app = QtGui.QApplication([])
        fileApp = CBFileDialogue(work_dir, mode=mode, filters=filters)
        model_file = fileApp.model_file
        fileApp.__dlg__.done(1)
        app.exit()
        if mode == 'open':
            return model_file
        else:
            return True

class ViewSVG(QtGui.QWidget):
    _fixColour = None
    _errColour = None
    _goodColour = None
    _appTitle = 'ViewSVG'

    def __init__(self, filename):
        super(ViewSVG, self).__init__()
        self.mousePos = self.cursor().pos()
        self._fixColour = QtGui.QColor(0,0,153,alpha=255)
        self._errColour = QtGui.QColor(255,0,0,alpha=255)
        self._goodColour = QtGui.QColor(0,100,0,alpha=255)
        self.filename = os.path.abspath(filename)
        print('\nViewing file: {}'.format(filename))
        self.initUI()


    def initUI(self):
        # create panels
        self.txtId = QtGui.QTextEdit()

        # set up menu and status bar
        menuBar = QtGui.QMenuBar()

        exitAction = QtGui.QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit and loose changes')
        exitAction.triggered.connect(QtGui.qApp.quit)

        fileMenu = menuBar.addMenu('&File')
        fileMenu.addAction(exitAction)

        self.statusBar = QtGui.QStatusBar()
        self.statusBar.showMessage('Ready ...')

        # webkit panel
        scene = QtGui.QGraphicsScene()
        view = QtGui.QGraphicsView(scene)

        br = QtSvg.QGraphicsSvgItem(self.filename).boundingRect()

        webview = QGraphicsWebView()
        #webview.load(QtCore.QUrl("C:\your_interactive_svg.svg"))
        webview.load(QtCore.QUrl(QtCore.QUrl.fromLocalFile(self.filename)))
        webview.setFlags(QtGui.QGraphicsItem.ItemClipsToShape)
        webview.setCacheMode(QtGui.QGraphicsItem.NoCache)
        webview.resize(br.width(), br.height())

        scene.addItem(webview)
        view.resize(br.width()+10, br.height()+10)
        #view.show()


        # do layout
        grid = QtGui.QGridLayout()
        grid.setSpacing(10)

        grid.addWidget(menuBar, 0, 0, 1, 2)
        #grid.addWidget(self.txtId, 1, 0, 1, 2)
        grid.addWidget(view, 1, 0, 4, 4)
        grid.addWidget(self.statusBar, 5, 0, 1, 4)

        self.setLayout(grid)
        self.setGeometry(self.mousePos.x()-75, self.mousePos.y()-75, 500, 640)
        self.setWindowTitle(self._appTitle)
        self.show()

def loadViewSVG(filename):
    app = QtGui.QApplication([])
    ex = ViewSVG(filename)
    app.exec_()


class ValueSlider(QtGui.QWidget):
    _fixColour = None
    _errColour = None
    _goodColour = None
    _appTitle = 'ValueSlider'

    def __init__(self):
        super(ValueSlider, self).__init__()
        self.mousePos = self.cursor().pos()
        self._fixColour = QtGui.QColor(0,0,153,alpha=255)
        self._errColour = QtGui.QColor(255,0,0,alpha=255)
        self._goodColour = QtGui.QColor(0,100,0,alpha=255)
        self.initUI()


    def initUI(self):

        # create panels
        #self.txtId = QtGui.QTextEdit()

        l1a = QtGui.QLabel(self)
        l1a.setText('Property')
        sld1 = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        sld1.setTickPosition(sld1.TicksBelow)

        sld1_min = -100
        sld1_max = 100

        sld1.setMinimum(sld1_min)
        sld1.setMaximum(sld1_max)
        sld1.setTickInterval((sld1_min-sld1_max)/10.0)
        sld1.setSingleStep(0.1)

        sld1.setFocusPolicy(QtCore.Qt.NoFocus)
        sld1.valueChanged[int].connect(self.changeValue)
        self.l1b = QtGui.QLabel(self)
        self.l1b.setText('0.0')

        # do layout
        grid = QtGui.QGridLayout()
        grid.setSpacing(10)

        grid.addWidget(l1a, 0, 0, 1, 1)
        grid.addWidget(self.l1b, 0, 1, 1, 1)
        grid.addWidget(sld1, 0, 2, 1, 5)


        #grid.addWidget(self.txtId, 1, 0, 1, 2)
        #grid.addWidget(menuBar, 0, 0)
        #grid.addWidget(self.statusBar, 0, 1)

        self.setLayout(grid)
        self.setGeometry(self.mousePos.x()-75, self.mousePos.y()-75, 280, 170)
        self.setWindowTitle(self._appTitle)
        self.show()


    def changeValue(self, value):
        getattr(self, 'l1b').setText('{}'.format(value))

def loadSlider():
    app = QtGui.QApplication([])
    ex = ValueSlider()
    app.exec_()


data = "<DATASTART><return>{}</return>"

if __name__ == '__main__':
    print(os.sys.argv)
    if os.sys.argv[1] == 'fileOpen':
        filename = fileDialogue(work_dir=None, mode='open', filters=None)
        print(data.format(filename))
        os.sys.exit(0)

    # subprocess.check_output(['python', '_qtloader.py', 'fileOpen']).split('<DATASTART>')[1].strip()


# template widget
"""
    class SmallAppBasicGrid(QtGui.QWidget):
        _fixColour = None
        _errColour = None
        _goodColour = None
        _appTitle = 'SmallAppBaseGrid'

        def __init__(self):
            super(SmallAppBasicGrid, self).__init__()
            self.mousePos = self.cursor().pos()
            self._fixColour = QtGui.QColor(0,0,153,alpha=255)
            self._errColour = QtGui.QColor(255,0,0,alpha=255)
            self._goodColour = QtGui.QColor(0,100,0,alpha=255)
            self.initUI()


        def initUI(self):

            # create panels
            self.txtId = QtGui.QTextEdit()

            # set up menu and status bar
            menuBar = QtGui.QMenuBar()

            exitAction = QtGui.QAction('&Exit', self)
            exitAction.setShortcut('Ctrl+Q')
            exitAction.setStatusTip('Exit and loose changes')
            exitAction.triggered.connect(QtGui.qApp.quit)

            fileMenu = menuBar.addMenu('&File')
            fileMenu.addAction(exitAction)

            self.statusBar = QtGui.QStatusBar()
            self.statusBar.showMessage('Ready ...')

            # do layout
            grid = QtGui.QGridLayout()
            grid.setSpacing(10)

            grid.addWidget(menuBar, 0, 0, 1, 2)
            grid.addWidget(self.txtId, 1, 0, 1, 2)
            grid.addWidget(self.statusBar, 2, 0, 1, 2)

            self.setLayout(grid)
            self.setGeometry(self.mousePos.x()-75, self.mousePos.y()-75, 500, 640)
            self.setWindowTitle(self._appTitle)
            self.show()

    def loadBasicApp(mod):
        app = QtGui.QApplication([])
        ex = SmallAppBasicGrid()
        app.exec_()
"""

    # template microGUI
"""
    class MicroGUI(QtGui.QWidget):
        _appTitle = 'MicroGUI'

        def __init__(self):
            super(MicroGUI, self).__init__()
            self.mousePos = self.cursor().pos()
            self.work_dir = work_dir
            self.initUI()


        def initUI(self):
            self.setAttribute(QtCore.Qt.WA_DeleteOnClose, True)
            self.setHidden(True)

            # action code (this example is for opening a file dialogue)
            self.model_file = str(QtGui.QFileDialog.getOpenFileName(self, 'Open file', os.getcwd()))



    def loadMicroGUI(*args):
        app = QtGui.QApplication([])
        mGUI = OpenFileDialogue()
        appTitle = mGUI._appTitle
        del mGUI, app
        return
"""

"""
import sys
from PyQt4 import QtCore, QtGui, QtSvg
from PyQt4.QtWebKit import QGraphicsWebView
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)

    scene = QtGui.QGraphicsScene()
    view = QtGui.QGraphicsView(scene)

    br = QtSvg.QGraphicsSvgItem("C:\your_interactive_svg.svg").boundingRect()

    webview = QGraphicsWebView()
    webview.load(QtCore.QUrl("C:\your_interactive_svg.svg"))
    webview.load(QtCore.QUrl(QtCore.QUrl.fromLocalFile("C:\your_interactive_svg.svg")))
    webview.setFlags(QtGui.QGraphicsItem.ItemClipsToShape)
    webview.setCacheMode(QtGui.QGraphicsItem.NoCache)
    webview.resize(br.width(), br.height())

    scene.addItem(webview)
    view.resize(br.width()+10, br.height()+10)
    view.show()
    sys.exit(app.exec_())
"""