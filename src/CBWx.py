"""
CBMPy: CBWx module
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
Last edit: $Author: bgoli $ ($Id: CBWx.py 575 2017-04-13 12:18:44Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os, time, random, math, re, webbrowser, locale, sys
# this is a hack that needs to be streamlined a bit
try:
    import cStringIO as csio
except ImportError:
    import io as csio
HAVE_WX = False
try:
    import wx
    import wx.grid
    import wx.html
    import wx.lib.ClickableHtmlWindow
    import wx.aui
    import wx.lib.colourdb
    from wx.lib.floatcanvas import NavCanvas, FloatCanvas, Resources
    from wx.lib import scrolledpanel
    import  wx.lib.mixins.listctrl  as  listmix
    from wx import ImageFromStream, BitmapFromImage
    HAVE_WX = True
    print('WX GUI tools available.')
except ImportError as ex:
    print('\nWX GUI tools not available.')
    print(ex)

HAVE_URLLIB2 = True
try:
    import urllib2
except:
    HAVE_URLLIB2 = False

if HAVE_WX:
    class ModelEditor(wx.Frame):
        PanelSize = (1100,600)
        RPwidth = 0.5
        RPheight = 1.0
        LPwidth = 0.5
        LPheight = 1.0
        mainSizer = None
        RGridCol = None
        RGridRow = None
        rlabels = ('Reaction','Name','Flux','d','LB','UB','RCost','Exch','Balanced')
        NOVAL = (None, '', 'None')
        RGrid = None
        pybox = None
        NoteB1 = None
        NoteB1_Psession = None
        NoteB1_Preaction = None
        NoteB1_Prelate = None
        NoteB1_Pgene = None
        FCanvas_met = None
        FCanvas_gen = None
        Colour = None
        Rinfbox = None
        _cmod_ = None
        _cbm_ = None
        MainPanel = None
        StatusBar = None
        OUT_dir = None
        OUT_file = None
        STS_OBJ_FUNC = None
        STS_OBJ_SENSE = None
        STS_OBJ_COND = None
        BUT_minsum = None
        BUT_optimise = None
        RGridOK = True
        __ActiveReaction = None
        __ScaledReducedCost = False
        __BoundaryDetection = True
        WX_COLOURS = None
        FC_LineStyles = None
        __C_BUFF__ = None
        RGdict = None
        SGdict = None
        GGdict = None
        Rmap = None
        Smap = None
        GPRmap = None
        PRGmap = None
        OnlyActiveReactions = False
        ZERO_TOL = 1.0e-11
        SEARCH_MODE = 0
        RTreeCache = None
        ReactionBalanceInfo = None
        FontCache = None
        SemanticSBMLClient = None

        def __init__(self, cmod):
            mystyle = wx.DEFAULT_DIALOG_STYLE
            mystyle = wx.DEFAULT_FRAME_STYLE
            fSize = wx.DisplaySize()
            fSize = fSize[0], fSize[1]-50
            wx.Frame.__init__(self, None, style=mystyle, name='frame1', pos=(0, 0), size=fSize)
            #self.SetSize(self.PanelSize)
            #self.SetMinSize(self.PanelSize)
            #self.ShowFullScreen(True)
            self.SetTitle('PySCeS-CBM Model Editor - editing: %s (%s)' % (cmod.getId(), cmod.getName()))
            fSize = self.GetSize()
            fSize = fSize[0], fSize[1]-20
            self.MainPanel = wx.Panel(self, -1, size=fSize)
            self.MainPanel.SetSize(fSize)
            self.MainPanel.SetMinSize(fSize)

            self.PanelSize = self.MainPanel.GetSize()
            self._cmod_ = cmod
            self.__C_BUFF__ = []
            self.RGridRow = cmod.getReactionIds()
            self.RGdict = {}
            self.SGdict = {}
            self.GGdict = {}
            self.Rmap = {}
            self.Smap = {}
            self.GPRmap = {}
            self.PRGmap = {}
            self.RTreeCache = []
            self.ReactionBalanceInfo = {}
            self.FontCache = {}

            self.BuildMainPanel()
            self.CreateRGrid()
            #self.LeftPanel.Refresh()
            #self.MainPanel.UpdateWindowUI()

            import cbmpy as cbm
            self._cbm_ = cbm
            self.OUT_dir = os.getcwd()
            self.UpdateModelStatus()
            self.CreateMaps()
            try:
                #self.RESTClient = cbm.CBNetDB.RESTClient()
                self.SemanticSBMLClient = cbm.CBNetDB.SemanticSBML()
            except:
                print('REST web serices not available')

        def EVT_SEARCH_MODE_SELECT(self, event):
            if event.GetId() == 211:
                self.SEARCH_MODE = 0
                self.searchSelectmenu.Check(211, True)
                self.searchSelectmenu.Check(212, False)
                self.SEARCH_OBJ.Enable(True)
            elif event.GetId() == 212:
                self.SEARCH_MODE = 1
                self.searchSelectmenu.Check(211, False)
                self.searchSelectmenu.Check(212, True)
                self.SEARCH_OBJ.Enable(False)
            print('Search mode selected:', self.SEARCH_MODE, event.GetId())

        def BuildMainPanel(self):
            # build some frame stuff
            # Setting up the menu.
            filemenu= wx.Menu()
            menuSave = filemenu.Append(101, "&Save Model"," Save to SBML", kind=wx.ID_SAVE)
            menuExport = filemenu.Append(102, "Save S&ession"," Export session as Python script", kind=wx.ID_SAVE)
            menuAbout= filemenu.Append(103, "&About"," Information about this program", kind=wx.ID_ABOUT)
            menuExit = filemenu.Append(104,"E&xit"," Terminate the program", kind=wx.ID_EXIT)

            self.searchSelectmenu = wx.Menu()
            self.CHECK_AUTOCOMP = self.searchSelectmenu.Append(211, 'AutoComplete', kind=wx.ITEM_CHECK)
            self.CHECK_SUBSEARCH = self.searchSelectmenu.Append(212, 'Subsearch', kind=wx.ITEM_CHECK)

            self.Bind(wx.EVT_MENU, self.EVT_SEARCH_MODE_SELECT, self.CHECK_AUTOCOMP)
            self.Bind(wx.EVT_MENU, self.EVT_SEARCH_MODE_SELECT, self.CHECK_SUBSEARCH)
            self.searchSelectmenu.Check(self.CHECK_AUTOCOMP.GetId(), True)

            optionsMenu = wx.Menu()
            optionsMenu.AppendMenu(201, '&Search mode', self.searchSelectmenu)

            actionMenu = wx.Menu()
            menuAnalyseBalances = actionMenu.Append(301, 'Run &Balance Checker', kind=wx.ID_DEFAULT)
            SemSBMLMenu = wx.Menu()
            menuSemSBMLMenu_name = SemSBMLMenu.Append(311, 'Name', kind=wx.ID_DEFAULT)
            menuSemSBMLMenu_id = SemSBMLMenu.Append(312, 'Id', kind=wx.ID_DEFAULT)
            actionMenu.AppendMenu(302, 'Run Semantic &SBML', SemSBMLMenu)
            self.Bind(wx.EVT_MENU, self.MENUAnalyseBalances, menuAnalyseBalances)
            self.Bind(wx.EVT_MENU, self.SemSBML_name, menuSemSBMLMenu_name)
            self.Bind(wx.EVT_MENU, self.SemSBML_id, menuSemSBMLMenu_id)



            # Creating the menubar and statusbar
            menuBar = wx.MenuBar()
            menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
            menuBar.Append(optionsMenu,"&Options") # Adding the "filemenu" to the MenuBar
            menuBar.Append(actionMenu,"&Action") # Adding the "filemenu" to the MenuBar
            self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
            # Menu Events.
            self.Bind(wx.EVT_MENU, self.MENUOnSave, menuSave)
            self.Bind(wx.EVT_MENU, self.MENUExport, menuExport)
            self.Bind(wx.EVT_MENU, self.MENUOnExit, menuExit)
            self.Bind(wx.EVT_MENU, self.MENUOnAbout, menuAbout)

            self.StatusBar = self.CreateStatusBar() # A StatusBar in the bottom of the window
            self.SetStatusText("Ready!")

            # create the sizers
            self.mainSizer = wx.BoxSizer(wx.HORIZONTAL)

            # Checkbox
            ##  self.insure = wx.CheckBox(self, label="Do you want Insured Shipment ?")
            ##  grid.Add(self.insure, pos=(4,0), span=(1,2), flag=wx.BOTTOM, border=5)
            ##  self.Bind(wx.EVT_CHECKBOX, self.EvtCheckBox, self.insure)

            # Radio Boxes
            ##  radioList = ['blue', 'red', 'yellow', 'orange', 'green', 'purple', 'navy blue', 'black', 'gray']
            ##  rb = wx.RadioBox(self, label="What color would you like ?", pos=(20, 210), choices=radioList,  majorDimension=3,
                        ##  style=wx.RA_SPECIFY_COLS)
            ##  grid.Add(rb, pos=(5,0), span=(1,2))
            ##  self.Bind(wx.EVT_RADIOBOX, self.EvtRadioBox, rb)

            # create the Lpanel
            leftSizer = wx.BoxSizer(wx.VERTICAL)
            self.LeftPanel = wx.Panel(self.MainPanel, -1)
            self.RGridGap = 80
            gridSize = (self.PanelSize[0]*self.LPwidth, self.PanelSize[1]-self.RGridGap)
            self.RGridInitSize = gridSize
            self.LeftPanel.SetMinSize(gridSize)

            # Define a dynamically updatable grid
            self.RGrid_scrollwindow = wx.ScrolledWindow(id=-1, size=gridSize,
                                                        name='RGrid_scrollwindow', parent=self.LeftPanel,
                                                        style=wx.HSCROLL | wx.VSCROLL)

            self.RGrid = wx.grid.Grid(id=-1, name='RGrid', size=gridSize,
                                      parent=self.RGrid_scrollwindow,
                                      style=0)
            self.RGridCol = len(self.rlabels)

            self.RGrid.CreateGrid(len(self._cmod_.reactions), self.RGridCol)
            #self.RGrid.SetSize(gridSize)
            #self.RGrid.SetMinSize(gridSize)
            ##  self.RGrid.SetMargins(5,5)
            self.RGrid.EnableScrolling(True, True)
            self.RGrid.SetScrollbars(10, 10, 10, 10)

            self.RGrid.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.EvtRGridCellSelect)
            self.RGrid.Bind(wx.grid.EVT_GRID_CELL_CHANGE, self.EvtRGridCellChange)

            # Some buttons
            buttonSizer = wx.BoxSizer(wx.HORIZONTAL)

            ButtonPanel = wx.Panel(self.LeftPanel, -1)
            #ButtonPanel.SetSize((20,30))
            #ButtonPanel.SetMinSize((20,30))
            self.BUT_optimise = wx.Button(ButtonPanel, label="Optimize")
            self.BUT_minsum = wx.Button(ButtonPanel, label="Min. SumAbsFlux")
            self.Bind(wx.EVT_BUTTON, self.EVT_BUT_optimise, self.BUT_optimise)
            self.Bind(wx.EVT_BUTTON, self.EVT_BUT_minsum, self.BUT_minsum)
            buttonSizer.Add(self.BUT_optimise, 1, wx.CENTER)
            buttonSizer.Add(self.BUT_minsum, 1, wx.CENTER)
            ButtonPanel.SetSizer(buttonSizer)

            # add status text controls
            statusSizer = wx.GridSizer(2,4,1,1)
            StatusPanel = wx.Panel(self.LeftPanel, -1, size=wx.Size(self.PanelSize[0]*self.LPwidth,-1))
            STS_OBJ_FUNC_LBL = wx.TextCtrl(StatusPanel, -1, style=wx.TE_READONLY | wx.TE_CENTER)
            STS_OBJ_FUNC_LBL.write('ObjValue')
            STS_OBJ_FUNC_LBL.SetBackgroundColour(wx.Colour(255,255,153))
            statusSizer.Add(STS_OBJ_FUNC_LBL)
            STS_OBJ_SENSE_LBL = wx.TextCtrl(StatusPanel, -1, style=wx.TE_READONLY | wx.TE_CENTER)
            STS_OBJ_SENSE_LBL.write('ObjSense')
            STS_OBJ_SENSE_LBL.SetBackgroundColour(wx.Colour(255,255,153))
            statusSizer.Add(STS_OBJ_SENSE_LBL)
            STS_OBJ_COND_LBL = wx.TextCtrl(StatusPanel, -1, style=wx.TE_READONLY | wx.TE_CENTER)
            STS_OBJ_COND_LBL.write('ObjStatus')
            STS_OBJ_COND_LBL.SetBackgroundColour(wx.Colour(255,255,153))
            statusSizer.Add(STS_OBJ_COND_LBL)

            SEARCH_OBJ_LBL = wx.TextCtrl(StatusPanel, -1, style=wx.TE_READONLY | wx.TE_CENTER)
            SEARCH_OBJ_LBL.write('SearchReactions')
            SEARCH_OBJ_LBL.SetMinSize(wx.Size(130,-1))
            SEARCH_OBJ_LBL.SetBackgroundColour(wx.Colour(102,255,255))
            statusSizer.Add(SEARCH_OBJ_LBL)

            ## the combobox Control
            #self.SEARCH_SELECT_COMB = wx.ComboBox(StatusPanel, size=wx.Size(130, -1), choices=['AutoComplete', 'SubString'], style=wx.CB_READONLY)
            #self.SEARCH_SELECT_COMB.SetStringSelection('AutoComplete')
            #statusSizer.Add(self.SEARCH_SELECT_COMB)
            ## self.Bind(wx.EVT_COMBOBOX, self.EvtComboBox, self.edithear)

            self.STS_OBJ_FUNC = wx.TextCtrl(StatusPanel, -1, style=wx.TE_READONLY | wx.TE_CENTER)
            statusSizer.Add(self.STS_OBJ_FUNC)
            self.STS_OBJ_SENSE = wx.TextCtrl(StatusPanel, -1, style=wx.TE_READONLY | wx.TE_CENTER)
            statusSizer.Add(self.STS_OBJ_SENSE)
            self.STS_OBJ_COND = wx.TextCtrl(StatusPanel, -1, style=wx.TE_READONLY | wx.TE_CENTER)
            statusSizer.Add(self.STS_OBJ_COND)

            def SearchSelectCallback(sList):
                print('Selection list', sList)
                self.UpdateReactionGraph(sList[0])
                self.updateInfoFromReactionName(sList[0])
                self.SelectGridRow(sList[0])


            try:
                self.SEARCH_OBJ = TextCtrlAutoComplete(StatusPanel, choices=self._cmod_.getReactionIds(),\
                                                       selectCallback=SearchSelectCallback, style=wx.TE_LEFT)
                self.SEARCH_OBJ.SetMinSize(wx.Size(130, -1))
                self.SEARCH_OBJ.Enable(True)
            except:
                self.SEARCH_OBJ = wx.TextCtrl(StatusPanel, -1, choices=self._cmod_.getReactionIds(), style=wx.TE_READONLY | wx.TE_CENTER)
                self.SEARCH_OBJ.SetMinSize(wx.Size(130,-1))
                self.SEARCH_OBJ.Enable(True)
            statusSizer.Add(self.SEARCH_OBJ)
            StatusPanel.SetSizer(statusSizer)
            self.StatusPanel = StatusPanel

            leftSizer.Add(StatusPanel)
            leftSizer.AddSpacer(2)
            leftSizer.Add(ButtonPanel)
            leftSizer.AddSpacer(2)
            leftSizer.Add(self.RGrid_scrollwindow)
            self.LeftPanel.SetSizer(leftSizer)

            # create the Rpanel notebook
            rightSizer = wx.BoxSizer(wx.VERTICAL)
            self.RightPanel = wx.Panel(self.MainPanel, -1)
            NoteB1size = (self.PanelSize[0]*self.RPwidth, self.PanelSize[1]*self.RPheight)
            self.NoteB1 = wx.Notebook(self.RightPanel, size=NoteB1size)
            self.NoteB1.SetSize(NoteB1size)
            self.NoteB1.SetMinSize(NoteB1size)
            ##  self.NoteB1.SetPadding((10,-1))
            # create panels
            RinfSize = wx.Size(NoteB1size[0]-8, NoteB1size[1]-25)
            self.NoteB1_Psession =  wx.Panel(self.NoteB1, -1)
            self.NoteB1_Preaction = wx.Panel(self.NoteB1, -1)
            self.NoteB1_Pspecies = wx.Panel(self.NoteB1, -1)
            self.NoteB1_Prelate = wx.Panel(self.NoteB1, -1)
            self.NoteB1_Pgene = wx.Panel(self.NoteB1, -1)
            #self.NoteB1_PRedit = wx.Panel(self.NoteB1, -1)
            self.NoteB1_PRedit = wx.lib.scrolledpanel.ScrolledPanel(self.NoteB1, -1, style=wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER|wx.HSCROLL|wx.VSCROLL)
            self.NoteB1_PRedit.SetAutoLayout(1)
            self.NoteB1_PRedit.SetupScrolling()

            #create pybox: python script window
            self.PyBox = wx.TextCtrl(self.NoteB1_Psession, -1, size=self.NoteB1.GetVirtualSize(), style=wx.TE_MULTILINE | wx.HSCROLL)
            PyBox_sizer = wx.BoxSizer(wx.VERTICAL)
            PyBox_sizer.Add(self.PyBox, -1, wx.EXPAND)
            #self.NoteB1_Psession.Fit(PyBox_sizer)
            self.PyBox.SetEditable(False)
            self.writeCmd("#############################\n# ")
            self.writeCmd("# PySCeS-CBM GUI generated command file ")
            self.writeCmd("# Please note that the model instance is: cmod ")
            self.writeCmd("#\n#############################\n")
            self.writeCmd("import cbmpy\n")
            self.writeCmd("# cmod = cbmpy.CBRead.readSBML2FBA(ModelFile, ModelDir)\n")
            self.PyBox.SetMinSize(self.NoteB1.GetVirtualSize())

            # create html panel for reaction information
            self.Rinfbox = wx.html.HtmlWindow(self.NoteB1_Preaction, -1, size=self.NoteB1.GetVirtualSize(), style=wx.html.HW_SCROLLBAR_AUTO)
            self.Rinfbox.SetBorders(0)
            Rinfbox_sizer = wx.BoxSizer(wx.VERTICAL)
            Rinfbox_sizer.Add(self.Rinfbox, 0, wx.EXPAND)
            self.Rinfbox.SetMinSize(self.NoteB1.GetVirtualSize())

            # create reaction editor
            self.FontCache.update({'SanSer12CS' : wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL, False, u'Comic Sans MS')})
            self.FontCache.update({'SanSer12' : wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL, False)})
            self.FontCache.update({'SanSer10' : wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL, False)})

            # create the box sizer
            Reditbox_sizer = wx.BoxSizer(wx.VERTICAL)

            # create the edit button
            #self.BUT_EnableEdit = wx.Button(self.NoteB1_PRedit, label="Edit Reaction")
            #self.Bind(wx.EVT_BUTTON, self.EVT_BUT_EnableEdit, self.BUT_EnableEdit)
            #Reditbox_sizer.Add(self.BUT_EnableEdit)

            # create the component grid
            #Reditgrid_sizer = wx.GridSizer(3,2,1,1)
            #grow = 10
            #gcol = 2
            #cntr = True
            #self.TEXTBOXES_Redit = []
            #for tp in range(grow*gcol):
                #cellHeight = -1
                #if cntr:
                    #tx = wx.TextCtrl(self.NoteB1_PRedit, style=wx.DEFAULT | wx.TE_READONLY | wx.TE_CENTER)
                    #tx.SetSize(wx.Size(RinfSize[1]/4, cellHeight))
                    #tx.SetMinSize(wx.Size(RinfSize[1]/4, cellHeight))
                    #tx.SetFont(self.FontCache['SanSer10'])
                    #cntr = False
                #else:
                    #tx = wx.TextCtrl(self.NoteB1_PRedit, style=wx.DEFAULT | wx.TE_LEFT | wx.EXPAND )
                    #tx.SetSize(wx.Size(RinfSize[1]/2, cellHeight))
                    #tx.SetMinSize(wx.Size(RinfSize[1]/2, cellHeight))
                    #tx.SetFont(self.FontCache['SanSer10'])
                    #cntr = True
                #tx.Disable()
                #tx.WriteText(str(tp+1))
                #tx.PID = str(tp+1)
                #Reditgrid_sizer.Add(tx, wx.EXPAND)
                #self.TEXTBOXES_Redit.append(tx)
            #Reditbox_sizer.Add(Reditgrid_sizer, wx.EXPAND) # add to boxsizer

            # create the annotate button
            #self.BUT_Annotate = wx.Button(self.NoteB1_PRedit, label="SemanticSBML")
            #self.Bind(wx.EVT_BUTTON, self.EVT_BUT_Annotate, self.BUT_Annotate)
            #Reditbox_sizer.Add(self.BUT_Annotate)

            # the big text box
            #self.TEXTBOX_Annotate = wx.TextCtrl(self.NoteB1_PRedit, style=wx.TE_MULTILINE | wx.EXPAND )
            #self.TEXTBOX_Annotate.SetMinSize(wx.Size(self.NoteB1.GetVirtualSize()[1], -1))
            #self.TEXTBOX_Annotate.SetFont(self.FontCache['SanSer10'])
            #self.TEXTBOX_Annotate.Disable()
            #Reditbox_sizer.Add(self.TEXTBOX_Annotate, wx.EXPAND) # add to boxsizer

            # add boxSizer to Panel
            self.NoteB1_PRedit.SetSizer(Reditbox_sizer, wx.DEFAULT)
            self.NoteB1_PRedit.SetMinSize(self.NoteB1.GetVirtualSize())

            # create html panel for reagent information
            self.Sinfbox = wx.lib.ClickableHtmlWindow.PyClickableHtmlWindow(self.NoteB1_Pspecies, -1, size=self.NoteB1.GetVirtualSize(), style=wx.html.HW_SCROLLBAR_AUTO)
            self.Sinfbox.SetBorders(0)
            Sinfbox_sizer = wx.BoxSizer(wx.VERTICAL)
            Sinfbox_sizer.Add(self.Sinfbox, -1, wx.EXPAND)
            self.Sinfbox.SetMinSize(self.NoteB1.GetVirtualSize())

            # Create canvas panels for graphical representations
            wx.lib.colourdb.updateColourDB()
            self.WX_COLOURS = wx.lib.colourdb.getColourList()
            self.FC_LineStyles = list(FloatCanvas.DrawObject.LineStyleList)

            # create a met/gen canvas and do voodoo to get it to fill the notebook panel
            self.MetCanvas = NavCanvas.NavCanvas(self.NoteB1_Prelate, -1, size=self.NoteB1.GetVirtualSize(), Debug=0, ProjectionFun=None, BackgroundColor="DARK SLATE BLUE")
            self.FCanvas_met = self.MetCanvas.Canvas
            self.FCanvas_met.SetProjectionFun(None)
            self.FCanvas_met.InitAll()
            self.GenCanvas = NavCanvas.NavCanvas(self.NoteB1_Pgene, -1, size=self.NoteB1.GetVirtualSize(), Debug=0, ProjectionFun=None, BackgroundColor="DARK SLATE BLUE")
            self.FCanvas_gen = self.GenCanvas.Canvas
            self.FCanvas_gen.SetProjectionFun(None)
            self.FCanvas_gen.InitAll()

            NC1_sizer = wx.BoxSizer(wx.VERTICAL)
            NC1_sizer.Add(self.MetCanvas, 1, wx.EXPAND)
            NC2_sizer = wx.BoxSizer(wx.VERTICAL)
            NC2_sizer.Add(self.GenCanvas, 1, wx.EXPAND)
            self.FCanvas_met.SetMinSize(self.NoteB1.GetVirtualSize())
            self.FCanvas_gen.SetMinSize(self.NoteB1.GetVirtualSize())

            ##  self.FCanvas_met.Bind(FloatCanvas.EVT_LEFT_DOWN, self.EVT_FC1_onClick)
            ##  self.FCanvas_gen.Bind(FloatCanvas.EVT_LEFT_DOWN, self.EVT_FC2_onClick)

            # add the pages to the notebook with the label to show on the tab
            self.NoteB1.AddPage(self.NoteB1_Psession, "Session")
            self.NoteB1.AddPage(self.NoteB1_Preaction, "Reaction")
            self.NoteB1.AddPage(self.NoteB1_Prelate, "Metabolism")
            # TODO
            self.NoteB1.AddPage(self.NoteB1_Pgene, "Genes")
            self.NoteB1.AddPage(self.NoteB1_PRedit, "ReacEdt")
            self.NoteB1.AddPage(self.NoteB1_Pspecies, "MIRIAM")

            rightSizer.Add(self.NoteB1)
            #self.RightPanel.Fit(rightSizer)

            # Add components to the main panel

            self.mainSizer.Add(self.LeftPanel, 1, wx.EXPAND)
            self.mainSizer.Add(self.RightPanel, 1, wx.EXPAND)
            self.MainPanel.SetSizerAndFit(self.mainSizer)

            ## Resize control
            self.Bind(wx.EVT_SIZE, self.EVT_FRAME_resize, self)

            ##GENERIC click event
            #self.Bind(wx.EVT_SIZE, self.EVT_onClick, self.MainPanel)

        def UpdateModelStatus(self):
            self.STS_OBJ_SENSE.Clear()
            self.STS_OBJ_COND.Clear()
            self.STS_OBJ_FUNC.Clear()
            self.STS_OBJ_COND.write(str(self._cmod_.SOLUTION_STATUS))
            if self._cmod_.SOLUTION_STATUS == 'LPS_OPT':
                self.STS_OBJ_COND.SetBackgroundColour(wx.Colour(255,255,255))
            else:
                self.STS_OBJ_COND.SetBackgroundColour(wx.Colour(255,0,51))
            if self._cmod_.getActiveObjective() != None:
                self.STS_OBJ_FUNC.write(str(self._cmod_.getActiveObjective().getValue()))
                self.STS_OBJ_SENSE.write(self._cmod_.getActiveObjective().operation)
            else:
                self.STS_OBJ_FUNC.write('None')
                self.STS_OBJ_SENSE.write('None')

        def UpdateReactionInfo(self,rid):
            print('Updating reaction info')
            r = self._cmod_.getReaction(rid)

            rs = '<h2>%s</h2><p>%s</p>' % (rid, r.getName())
            rs += '<h3>Equation</h3><p>%s</p>' % self.GetEquation(r)
            rs += '<p><br/></p>'
            ##  rs += '<h3>Details</h3>'
            props = '<tr><th colspan="2"><strong>Properties</strong></th></tr>'
            props += '<tr><td>%s</td><td>%s</td></tr><tr><td>%s</td><td>%s</td></tr>' % ('Reversible',r.reversible,'Exchange',r.is_exchange)
            props += '<tr><th colspan="2"><strong>Annotations</strong></th></tr>'
            RA = r.getAnnotations()
            for a in RA:
                props += '<tr><td>%s</td><td>%s</td></tr>' % (a, RA[a])
            ##  rs += '<table border="1" cellpadding="5" width="70%s">%s</table>' % ('\\%', props)
            rs += '<table border="1" cellpadding="5">%s</table>' % (props)

            if self.ReactionBalanceInfo != None:
                if rid in self.ReactionBalanceInfo:
                    rs += '<h3>Balancing information</h3>'
                    bi = self.ReactionBalanceInfo[rid]
                    #print bi
                    bal = '<tr><td colspan="2">Charge balanced: %s {\'charge\', %s}</td></tr>' % (bi['charge_balanced'], bi['charge'])
                    bal += '<tr><td colspan="2">Element balanced: %s %s</td></tr>' % (bi['element_balanced'], bi['elements'])
                    for rre in bi['stuff']:
                        sbal = ''
                        coeff = None
                        for det in range(len(rre)):
                            out = ''
                            if det == 0:
                                out = rre[det]
                            if det == 1:
                                coeff = rre[det]
                                #out = re[det]
                            elif det == 3:
                                out = ''
                                if rre[2] == '' or rre[3] == None:
                                    out = 'Unknown'
                                else:
                                    for e in rre[det]:
                                        out += '%.1f %s, ' % (coeff*e[1], e[0])
                                    out = out[:-2]
                            if out != '':
                                sbal += '<td>%s</td>' % out
                        bal += '<tr>%s</tr>' % (sbal)
                    rs += '<table border="1" cellpadding="5">%s</table>' % (bal)
            props = ''
            rs += '<h3>Reagents</h3>'
            for sid in r.getSpeciesIds():
                s = self._cmod_.getSpecies(sid)
                props += '<h4>%s</h4><table border="1" cellpadding="5">' % sid
                ##  props += '<tr><th colspan="2"><strong>%s</strong></th></tr>' % sid
                props += '<tr><td>%s</td><td>%s</td></tr><tr><td>%s</td><td>%s</td></tr>' % ('Name',s.getName(),'Compartment',s.compartment)
                props += '<tr><td>%s</td><td>%s</td></tr>' % ('Fixed',s.is_boundary)
                props += '<tr><td>%s</td><td>%s</td></tr>' % ('Coefficient', r.getReagentWithSpeciesRef(sid).coefficient)
                ##  props += '<tr><th colspan="2">Annotations</th></tr>'
                props += '<tr><td>%s</td><td>%s</td></tr><tr><td>%s</td><td>%s</td></tr>' % ('ChemFormula',s.chemFormula,'Charge',s.charge)
                ##  SA = s.getAnnotations()
                ##  for a in RA:
                    ##  props += '<tr><td>%s</td><td>%s</td></tr>' % (a, SA[a])
                ##  rs += '<table border="1" cellpadding="5" width="70%s">%s</table>' % ('\\%', props)
                props += '</table>'
            rs += props
                ##  rs += '<table border="1" cellpadding="3" width="80%s"><tr><th>Property</th><th>Value</th></tr>%s</table>' % ('\\%', props)


            self.Rinfbox.SetPage("<html><body>%s</body></html>" % rs)

        def SelectGridRow(self, rid):
            self.RGrid.SelectRow(self.RGridRow.index(rid))
            self.RGrid.MakeCellVisible(self.RGridRow.index(rid), 0)

        def UpdateReactionGraph(self, rid):
            R = self._cmod_.getReaction(rid)
            subs = []
            prods = []
            self.FCanvas_met.ClearAll()
            STB_LW = 2
            STB_PAD = 0.15
            STB_SZ = 0.25
            STB_LC = 'Black'
            GC = (0,0)

            FC_OBJ = []
            if rid not in self.RGdict:
                STB = FloatCanvas.ScaledTextBox(rid, GC, Size=STB_SZ, PadSize=STB_PAD,\
                                                LineWidth=STB_LW, LineColor=STB_LC,\
                                                Family=wx.FONTFAMILY_DEFAULT, Weight=wx.BOLD,\
                                                BackgroundColor='White')
                self.__STB_center__(STB)
                STB.Name = rid
                STB.HitFill = True
                STB.HitLineWidth = 5
                FC_OBJ.append(STB)

                for rr in R.reagents:
                    if rr.coefficient < 0:
                        subs.append(rr)
                    elif rr.coefficient > 0:
                        prods.append(rr)
                radius = None
                if len(R.reagents) == 1:
                    radius = 2
                    cxy = [(0,radius)]
                elif len(R.reagents) == 2:
                    radius = 2
                    cxy = circlePoints(totalPoints=len(R.reagents), startAngle=0, arc=360, circleradius=radius, centerxy=GC, direction='forward', evenDistribution=True)
                else:
                    if len(R.reagents) < 12:
                        radius = 3
                    else:
                        radius = 10
                    cxy = circlePoints(totalPoints=len(R.reagents), startAngle=270, arc=360, circleradius=radius, centerxy=GC, direction='forward', evenDistribution=True)

                scntr = 0
                rcntr = 0
                for s in range(len(subs)):
                    if self._cmod_.getSpecies(subs[s].species_ref).is_boundary:
                        STB_LC = 'Red'
                    else:
                        STB_LC = 'Black'
                    if abs(subs[s].coefficient) != 1:
                        cf = '{%.1f} ' % abs(subs[s].coefficient)
                    else:
                        cf = ''
                    STB = FloatCanvas.ScaledTextBox('%s%s' % (cf,subs[s].species_ref), cxy[rcntr], Size=STB_SZ, PadSize=STB_PAD, LineWidth=STB_LW, LineColor=STB_LC, Family=wx.SWISS, BackgroundColor='Green')
                    STB.Name = subs[s].species_ref
                    self.__STB_center__(STB)
                    STB.HitFill = True
                    STB.HitLineWidth = 5
                    FC_OBJ.append(STB)
                    scntr += 1
                    rcntr += 1
                pcntr = 0
                for p in range(len(prods)):
                    if self._cmod_.getSpecies(prods[p].species_ref).is_boundary:
                        STB_LC = 'Red'
                    else:
                        STB_LC = 'Black'
                    if abs(prods[p].coefficient) != 1:
                        cf = '{%.1f} ' % abs(prods[p].coefficient)
                    else:
                        cf = ''
                    STB = FloatCanvas.ScaledTextBox('%s%s' % (cf, prods[p].species_ref), cxy[rcntr], Size=STB_SZ, PadSize=STB_PAD, LineWidth=STB_LW, LineColor=STB_LC, Family=wx.SWISS, BackgroundColor='Yellow')
                    STB.Name = prods[p].species_ref
                    self.__STB_center__(STB)
                    STB.HitFill = True
                    STB.HitLineWidth = 5
                    FC_OBJ.append(STB)
                    pcntr += 1
                    rcntr += 1
                self.RGdict.update({rid : {'obj' : FC_OBJ}})
                ##  print 'Drawing from scratch'
            else:
                ##  print 'Using cached graph'
                FC_OBJ = self.RGdict[rid]['obj']
            self.FCanvas_met.AddObjects(FC_OBJ)
            for o in FC_OBJ:
                o.Bind(FloatCanvas.EVT_FC_LEFT_DOWN, self.HIT_STB_onClick)
            self.FCanvas_met.Draw(True)
            self.FCanvas_met.ZoomToBB()

        def __STB_center__(self, stb):
            box = stb.GetBoxRect()
            stb.SetPoint( ( box[0][0] - abs(box[1][0])/2.0, box[0][1] - abs(box[1][1])/2.0) )
            return stb.GetBoxRect()

        def UpdateSpeciesGraph(self, sid):
            #print sid, self.Smap[sid]
            STB_LW = 2
            STB_PAD = 0.15
            STB_SZ = 0.25
            STB_LC = 'Black'
            GC = (0,0)
            radius = None

            rList = self.Smap[sid]
            if self.OnlyActiveReactions:
                rList = [j for j in self.Smap[sid] if round(abs(self._cmod_.getReaction(j).getValue()),6) != 0.0]

            self.FCanvas_met.ClearAll()
            FC_OBJ = []
            if sid not in self.SGdict:
                orphan = False
                print('ReactionList', rList)
                if len(rList) == 0:
                    radius = 2
                    cxz = GC
                elif len(rList) == 1:
                    radius = 2
                    cxy = [(0,radius)]
                elif len(rList) == 2:
                    radius = 2
                    cxy = circlePoints(totalPoints=len(rList), startAngle=0, arc=360, circleradius=radius, centerxy=GC, direction='forward', evenDistribution=True)
                else:
                    if len(rList) < 12:
                        radius = 3
                    elif len(rList) < 24:
                        radius = 6
                    elif len(rList) < 36:
                        radius = 9
                    else:
                        radius = 15
                    cxy = circlePoints(totalPoints=len(rList), startAngle=270, arc=360, circleradius=radius, centerxy=GC, direction='forward', evenDistribution=True)
                    print(cxy)

                STB = FloatCanvas.ScaledTextBox(sid, GC, Size=STB_SZ, PadSize=STB_PAD, LineWidth=STB_LW, LineColor=STB_LC, Family=wx.NORMAL, Weight=wx.BOLD, BackgroundColor='White')
                self.__STB_center__(STB)
                STB.Name = sid
                STB.HitFill = True
                STB.HitLineWidth = 5
                FC_OBJ.append(STB)

                rcntr = 0
                for R in range(len(rList)):
                    cf = ''
                    STB = FloatCanvas.ScaledTextBox('%s%s' % (cf, rList[R]), cxy[R], Size=STB_SZ, PadSize=STB_PAD, LineWidth=STB_LW, LineColor=STB_LC, Family=wx.SWISS, BackgroundColor='Yellow')
                    STB.Name = self.Smap[sid][R]
                    self.__STB_center__(STB)
                    STB.HitFill = True
                    STB.HitLineWidth = 5
                    FC_OBJ.append(STB)
                    rcntr += 1
                self.SGdict.update({sid : {'obj' : FC_OBJ}})
                ##  print 'Drawing from scratch'
            else:
                ##  print 'Using cached graph'
                FC_OBJ = self.SGdict[sid]['obj']
            self.FCanvas_met.AddObjects(FC_OBJ)
            for o in FC_OBJ:
                o.Bind(FloatCanvas.EVT_FC_LEFT_DOWN, self.HIT_STB_onClick)
            self.FCanvas_met.Draw(True)
            self.FCanvas_met.ZoomToBB()

        def CreateMaps(self):
            for S in self._cmod_.species:
                self.Smap.update({S.getId() : S.isReagentOf()})
            for R in self._cmod_.reactions:
                self.Rmap.update({R.getId() : R.getSpeciesIds()})
            self.GPRmap = self._cmod_.getAllGeneProteinAssociations()
            self.PRGmap = self._cmod_.getAllProteinGeneAssociations()

        def GetEquation(self, R):
            sub = ''
            prod = ''
            for r in R.reagents:
                coeff = abs(r.coefficient)
                if r.role == 'substrate':
                    if coeff == 1.0:
                        sub += '%s + ' % (r.species_ref)
                    else:
                        sub += '{%s} %s + ' % (coeff, r.species_ref)
                else:
                    if coeff == 1.0:
                        prod += '%s + ' % (r.species_ref)
                    else:
                        prod += '{%s} %s + ' % (coeff, r.species_ref)
            if R.reversible:
                eq = '%s = %s' % (sub[:-3], prod[:-2])
            else:
                eq = '%s > %s' % (sub[:-3], prod[:-2])
            return eq

        def CreateRGrid(self):
            assert len(self.rlabels) == self.RGridCol, '\nlabels != #col'
            grid = self.RGrid
            reactions = self._cmod_.reactions
            fluxObjs = []
            if self._cmod_.getActiveObjective() != None:
                fluxObjs = self._cmod_.getActiveObjective().getFluxObjectiveReactions()
            for c in range(self.RGridCol):
                grid.SetColLabelValue(c, self.rlabels[c])
            for r in range(len(reactions)):
                for c in range(self.RGridCol):
                    ##  print r,c
                    ##  print reactions[r].getId()
                    if self.rlabels[c] == 'LB':
                        ##  grid.SetColFormatFloat(c, 10, 3)
                        grid.SetCellValue(r, c, '%s' % self._cmod_.getReactionLowerBound(reactions[r].getId()))
                    elif self.rlabels[c] == 'UB':
                        ##  grid.SetColFormatFloat(c, 10, 3)
                        grid.SetCellValue(r, c, '%s' % self._cmod_.getReactionUpperBound(reactions[r].getId()))
                    elif self.rlabels[c] == 'd':
                        grid.SetCellValue(r, c, '%s' % '   ')
                    elif self.rlabels[c] == 'Balanced':
                        grid.SetReadOnly(r,c,True)
                        txt = reactions[r].is_balanced
                        grid.SetCellValue(r, c, '  %s' % str(txt) )
                        if txt:
                            grid.SetCellBackgroundColour(r, c, wx.Colour(255,193,96))
                    elif self.rlabels[c] == 'Flux':
                        Rval = reactions[r].getValue()
                        if Rval != None:
                            Rval = round(Rval, 10)
                        else:
                            Rval = 'None'
                        grid.SetCellValue(r, c, '  %s' % Rval)
                        grid.SetReadOnly(r,c,True)
                        # colour by sign
                        if Rval == None or Rval == '' or Rval == 'None' or Rval == 0.0:
                            grid.SetCellBackgroundColour(r,c,wx.Colour(255,255,255))
                        elif Rval < 0.0:
                            grid.SetCellBackgroundColour(r,c,wx.Colour(255,204,204))
                        elif Rval > 0.0:
                            grid.SetCellBackgroundColour(r,c,wx.Colour(153,255,153))
                        # boundary detection
                    elif self.rlabels[c] == 'Exch':
                        grid.SetReadOnly(r,c,True)
                        grid.SetCellValue(r, c, str(reactions[r].is_exchange))
                    elif self.rlabels[c] == 'RCost':
                        grid.SetReadOnly(r,c,True)
                        if not self.__ScaledReducedCost:
                            if reactions[r].reduced_cost != None:
                                rcval = '%2.3e' % reactions[r].reduced_cost
                            else:
                                rcval = 'None'
                            grid.SetCellValue(r, c, rcval)
                    elif self.rlabels[c] == 'Reaction':
                        grid.SetCellBackgroundColour(r,c,wx.Colour(198,226,255))
                        grid.SetCellValue(r, c, str(reactions[r].getId()))
                        grid.SetReadOnly(r,c,True)
                    elif self.rlabels[c] == 'Name':
                        ##  grid.SetCellBackgroundColour(r,c,wx.Colour(255,255,153))
                        grid.SetCellValue(r, c, str(reactions[r].getName()))
                    else:
                        grid.SetCellValue(r, c, '')
                    if reactions[r].getId() in fluxObjs:
                        grid.SetCellBackgroundColour(r,c,wx.Colour(255,255,153))
            for c in range(self.RGridCol):
                if self.rlabels[c] not in ('Name'):
                    grid.AutoSizeColumn(c, True)
                else:
                    grid.SetColMinimalWidth(c, 40)

        def UpdateRGridData(self):
            """
            Updates Rgrid numeric data
            """
            grid = self.RGrid
            reactions = self._cmod_.reactions
            if self._cmod_.getActiveObjective() != None:
                fluxObjs = self._cmod_.getActiveObjective().getFluxObjectiveReactions()
            else:
                fluxObjs = []

            for r in range(len(reactions)):
                for c in range(self.RGridCol):
    ##                 if self.rlabels[c] == 'LB':
    ##                     grid.SetCellValue(r, c, '%s' % self._cmod_.getFluxBoundByReactionID(reactions[r].getId(), 'lower').getValue())
    ##                 elif self.rlabels[c] == 'UB':
    ##                     grid.SetCellValue(r, c, '%s' % self._cmod_.getFluxBoundByReactionID(reactions[r].getId(), 'upper').getValue())
                    if self.rlabels[c] == 'Flux':
                        Rval = reactions[r].getValue()
                        try:
                            Rval_curr = float(grid.GetCellValue(r, c))
                        except:
                            Rval_curr = None
                        if Rval != None:
                            Rval = round(Rval, 10)
                        else:
                            Rval = 'None'
                        grid.SetCellValue(r, c, '  %s' % Rval)
                        # colour delta
                        D_col_idx = self.rlabels.index('d')
                        if Rval in self.NOVAL or Rval_curr in self.NOVAL or abs(Rval - Rval_curr) <= self.ZERO_TOL:
                            grid.SetCellValue(r, D_col_idx, '  ')
                            grid.SetCellTextColour(r, D_col_idx, wx.Colour(255,255,255))
                        elif abs(round(Rval,10)) < abs(round(Rval_curr,10)):
                            grid.SetCellValue(r, D_col_idx, ' \/ ')
                            grid.SetCellTextColour(r, D_col_idx, wx.Colour(255,0,0))
                        elif abs(round(Rval,10)) > abs(round(Rval_curr,10)):
                            grid.SetCellValue(r, D_col_idx, ' /\ ')
                            grid.SetCellTextColour(r, D_col_idx, wx.Colour(0,204,0))
                        # colour by sign J
                        if Rval in self.NOVAL or Rval == 0.0:
                            grid.SetCellBackgroundColour(r,c,wx.Colour(255,255,255))
                        elif Rval < 0.0:
                            grid.SetCellBackgroundColour(r,c,wx.Colour(255,204,204))
                        elif Rval > 0.0:
                            grid.SetCellBackgroundColour(r,c,wx.Colour(153,255,153))

                        # boundary detection
                        if self.__BoundaryDetection:
                            LB = self._cmod_.getReactionLowerBound(reactions[r].getId())
                            UB = self._cmod_.getReactionUpperBound(reactions[r].getId())
                            try:
                                if abs(Rval - round(LB, 10)) <= self.ZERO_TOL:
                                    grid.SetCellBackgroundColour(r,self.rlabels.index('LB'), wx.Colour(255,204,204))
                                else:
                                    grid.SetCellBackgroundColour(r,self.rlabels.index('LB'), wx.Colour(255,255,255))
                            except:
                                print('INFO: LowerBound detector failed')
                            try:
                                if abs(Rval - round(UB, 10)) <= self.ZERO_TOL:
                                    grid.SetCellBackgroundColour(r,self.rlabels.index('UB'), wx.Colour(153,255,153))
                                else:
                                    grid.SetCellBackgroundColour(r,self.rlabels.index('UB'), wx.Colour(255,255,255))
                            except:
                                print('INFO: UpperBound detector failed')
                    elif self.rlabels[c] == 'RCost':
                        if not self.__ScaledReducedCost:
                            if reactions[r].reduced_cost != None:
                                rcval = '%2.3e' % reactions[r].reduced_cost
                            else:
                                rcval = 'None'
                            grid.SetCellValue(r, c, rcval)
                    elif self.rlabels[c] == 'Balanced':
                        bval = reactions[r].is_balanced
                        grid.SetCellValue(r, c, str(bval))
                        if bval != None and not bval:
                            grid.SetCellBackgroundColour(r, c, wx.Colour(255,193,96))

                    if reactions[r].getId() in fluxObjs:
                        grid.SetCellBackgroundColour(r,c,wx.Colour(255,255,153))
            grid.ForceRefresh()

        ##  def EvtGridCellSelect(self, event):
            ##  row = event.GetRow()
            ##  col = event.GetCol()
            ##  self.cell_selected_value = self.RGrid.GetCellValue(row, col)
            ##  print 'selected value', self.cell_selected_value
            ##  event.Skip()

        def writeCmd(self, txt):
            print('cmd: '+txt)
            self.PyBox.write(txt+'\n')


        def MENUOnAbout(self,e):
            # Create a message dialog box
            dlg = wx.MessageDialog(self, " PySCes-CBM model editor\n(C) Brett G. Olivier, Amsterdam 2012", "About PySCes-CBM model editor", wx.OK)
            dlg.ShowModal() # Shows it
            dlg.Destroy() # finally destroy it when finished.

        def MENUOnExit(self,e):
            self.Close(True)  # Close the frame.

        def MENUExport(self,e):
            dlg = wx.FileDialog(self, "Enter filename", self.OUT_dir, "", "Python files (*.py)|*.py|All files (*.*)|*.*", wx.FD_SAVE)
            if dlg.ShowModal() == wx.ID_OK:
                OUT_file = dlg.GetFilename()
                self.OUT_dir = dlg.GetDirectory()
                if OUT_file[-3:] != '.py':
                    OUT_file += '.py'
                F = file(os.path.join(self.OUT_dir, OUT_file), 'w')
                F.write(self.PyBox.GetValue())
                F.flush()
                F.close()
            dlg.Destroy()

        def MENUOnSave(self,e):
            """ Open a file"""
            dlg = wx.FileDialog(self, "Enter filename", self.OUT_dir, "", "SBML files (*.xml)|*.xml|All files (*.*)|*.*", wx.FD_SAVE)
            if dlg.ShowModal() == wx.ID_OK:
                self.OUT_file = dlg.GetFilename()
                self.OUT_dir = dlg.GetDirectory()
                if self.OUT_file[-4:] != '.xml':
                    self.OUT_file += '.xml'
                self._cmod_.inputfile_id = self.OUT_file
                self._cbm_.CBWrite.writeSBML2FBA(self._cmod_, self.OUT_file, directory=self.OUT_dir, sbml_level_version=None)
            dlg.Destroy()

        def MENUAnalyseBalances(self, e):
            ''' Check the reaction balances'''
            wait = wx.BusyCursor()
            rids = self._cmod_.getReactionIds()
            self.ReactionBalanceInfo = self._cbm_.CBTools.checkReactionBalanceElemental(self._cmod_, Rid=rids)
            for r in range(len(self._cmod_.reactions)):
                for c in range(self.RGridCol):
                    if self.rlabels[c] == 'Balanced':
                        bval = self._cmod_.reactions[r].is_balanced
                        self.RGrid.SetCellValue(r, c, str(bval))
                        if bval != None and not bval:
                            self.RGrid.SetCellBackgroundColour(r, c, wx.Colour(255,193,96))
            del wait
            wx.MessageBox('Balance Check Complete', 'Info', wx.OK | wx.ICON_INFORMATION)

        def chkFloat(self, val):
            try:
                x = float(val)
                self.StatusBar.SetStatusText('')
                return True
            except:
                self.StatusBar.SetStatusText('Invalid input %s was not a float!' % val)
                return False

        def EvtRGridCellSelect(self, event):
            print("Cell select event")
            row = event.GetRow()
            if row != self.__ActiveReaction:
                cval = self.RGrid.GetCellValue(row, 0)
                self.updateInfoFromReactionName(cval)
            event.Skip()

        def updateInfoFromReactionName(self, cval):
            self.UpdateReactionInfo(cval)
            #self.UpdateSpeciesInfoForReaction(cval)
            self.UpdateReactionGraph(cval)
            self.__ActiveReaction = cval

        def EvtRGridCellChange(self, event):
            row = event.GetRow()
            col = event.GetCol()
            clbl = self.RGrid.GetColLabelValue(col)
            #All cells have a value, regardless of the editor.
            print('Changed cell: ({}, {})'.format(row, col))
            rid = self.RGrid.GetCellValue(row, 0)
            cell_val = self.RGrid.GetCellValue(row, col)
            print('Row/Col: {} | {} | {}'.format(rid, clbl, cell_val))
            # update reaction info
            self.UpdateReactionInfo(rid)
            print(cell_val, type(cell_val))
            cell_val = str(cell_val)
            if clbl == 'LB':
                if self.chkFloat(cell_val):
                    print(cell_val, type(cell_val))
                    print('\nLB old', self._cmod_.getReactionLowerBound(rid))
                    self._cmod_.setReactionLowerBound(rid, cell_val)
                    print('LB new', self._cmod_.getReactionLowerBound(rid))
                    if cell_val == 'inf' or cell_val == '-inf':
                        self.writeCmd("cmod.setReactionLowerBound('%s', '%s')" % (rid, cell_val))
                    else:
                        self.writeCmd("cmod.setReactionLowerBound('%s', %s)" % (rid, cell_val))
                else:
                    self.RGrid.SetCellValue(row, col, str(self._cmod_.getReactionLowerBound(rid)))
            elif clbl == 'UB':
                if self.chkFloat(cell_val):
                    print('\nUB old', self._cmod_.getReactionUpperBound(rid))
                    self._cmod_.setReactionUpperBound(rid, cell_val)
                    print('UB new', self._cmod_.getReactionUpperBound(rid))
                    if cell_val == 'inf' or cell_val == '-inf':
                        self.writeCmd("cmod.setReactionUpperBound('%s', '%s')" % (rid, cell_val))
                    else:
                        self.writeCmd("cmod.setReactionUpperBound('%s', %s)" % (rid, cell_val))
                else:
                    self.RGrid.SetCellValue(row, col, str(self._cmod_.getReactionUpperBound(rid)))
            elif clbl == 'Reaction':
                pass
            elif clbl == 'Name':
                print('\nName old', self._cmod_.getReaction(rid).getName())
                self._cmod_.getReaction(rid).setName(cell_val)
                print('\nName new', self._cmod_.getReaction(rid).getName())
                self.writeCmd("cmod.getReaction('%s').setName('%s')" % (rid, cell_val))
            event.Skip()

        #def EVT_BUT_EnableEdit(self,event):
            #print 'EVT_BUT_EnableEdit'
            #for tb in self.TEXTBOXES_Redit:
                #tb.Enable()
            #event.Skip()

        def SemSBML_id(self,event):
            if self.__ActiveReaction != None:
                R = self._cmod_.getReaction(self.__ActiveReaction)
                self.SelectGridRow(R.getId())
                searchString = R.getId()
                self.CallSemanticSBML(searchString)
            else:
                self.Sinfbox.SetPage("<html><body>'<h1>Please select a reaction!</h1>'</body></html>")
            self.NoteB1.ChangeSelection(5)
            event.Skip()

        def SemSBML_name(self,event):
            if self.__ActiveReaction != None:
                R = self._cmod_.getReaction(self.__ActiveReaction)
                self.SelectGridRow(R.getId())
                searchString = R.getName()
                self.CallSemanticSBML(searchString)
            else:
                self.Sinfbox.SetPage("<html><body>'<h1>Please select a reaction!</h1>'</body></html>")
            self.NoteB1.ChangeSelection(5)
            event.Skip()

        def CallSemanticSBML(self, searchString):
            print('MENUAnalyseSemanticSBML')
            #self.TEXTBOX_Annotate.Enable()
            rs = ''
            if HAVE_URLLIB2:
                wait = wx.BusyCursor()
                #self.TEXTBOX_Annotate.Clear()

                rs = '<h2>SemanticSBML query</h2><br/><h3>%s</h3>' % searchString
                # semanticSBML REST web services
                site_root = "www.semanticsbml.org"
                reply_mode = '.xml'
                try:
                    self.SemanticSBMLClient.Connect(site_root)
                    searchString = self.SemanticSBMLClient.URLEncode(searchString)
                    RESTquery = "/semanticSBML/annotate/search%s?q=%s" % (reply_mode, searchString.strip().replace(' ','+'))
                    print(RESTquery)
                    data1 = self.SemanticSBMLClient.Get(RESTquery)
                    data1 = self.SemanticSBMLClient.URLDecode(data1)
                    self.SemanticSBMLClient.Close()
                    self.writeCmd('SemanticSBMLClient = cbm.CBNetDB.SemanticSBML()')
                    self.writeCmd('SemanticSBMLClient.Connect(\"%s\")' % (site_root))
                    self.writeCmd('data1 = SemanticSBMLClient.Get(\"%s\")' % (RESTquery))
                    self.writeCmd('SemanticSBMLClient.Close()')
                    print(self.SemanticSBMLClient.GetLog())

                    item_re = re.compile('<item>.+?</item>')
                    items = re.findall(item_re, data1)
                    items = [i.replace('<item>','').replace('</item>','').strip() for i in items]
                    #self.TEXTBOX_Annotate.write('SemanticSBML results for: %s\n\n' % (urllib2.unquote(searchString)))
                    rs += '<table border="1" cellpadding="5">'
                    for i in items:
                        #url = self.SemanticSBMLClient.MiriamURN2IdentifiersURL(i)
                        url = i
                        rs += '<tr><td>%s</td><td><a href=%s>%s</a></td></tr>' % (i, url, url)
                        print(i)
                        #self.TEXTBOX_Annotate.write('%s\n' % i)<a href="/Education/">NCBI Education</a>
                    rs += '</table border="1" cellpadding="5">'
                    del wait


                except Exception as ex:
                    rs += '<h1>Error connecting to: %s</h1>' % (site_root)
                    #self.TEXTBOX_Annotate.WriteText('\n******************************************************\n')
                    #self.TEXTBOX_Annotate.WriteText('\n* Error connecting to: %s *\n\n' % (site_root))
                    #self.TEXTBOX_Annotate.WriteText('******************************************************\n\n')
                    print(ex)
            elif HAVE_URLLIB2:
                #self.TEXTBOX_Annotate.WriteText('Please select a reaction!')
                rs += '<h1>Please select a reaction!</h1>'
            else:
                #self.TEXTBOX_Annotate.WriteText('HTTPLIB not available')
                rs += '<h1>HTTPLIB not available!</h1>'
            self.Sinfbox.SetPage("<html><body>%s</body></html>" % rs)



        def EVT_BUT_optimise(self, event):
            wait = wx.BusyCursor()
            self._cbm_.CBSolver.analyzeModel(self._cmod_)
            self.UpdateRGridData()
            self.UpdateModelStatus()
            self.writeCmd("cbmpy.CBSolver.analyzeModel(cmod)")
            del wait
            event.Skip()

        def EVT_BUT_minsum(self,event):
            wait = wx.BusyCursor()
            self._cbm_.CBSolver.cplx_MinimizeSumOfAbsFluxes(self._cmod_, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower', optPercentage=100.0, work_dir=None, quiet=False, debug=False, objective_coefficients={}, return_lp_obj=False)
            self.UpdateRGridData()
            self.UpdateModelStatus()
            self.writeCmd("cbmpy.CBSolver.cplx_MinimizeSumOfAbsFluxes(cmod, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower', optPercentage=100.0, objective_coefficients={}, return_lp_obj=False)")
            del wait
            event.Skip()
            ##  cbmpy.CBSolver.cplx_MinimizeSumOfAbsFluxes

        def PrintCoords(self,event):
            print("coords are: {}".format(event.Coords), end=" ")
            print("pixel coords are: {}\n".format(event.GetPosition()), end=" ")

        def EVT_FRAME_resize(self,event):
            #print "Pixel coords are: %s"%(event.GetPosition(),)
            print('Panel size', self.MainPanel.GetSize())
            #print 'Frame size', self.GetSize()
            #print 'Panel minsize', self.MainPanel.GetMinSize()
            #print 'Frame minsize', self.GetMinSize()
            #print 'Panel virtualsize', self.MainPanel.GetVirtualSize()
            #print 'Frame virtualsize', self.GetVirtualSize()
            #print 'mainSizer minsize', self.mainSizer.GetMinSize()
            self.FRAME_resize()
            #print '\n'

        def FRAME_resize(self):
            # get the frame size and adjust the main panel to it

            self.PanelSize = self.GetVirtualSize()
            self.MainPanel.SetSize(self.GetVirtualSize())

            #LPsize = wx.Size(self.PanelSize[0]*self.LPwidth, self.PanelSize[1]*self.LPheight)
            #RPsize = wx.Size(self.PanelSize[0]*self.RPwidth, self.PanelSize[1]*self.RPheight)
            #LPsize = wx.Size(self.PanelSize[0]*self.LPwidth, self.RGridInitSize[1])
            LPsize = wx.Size(self.PanelSize[0]*self.LPwidth, self.PanelSize[1]-self.RGridGap)
            RPsize = wx.Size(self.PanelSize[0]*self.RPwidth, self.PanelSize[1])
            # using the new panels to resize the components
            # self.StatusPanel.SetMinSize(wx.Size(LPsize[0],self.StatusPanel.GetVirtualSize()[1]))
            self.RGrid_scrollwindow.SetSize(LPsize)
            self.RGrid.SetSize(LPsize)
            self.NoteB1.SetSize(RPsize)
            RinfSize = wx.Size(self.PanelSize[0]*self.RPwidth-8, self.PanelSize[1]*self.RPheight-25)

            self.PyBox.SetSize(RinfSize)
            self.Rinfbox.SetSize(RinfSize)
            self.Sinfbox.SetSize(RinfSize)
            self.MetCanvas.SetSize(RinfSize)
            self.GenCanvas.SetSize(RinfSize)
            self.RightPanel.SetSize(RPsize)
            #self.LeftPanel.Layout()
            self.MainPanel.CenterOnParent()


        def EVT_onClick(self,event):
            print("Pixel coords are: {}".format(event.GetPosition()), end=" ")
            #print 'Panel size', self.MainPanel.GetSize()
            #print 'Frame size', self.GetSize()
            #print 'Panel minsize', self.MainPanel.GetMinSize()
            #print 'Frame minsize', self.GetMinSize()
            #print 'Panel virtualsize', self.MainPanel.GetVirtualSize()
            #print 'Frame virtualsize', self.GetVirtualSize()
            #print '\nmainSizer minsize', self.mainSizer.GetMinSize()
            #print '\n'


        def __C_BUFF_ADD__(self, obj):
            self.__C_BUFF__.append(obj)
            if len(self.__C_BUFF__) > 10:
                return self.__C_BUFF__.pop(0)
            else:
                return None


        def HIT_STB_onClick(self, Object):
            print('')
            print(repr(Object))
            print(Object.Name + " got Hit with Left")
            if Object.Name in self.Smap:
                self.UpdateSpeciesGraph(Object.Name)
            if Object.Name in self.Rmap:
                self.UpdateReactionGraph(Object.Name)
                self.updateInfoFromReactionName(Object.Name)
                self.RGrid.SelectRow(self.RGridRow.index(Object.Name))
                self.RGrid.MakeCellVisible(self.RGridRow.index(Object.Name), 0)



            ##  Object.SetLineColor('Green')
            ##  self.FCanvas_met.Draw(True)
            ##  self.FCanvas_met.ZoomToBB()
            ##  wx.CallAfter(self.FCanvas_met.ZoomToBB)
            ##  wx.CallAfter(self.FCanvas_met.Draw)


    class HtmlWindowMod(wx.html.HtmlWindow):
        """
        Overrides 'OnLinkClicked' to open links in external browser

        """
        def __init__(self, *args, **kwargs):
            wx.html.HtmlWindow.__init__(*args, **kwargs)
            if "gtk2" in wx.PlatformInfo:
                self.SetStandardFonts()

        def OnLinkClicked(self, link):
            wx.LaunchDefaultBrowser(link.GetHref())

    def circlePoints(totalPoints=4, startAngle=0, arc=360, circleradius=1, centerxy=(0,0), direction='forward', evenDistribution=True):
        """
        Returns a list of points evenly spread around a circle:

         - *totalPoints* how many points
         - *startAngle* where to start
         - *arc* how far to go
         - *circleradius* radius
         - *centerxy* origin
         - *direction* 'forward' or 'backward'
         - *evenDistribution* True/False

        This code has been adapted from the Flash example that can be found here:
        http://www.lextalkington.com/blog/2009/12/generate-points-around-a-circles-circumference/
        """

        ##  totalPoints = 4
        ##  startAngle = 270
        ##  arc = 180
        ##  direction = 'forward'
        ##  evenDistribution = True

        ##  circleradius = 1
        ##  centerx = 0
        ##  centery = 0

        roundfact = 10
        mpi = math.pi/180.0
        startRadians = startAngle * mpi
        incrementAngle = float(arc)/float(totalPoints)
        incrementRadians = incrementAngle * mpi

        if arc < 360:
            # this spreads the points out evenly across the arc
            if evenDistribution:
                incrementAngle = float(arc)/(float(totalPoints-1))
                incrementRadians = incrementAngle * mpi
            else:
                incrementAngle = float(arc)/float(totalPoints)
                incrementRadians = incrementAngle * mpi;

        cxy = []
        for p in range(totalPoints):
            xp = centerxy[0] + math.sin(startRadians) * circleradius
            yp = centerxy[1] + math.cos(startRadians) * circleradius
            if direction == 'forward':
                startRadians += incrementRadians
            else:
                startRadians -= incrementRadians
            cxy.append((round(xp,roundfact), round(yp, roundfact)))
            ##  print (round(xp,roundfact), round(yp, roundfact))
        return cxy


    '''

    wxPython Custom Widget Collection 20060207
    Written By: Edward Flick (eddy -=at=- cdf-imaging -=dot=- com)
                Michele Petrazzo (michele -=dot=- petrazzo -=at=- unipex -=dot=- it)
                Will Sadkin (wsadkin-=at=- nameconnector -=dot=- com)
    Copyright 2006 (c) CDF Inc. ( http://www.cdf-imaging.com )
    Contributed to the wxPython project under the wxPython project's license.

    '''


    #----------------------------------------------------------------------
    def getSmallUpArrowData():
        return \
    '\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x10\x00\x00\x00\x10\x08\x06\
    \x00\x00\x00\x1f\xf3\xffa\x00\x00\x00\x04sBIT\x08\x08\x08\x08|\x08d\x88\x00\
    \x00\x00<IDAT8\x8dcddbf\xa0\x040Q\xa4{h\x18\xf0\xff\xdf\xdf\xffd\x1b\x00\xd3\
    \x8c\xcf\x10\x9c\x06\xa0k\xc2e\x08m\xc2\x00\x97m\xd8\xc41\x0c \x14h\xe8\xf2\
    \x8c\xa3)q\x10\x18\x00\x00R\xd8#\xec\xb2\xcd\xc1Y\x00\x00\x00\x00IEND\xaeB`\
    \x82'

    def getSmallUpArrowBitmap():
        return BitmapFromImage(getSmallUpArrowImage())

    def getSmallUpArrowImage():
        stream = csio.StringIO(getSmallUpArrowData())
        return ImageFromStream(stream)


    def getSmallDnArrowData():
        return \
    "\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x10\x00\x00\x00\x10\x08\x06\
    \x00\x00\x00\x1f\xf3\xffa\x00\x00\x00\x04sBIT\x08\x08\x08\x08|\x08d\x88\x00\
    \x00\x00HIDAT8\x8dcddbf\xa0\x040Q\xa4{\xd4\x00\x06\x06\x06\x06\x06\x16t\x81\
    \xff\xff\xfe\xfe'\xa4\x89\x91\x89\x99\x11\xa7\x0b\x90%\ti\xc6j\x00>C\xb0\x89\
    \xd3.\x10\xd1m\xc3\xe5*\xbc.\x80i\xc2\x17.\x8c\xa3y\x81\x01\x00\xa1\x0e\x04e\
    ?\x84B\xef\x00\x00\x00\x00IEND\xaeB`\x82"

    def getSmallDnArrowBitmap():
        return BitmapFromImage(getSmallDnArrowImage())

    def getSmallDnArrowImage():
        stream = csio.StringIO(getSmallDnArrowData())
        return ImageFromStream(stream)
    #----------------------------------------------------------------------

    class myListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin):
        def __init__(self, parent, ID=-1, pos=wx.DefaultPosition,
                     size=wx.DefaultSize, style=0):
            wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
            listmix.ListCtrlAutoWidthMixin.__init__(self)

    class TextCtrlAutoComplete (wx.TextCtrl, listmix.ColumnSorterMixin ):

        def __init__ ( self, parent, colNames=None, choices = None,
                       multiChoices=None, showHead=True, dropDownClick=True,
                       colFetch=-1, colSearch=0, hideOnNoMatch=True,
                       selectCallback=None, entryCallback=None, matchFunction=None,
                       **therest) :
            '''
            Constructor works just like wx.TextCtrl except you can pass in a
            list of choices.  You can also change the choice list at any time
            by calling setChoices.
            '''

            if 'style' in therest:
                therest['style']=wx.TE_PROCESS_ENTER | therest['style']
            else:
                therest['style']=wx.TE_PROCESS_ENTER

            wx.TextCtrl.__init__(self, parent, **therest )

            #Some variables
            self._dropDownClick = dropDownClick
            self._colNames = colNames
            self._multiChoices = multiChoices
            self._showHead = showHead
            self._choices = choices
            self._lastinsertionpoint = 0
            self._hideOnNoMatch = hideOnNoMatch
            self._selectCallback = selectCallback
            self._entryCallback = entryCallback
            self._matchFunction = matchFunction

            self._screenheight = wx.SystemSettings.GetMetric( wx.SYS_SCREEN_Y )

            #sort variable needed by listmix
            self.itemDataMap = dict()

            #Load and sort data
            if not (self._multiChoices or self._choices):
                raise ValueError("Pass me at least one of multiChoices OR choices")

            #widgets
            self.dropdown = wx.PopupWindow( self )

            #Control the style
            flags = wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.LC_SORT_ASCENDING
            if not (showHead and multiChoices) :
                flags = flags | wx.LC_NO_HEADER

            #Create the list and bind the events
            self.dropdownlistbox = myListCtrl( self.dropdown, style=flags,
                                               pos=wx.Point( 0, 0) )

            #initialize the parent
            if multiChoices: ln = len(multiChoices)
            else: ln = 1
            #else: ln = len(choices)
            listmix.ColumnSorterMixin.__init__(self, ln)

            #load the data
            if multiChoices: self.SetMultipleChoices (multiChoices, colSearch=colSearch, colFetch=colFetch)
            else: self.SetChoices ( choices )

            gp = self
            while ( gp != None ) :
                gp.Bind ( wx.EVT_MOVE , self.onControlChanged, gp )
                gp.Bind ( wx.EVT_SIZE , self.onControlChanged, gp )
                gp = gp.GetParent()

            self.Bind( wx.EVT_KILL_FOCUS, self.onControlChanged, self )
            self.Bind( wx.EVT_TEXT , self.onEnteredText, self )
            self.Bind( wx.EVT_KEY_DOWN , self.onKeyDown, self )

            #If need drop down on left click
            if dropDownClick:
                self.Bind ( wx.EVT_LEFT_DOWN , self.onClickToggleDown, self )
                self.Bind ( wx.EVT_LEFT_UP , self.onClickToggleUp, self )

            self.dropdown.Bind( wx.EVT_LISTBOX , self.onListItemSelected, self.dropdownlistbox )
            self.dropdownlistbox.Bind(wx.EVT_LEFT_DOWN, self.onListClick)
            self.dropdownlistbox.Bind(wx.EVT_LEFT_DCLICK, self.onListDClick)
            self.dropdownlistbox.Bind(wx.EVT_LIST_COL_CLICK, self.onListColClick)

            #self.il = wx.ImageList(16, 16)

            #self.sm_dn = self.il.Add(getSmallDnArrowBitmap())
            #self.sm_up = self.il.Add(getSmallUpArrowBitmap())

            #self.dropdownlistbox.SetImageList(self.il, wx.IMAGE_LIST_SMALL)
            #self._ascending = True


        #-- methods called from mixin class
        def GetSortImages(self):
            return (self.sm_dn, self.sm_up)

        def GetListCtrl(self):
            return self.dropdownlistbox

        # -- event methods
        def onListClick(self, evt):
            toSel, flag = self.dropdownlistbox.HitTest( evt.GetPosition() )
            #no values on poition, return
            if toSel == -1: return
            self.dropdownlistbox.Select(toSel)

        def onListDClick(self, evt):
            self._setValueFromSelected()

        def onListColClick(self, evt):
            col = evt.GetColumn()

            #reverse the sort
            if col == self._colSearch:
                self._ascending = not self._ascending

            self.SortListItems( evt.GetColumn(), ascending=self._ascending )
            self._colSearch = evt.GetColumn()
            evt.Skip()

        def onEnteredText(self, event):
            text = event.GetString()

            if self._entryCallback:
                self._entryCallback()

            if not text:
                # control is empty; hide dropdown if shown:
                if self.dropdown.IsShown():
                    self._showDropDown(False)
                event.Skip()
                return


            found = False
            if self._multiChoices:
                #load the sorted data into the listbox
                dd = self.dropdownlistbox
                choices = [dd.GetItem(x, self._colSearch).GetText()
                           for x in xrange(dd.GetItemCount())]
            else:
                choices = self._choices

            for numCh, choice in enumerate(choices):
                if self._matchFunction and self._matchFunction(text, choice):
                    found = True
                elif choice.lower().startswith(text.lower()) :
                    found = True
                if found:
                    self._showDropDown(True)
                    item = self.dropdownlistbox.GetItem(numCh)
                    toSel = item.GetId()
                    self.dropdownlistbox.Select(toSel)
                    break

            if not found:
                self.dropdownlistbox.Select(self.dropdownlistbox.GetFirstSelected(), False)
                if self._hideOnNoMatch:
                    self._showDropDown(False)

            self._listItemVisible()

            event.Skip ()

        def onKeyDown ( self, event ) :
            """ Do some work when the user press on the keys:
                up and down: move the cursor
                left and right: move the search
            """
            skip = True
            sel = self.dropdownlistbox.GetFirstSelected()
            visible = self.dropdown.IsShown()

            KC = event.GetKeyCode()
            if KC == wx.WXK_DOWN :
                if sel < (self.dropdownlistbox.GetItemCount () - 1) :
                    self.dropdownlistbox.Select ( sel+1 )
                    self._listItemVisible()
                self._showDropDown ()
                skip = False
            elif KC == wx.WXK_UP :
                if sel > 0 :
                    self.dropdownlistbox.Select ( sel - 1 )
                    self._listItemVisible()
                self._showDropDown ()
                skip = False
            elif KC == wx.WXK_LEFT :
                if not self._multiChoices: return
                if self._colSearch > 0:
                    self._colSearch -=1
                self._showDropDown ()
            elif KC == wx.WXK_RIGHT:
                if not self._multiChoices: return
                if self._colSearch < self.dropdownlistbox.GetColumnCount() -1:
                    self._colSearch += 1
                self._showDropDown()

            if visible :
                if event.GetKeyCode() == wx.WXK_RETURN :
                    self._setValueFromSelected()
                    skip = False
                if event.GetKeyCode() == wx.WXK_ESCAPE :
                    self._showDropDown( False )
                    skip = False
            if skip :
                event.Skip()

        def onListItemSelected (self, event):
            self._setValueFromSelected()
            event.Skip()

        def onClickToggleDown(self, event):
            self._lastinsertionpoint = self.GetInsertionPoint()
            event.Skip ()

        def onClickToggleUp ( self, event ) :
            if ( self.GetInsertionPoint() == self._lastinsertionpoint ) :
                self._showDropDown ( not self.dropdown.IsShown() )
            event.Skip ()

        def onControlChanged(self, event):
            self._showDropDown( False )
            event.Skip()


        # -- Interfaces methods
        def SetMultipleChoices(self, choices, colSearch=0, colFetch=-1):
            ''' Set multi-column choice
            '''
            self._multiChoices = choices
            self._choices = None
            if not isinstance(self._multiChoices, list):
                self._multiChoices = [ x for x in self._multiChoices]

            flags = wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.LC_SORT_ASCENDING
            if not self._showHead:
                flags |= wx.LC_NO_HEADER
            self.dropdownlistbox.SetWindowStyleFlag(flags)

            #prevent errors on "old" systems
            if sys.version.startswith("2.3"):
                self._multiChoices.sort(lambda x, y: cmp(x[0].lower(), y[0].lower()))
            else:
                self._multiChoices.sort(key=lambda x: locale.strxfrm(x[0]).lower() )

            self._updateDataList(self._multiChoices)

            lChoices = len(choices)
            if lChoices < 2:
                raise ValuError( "You have to pass me a multi-dimension list")

            for numCol, rowValues in enumerate(choices[0]):

                if self._colNames: colName = self._colNames[numCol]
                else: colName = "Select %i" % numCol

                self.dropdownlistbox.InsertColumn(numCol, colName)

            for numRow, valRow in enumerate(choices):

                for numCol, colVal in enumerate(valRow):
                    if numCol == 0:
                        index = self.dropdownlistbox.InsertImageStringItem(sys.maxint, colVal, -1)
                    self.dropdownlistbox.SetStringItem(index, numCol, colVal)
                    self.dropdownlistbox.SetItemData(index, numRow)

            self._setListSize()
            self._colSearch = colSearch
            self._colFetch = colFetch

        def SetChoices(self, choices):
            '''
            Sets the choices available in the popup wx.ListBox.
            The items will be sorted case insensitively.
            '''
            self._choices = choices
            self._multiChoices = None
            flags = wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.LC_SORT_ASCENDING | wx.LC_NO_HEADER
            self.dropdownlistbox.SetWindowStyleFlag(flags)

            if not isinstance(choices, list):
                self._choices = [ x for x in choices]

            #prevent errors on "old" systems
            if sys.version.startswith("2.3"):
                self._choices.sort(lambda x, y: cmp(x.lower(), y.lower()))
            else:
                self._choices.sort(key=lambda x: locale.strxfrm(x).lower())

            self._updateDataList(self._choices)

            self.dropdownlistbox.InsertColumn(0, "")

            for num, colVal in enumerate(self._choices):
                index = self.dropdownlistbox.InsertImageStringItem(sys.maxint, colVal, -1)

                self.dropdownlistbox.SetStringItem(index, 0, colVal)
                self.dropdownlistbox.SetItemData(index, num)

            self._setListSize()

            # there is only one choice for both search and fetch if setting a single column:
            self._colSearch = 0
            self._colFetch = -1

        def GetChoices(self):
            if self._choices:
                return self._choices
            else:
                return self._multiChoices

        def SetSelectCallback(self, cb=None):
            self._selectCallback = cb

        def SetEntryCallback(self, cb=None):
            self._entryCallback = cb

        def SetMatchFunction(self, mf=None):
            self._matchFunction = mf


        #-- Internal methods
        def _setValueFromSelected( self ) :
            '''
            Sets the wx.TextCtrl value from the selected wx.ListCtrl item.
            Will do nothing if no item is selected in the wx.ListCtrl.
            '''
            sel = self.dropdownlistbox.GetFirstSelected()
            if sel > -1:
                if self._colFetch != -1: col = self._colFetch
                else: col = self._colSearch

                itemtext = self.dropdownlistbox.GetItem(sel, col).GetText()
                if self._selectCallback:
                    dd = self.dropdownlistbox
                    values = [dd.GetItem(sel, x).GetText()
                              for x in xrange(dd.GetColumnCount())]
                    self._selectCallback( values )

                self.SetValue (itemtext)
                self.SetInsertionPointEnd ()
                self.SetSelection ( -1, -1 )
                self._showDropDown ( False )


        def _showDropDown ( self, show = True ) :
            '''
            Either display the drop down list (show = True) or hide it (show = False).
            '''
            if show :
                size = self.dropdown.GetSize()
                width, height = self . GetSizeTuple()
                x, y = self . ClientToScreenXY ( 0, height )
                if size.GetWidth() != width :
                    size.SetWidth(width)
                    self.dropdown.SetSize(size)
                    self.dropdownlistbox.SetSize(self.dropdown.GetClientSize())
                if (y + size.GetHeight()) < self._screenheight :
                    self.dropdown . SetPosition ( wx.Point(x, y) )
                else:
                    self.dropdown . SetPosition ( wx.Point(x, y - height - size.GetHeight()) )
            self.dropdown.Show ( show )

        def _listItemVisible( self ) :
            '''
            Moves the selected item to the top of the list ensuring it is always visible.
            '''
            toSel =  self.dropdownlistbox.GetFirstSelected ()
            if toSel == -1: return
            self.dropdownlistbox.EnsureVisible( toSel )

        def _updateDataList(self, choices):
            #delete, if need, all the previous data
            if self.dropdownlistbox.GetColumnCount() != 0:
                self.dropdownlistbox.DeleteAllColumns()
                self.dropdownlistbox.DeleteAllItems()

            #and update the dict
            if choices:
                for numVal, data in enumerate(choices):
                    self.itemDataMap[numVal] = data
            else:
                numVal = 0
            self.SetColumnCount(numVal)

        def _setListSize(self):
            if self._multiChoices:
                choices = self._multiChoices
            else:
                choices = self._choices

            longest = 0
            for choice in choices :
                longest = max(len(choice), longest)

            longest += 6
            itemcount = min( len( choices ) , 20 ) + 2
            charheight = self.dropdownlistbox.GetCharHeight()
            charwidth = self.dropdownlistbox.GetCharWidth()
            self.popupsize = wx.Size( charwidth*longest, charheight*itemcount )
            self.dropdownlistbox.SetSize ( self.popupsize )
            self.dropdown.SetClientSize( self.popupsize )

    class MyAUIFrame(wx.Frame):
        def __init__(self, *args, **kwargs):
            wx.Frame.__init__(self, *args, **kwargs)

            self.mgr = wx.aui.AuiManager(self)

            leftpanel = wx.Panel(self, -1, size = (200, 150))
            rightpanel = wx.Panel(self, -1, size = (200, 150))
            bottompanel = wx.Panel(self, -1, size = (200, 150))
            toppanel = wx.Panel(self, -1, size = (200, 150))
            centerpanel = wx.Panel(self, -1, size = (200, 150))

            wx.TextCtrl(rightpanel, -1, 'rightpanel')
            wx.TextCtrl(leftpanel, -1, 'leftpanel')
            wx.TextCtrl(toppanel, -1, 'toppanel')
            wx.TextCtrl(bottompanel, -1, 'bottompanel')
            wx.TextCtrl(centerpanel, -1, 'centerpanel')

            self.mgr.AddPane(leftpanel, wx.aui.AuiPaneInfo().Left().Layer(1))
            self.mgr.AddPane(rightpanel, wx.aui.AuiPaneInfo().Right().Layer(1))
            self.mgr.AddPane(bottompanel, wx.aui.AuiPaneInfo().Bottom().Layer(2))
            self.mgr.AddPane(toppanel, wx.aui.AuiPaneInfo().Top().Layer(2))
            self.mgr.AddPane(centerpanel, wx.aui.AuiPaneInfo().Center().Layer(1))
            self.Maximize()

            self.mgr.Update()

    class MyAUIApp(wx.App):
        def OnInit(self):
            frame = MyAUIFrame(None, -1, 'CBMPy Gen2 GUI')
            frame.Show()
            self.SetTopWindow(frame)
            return 1

    def runModelEditor(mod):
        app = wx.App(False)
        frame = ModelEditor(mod)
        frame.Show(True)
        app.SetTopWindow(frame)
        app.MainLoop()

    def runMyAUIApp():
        app = MyAUIApp(0)
        app.MainLoop()
