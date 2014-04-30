from PyQt4 import uic
from os import path
import sys
from PyQt4.QtGui import QMainWindow, QWidget, QSizePolicy, QMenu, QInputDialog
from PyQt4 import QtCore
from matplotlib.figure import Figure
from PyQt4.QtCore import QModelIndex
with open('MainWindow.ui', 'r') as infile:
    with open('ui_MainWindow.py', 'w') as outfile:
        uic.compileUi(infile, outfile)
import ui_MainWindow

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from DataModel.TreeModel import Node, TreeModel
from DataModel.Database import RMFClasses, RMFParams

class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = ui_MainWindow.Ui_MainWindow()
        self.ui.setupUi(self)
#       self.main_frame = QWidget()
#         
#         # Create the mpl Figure and FigCanvas objects. 
#         # 5x4 inches, 100 dots-per-inch
#         #
#         self.dpi = 100
#         self.fig = Figure((5.0, 4.0), dpi=self.dpi)
#         self.canvas = FigureCanvas(self.fig)
#         self.canvas.setParent(self.main_frame)
#         self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
#         self.canvas.updateGeometry()
#         
#         # Since we have only one plot, we can use add_axes 
#         # instead of add_subplot, but then the subplot
#         # configuration tool in the navigation toolbar wouldn't
#         # work.
#         #
#         self.axes = self.fig.add_subplot(111)
#         self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
#         self.ui.vlay.addWidget(self.main_frame)
        
        self.database = RMFParams()
        self.classbase = RMFClasses()
#         class1 = RMFClasses.create(class__='TEST')
#         record1 = RMFParams.create(name='Test 1', massfile='', tabfile='',
#                                    constants='', type=class1, class__='TEST', module='')
        
        self.currentWr = None
        
        self.rootNode = Node("root")
        self.model = TreeModel(self.rootNode, self.classbase, self.database)
        self.ui.treeView.setModel(self.model)
        self.ui.treeView.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.treeView.customContextMenuRequested.connect(self.showTreeMenu)        
        self.ui.treeView.selectionModel().selectionChanged.connect(self.onSelectionChanged)
        
        
        
    def onSelectionChanged(self, selected, deselected):
        index = selected.indexes()[0]
        if index.internalPointer().parent() != self.model._rootNode:
            print "It is not a Class"
        
        self.currentWr = self.getWrapper(index)
        
        
    def getWrapper(self, index=QModelIndex()):
        class__=RMFParams.get()
        
    def showTreeMenu(self, pos):
        index = self.ui.treeView.indexAt(pos)
        
#         if not index.isValid():
#             return
        
        self.treeMenu = QMenu(self)
        self.treeMenu.addAction("Add Class", self.treeMenuAddClass)
        if index.isValid():
            self.treeMenu.addAction("Add Params", lambda:self.treeMenuAddParams(index))
        self.treeMenu.exec_(self.mapToGlobal(pos))
        
    def treeMenuAddClass(self):
        print 'AddClass'
        text, ok = QInputDialog.getText(self, 'Input', 'Enter class name:')
        if ok:
            self.model.insertClass(text, 0)
        self.ui.treeView.expandAll()
            
        
    def treeMenuAddParams(self, parentIndex):
        print 'AddParams'
        text, ok = QInputDialog.getText(self, 'Input', 'Enter params name:')
        if ok:
            self.model.insertParams(text, parentIndex)
        self.ui.treeView.expandAll()
        
       
       
        
        